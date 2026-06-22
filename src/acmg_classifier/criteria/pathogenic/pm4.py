"""PM4 -- protein length change due to in-frame indel or stop-loss."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


def _load_pm4_columns(tsv_path) -> tuple[frozenset[str], dict[str, int]]:
    """Load PM4 per-gene config from ``disease_prevalence.tsv``:

    * the set of genes whose VCEP declined PM4 (``pm4`` == ``not_applicable``);
    * ``pm4_supporting_max_aa`` — the in-frame-indel size (amino acids) at or
      below which PM4 downgrades from the default Moderate to Supporting.

    A missing file/column degrades to "no gene declined" / "no size downgrade"."""
    import csv
    declined: set[str] = set()
    max_aa: dict[str, int] = {}
    try:
        if not tsv_path.exists():
            return frozenset(), {}
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                if (row.get("pm4") or "").strip().lower() == "not_applicable":
                    declined.add(gene)
                raw = (row.get("pm4_supporting_max_aa") or "").strip()
                if raw:
                    try:
                        max_aa[gene] = int(raw)
                    except ValueError:
                        pass
    except OSError:
        return frozenset(), {}
    return frozenset(declined), max_aa


class PM4Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._not_applicable, self._supporting_max_aa = _load_pm4_columns(
            cfg.disease_prevalence_tsv
        )
        from acmg_classifier.criteria.pm4_regions import PM4Regions
        self._regions = PM4Regions(cfg.pm4_regions_tsv)
        from acmg_classifier.local_db.conservation import PhyloPReader
        self._phylop = PhyloPReader(cfg.phylop_bigwig)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.PM4, "No consequence")

        # VCEP gate: cancer panels (BRCA1/2, MMR, TP53, APC, PALB2) and the
        # PI3K-pathway specs decline PM4 — withhold it for those genes.
        if pc.gene_symbol in self._not_applicable:
            return CriteriaResult.not_met(
                ACMGCriterion.PM4, f"{pc.gene_symbol}: VCEP designates PM4 not applicable"
            )

        # PM4 fires for changes that alter protein length without truncation
        # (frameshift/stop-gain are PVS1). For in-frame indels inside a
        # repetitive region, the length change is expected/tolerated and the
        # variant is better captured by BP3 — explicitly redirect rather than
        # silently failing so the audit trail reflects the reasoning.
        gene = pc.gene_symbol
        is_indel = pc.consequence in (
            ConsequenceType.INFRAME_INSERTION, ConsequenceType.INFRAME_DELETION,
        )
        is_stoploss = pc.consequence == ConsequenceType.STOP_LOST

        # Nucleotide-conservation PM4 (ABCA4): a synonymous or missense variant
        # affecting >1 nucleotide at a highly-conserved position (phyloP >= cutoff)
        # earns PM4_Moderate. A single-nucleotide change (SNV) does NOT qualify —
        # the ABCA4 VCEP restricts this rule to multi-nucleotide events. Skipped
        # when phyloP is unavailable.
        nt_cut = self._regions.nt_phylop(gene)
        if nt_cut is not None and pc.consequence in (
            ConsequenceType.SYNONYMOUS, ConsequenceType.MISSENSE,
        ):
            best = _max_phylop(self._phylop, variant)
            n_nt = _changed_nt(variant)
            if best is not None and best >= nt_cut and n_nt > 1:
                return CriteriaResult.met(
                    ACMGCriterion.PM4, strength=CriterionStrength.MODERATE,
                    evidence=(f"{gene}: {pc.consequence.value} at conserved nucleotide "
                              f"(phyloP {best:.2f} >= {nt_cut:g}, {n_nt} nt)"),
                )

        if not (is_indel or is_stoploss):
            return CriteriaResult.not_met(ACMGCriterion.PM4, "Not an in-frame indel or stop-loss")

        if annotation.repeat and annotation.repeat.in_repeat:
            return CriteriaResult.not_met(
                ACMGCriterion.PM4,
                f"In-frame indel in repeat ({annotation.repeat.repeat_class}); use BP3",
            )

        # Stop-loss: a VCEP may set a stop-loss-specific strength (or decline it —
        # CYP1B1, where stop-loss is not a disease mechanism). Stop-loss is not
        # size-scoped, so it keeps the default Moderate when no override exists.
        if is_stoploss:
            sl = self._regions.stoploss_strength(gene)
            if sl == "not_applicable":
                return CriteriaResult.not_met(
                    ACMGCriterion.PM4, f"{gene}: VCEP PM4 not applicable to stop-loss"
                )
            if isinstance(sl, CriterionStrength):
                return CriteriaResult.met(
                    ACMGCriterion.PM4, strength=sl,
                    evidence=f"stop_lost ({gene}: VCEP PM4_{sl.value})",
                )
            return CriteriaResult.met(
                ACMGCriterion.PM4, evidence="stop_lost outside repeat region",
            )

        # In-frame indel. Size-based Supporting (single / <3 aa) is eligible when
        # the VCEP set a size cutoff (pm4_supporting_max_aa).
        max_aa = self._supporting_max_aa.get(gene)
        aa = _indel_aa_change(variant)
        size_eligible = max_aa is not None and aa is not None and aa <= max_aa

        # Gene-specific PM4 region rule (Strong residues, allow/deny regions,
        # region default, conservation / deletion-content gates), if any.
        if self._regions.has_gene(gene):
            # Conservation gate (RPE65/CTLA4/PIK3R1): PM4 only when the indel sits
            # at a conserved position (PhyloP > cutoff). Skipped (PM4 proceeds)
            # when phyloP is unavailable — graceful degradation, like BP7.
            cutoff = self._regions.conserved_phylop(gene)
            if cutoff is not None and _not_conserved(self._phylop, variant, cutoff):
                return CriteriaResult.not_met(
                    ACMGCriterion.PM4,
                    f"{gene}: in-frame indel not at a conserved position "
                    f"(phyloP <= {cutoff})",
                )
            # Deletion-content gate (SCID panel): a DELETION earns PM4 only if its
            # deleted genomic span contains a known ClinVar P/LP (Moderate) or VUS
            # (Supporting) variant.
            if (self._regions.requires_deletion_content(gene)
                    and pc.consequence == ConsequenceType.INFRAME_DELETION):
                from acmg_classifier.local_db.clinvar_sqlite import (
                    query_pm4_deletion_content,
                )
                end = variant.pos + max(0, len(variant.ref) - 1)
                content = query_pm4_deletion_content(
                    self._cfg.clinvar_sqlite, gene, variant.chrom, variant.pos, end,
                )
                if content == "pathogenic":
                    return CriteriaResult.met(
                        ACMGCriterion.PM4, strength=CriterionStrength.MODERATE,
                        evidence=f"{gene}: deleted region contains a ClinVar P/LP variant",
                    )
                if content == "vus":
                    return CriteriaResult.met(
                        ACMGCriterion.PM4, strength=CriterionStrength.SUPPORTING,
                        evidence=f"{gene}: deleted region contains a ClinVar VUS",
                    )
                return CriteriaResult.not_met(
                    ACMGCriterion.PM4,
                    f"{gene}: deleted region contains no ClinVar P/LP or VUS variant",
                )
            rstr = self._regions.indel_strength(gene, pc.protein_position)
            if rstr == "not_met":
                # A small indel still earns Supporting even where Moderate is
                # withheld (the size-Supporting tier is region-independent).
                if size_eligible:
                    return CriteriaResult.met(
                        ACMGCriterion.PM4, strength=CriterionStrength.SUPPORTING,
                        evidence=f"in-frame indel of <={max_aa} aa ({gene}: PM4_Supporting)",
                    )
                return CriteriaResult.not_met(
                    ACMGCriterion.PM4,
                    f"{gene}: in-frame indel outside PM4 region / in denied region",
                )
            if isinstance(rstr, CriterionStrength):
                if rstr == CriterionStrength.MODERATE and size_eligible:
                    return CriteriaResult.met(
                        ACMGCriterion.PM4, strength=CriterionStrength.SUPPORTING,
                        evidence=f"in-frame indel of <={max_aa} aa ({gene}: PM4_Supporting)",
                    )
                return CriteriaResult.met(
                    ACMGCriterion.PM4, strength=rstr,
                    evidence=f"in-frame indel ({gene}: VCEP PM4_{rstr.value} region)",
                )
            # rstr is None → no region default → fall through to the flat default.

        if size_eligible:
            return CriteriaResult.met(
                ACMGCriterion.PM4, strength=CriterionStrength.SUPPORTING,
                evidence=f"in-frame indel of <={max_aa} aa ({gene}: VCEP PM4_Supporting)",
            )
        return CriteriaResult.met(
            ACMGCriterion.PM4, evidence="in-frame indel outside repeat region",
        )


def _max_phylop(phylop, variant: VariantRecord) -> float | None:
    """Max phyloP over the variant's reference span, or None when phyloP is
    unavailable / unscored at every position."""
    if not phylop.is_available():
        return None
    best: float | None = None
    for offset in range(max(1, len(variant.ref))):
        score = phylop.value(variant.chrom, variant.pos + offset)
        if score is not None and (best is None or score > best):
            best = score
    return best


def _changed_nt(variant: VariantRecord) -> int:
    """Number of substituted nucleotides (1 for an SNV). For an equal-length
    substitution it is the count of differing positions; otherwise the longer of
    ref/alt (a conservative upper bound)."""
    ref, alt = variant.ref, variant.alt
    if len(ref) == len(alt):
        return sum(1 for r, a in zip(ref, alt) if r != a) or 1
    return max(len(ref), len(alt))


def _not_conserved(phylop, variant: VariantRecord, cutoff: float) -> bool:
    """True when the variant's span is confidently NOT conserved (max phyloP over
    the ref span <= cutoff), so a conservation-gated PM4 must be withheld. False
    when conserved OR when phyloP is unavailable (gate skipped → PM4 proceeds)."""
    best = _max_phylop(phylop, variant)
    if best is None:
        return False  # phyloP unavailable / unscored → cannot disprove conservation
    return best <= cutoff


def _indel_aa_change(variant: VariantRecord) -> int | None:
    """Net amino-acid length change of an in-frame indel = |len(ref)-len(alt)|/3.
    Returns ``None`` when the length difference is not a whole number of codons
    (defensive — an in-frame indel should always be a multiple of 3)."""
    diff = abs(len(variant.ref) - len(variant.alt))
    if diff == 0 or diff % 3 != 0:
        return None
    return diff // 3
