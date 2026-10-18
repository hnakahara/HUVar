"""PS1 — same amino acid change as established pathogenic variant."""
from __future__ import annotations
import re
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


SPLICE_TAG = "[PS1 splice T2:"
_CANON_HGVS = re.compile(r"c\.[-*]?\d+([+-])([12])(?![0-9])[ACGT]?")


def canonical_offset(pc) -> tuple[str, int] | None:
    """``(sign, offset)`` of a variant at a canonical ±1,2 splice position
    (from the MANE HGVS c. name), or ``None``."""
    m = _CANON_HGVS.search((pc.hgvs_c or "").split(":")[-1])
    if m:
        return m.group(1), int(m.group(2))
    return None


def _no_protein_change(h) -> bool:
    p = (h.hgvs_p or "").split(":")[-1]
    return not p or p.endswith("=") or p.endswith("=)")


_RANK = {
    CriterionStrength.SUPPORTING: 1,
    CriterionStrength.MODERATE: 2,
    CriterionStrength.STRONG: 3,
}


def _ps1_strength(hits) -> CriterionStrength:
    """PS1 strength from the comparator's ClinVar classification (ClinGen SVI).

    Strong when at least one same-change comparator is classified Pathogenic
    (or the P-containing aggregate 'Pathogenic/Likely pathogenic'); Moderate
    when every comparator is only Likely pathogenic. Without this tiering PS1
    always fired at its Strong default, over-weighting LP-only comparators.
    """
    has_pathogenic = any(
        "pathogenic" in (h.clinical_significance or "").lower()
        and (h.clinical_significance or "").lower() != "likely pathogenic"
        for h in hits
    )
    return CriterionStrength.STRONG if has_pathogenic else CriterionStrength.MODERATE


def _cap(strength: CriterionStrength, ceiling: CriterionStrength | None) -> CriterionStrength:
    """Clamp a PS1 strength to the gene's VCEP ceiling (e.g. RMRP → Supporting)."""
    if ceiling is None or _RANK[strength] <= _RANK[ceiling]:
        return strength
    return ceiling


def _aa_ref_alt(amino_acids: str | None) -> tuple[str, str] | None:
    """The 1-letter (ref, alt) from a VEP ``amino_acids`` "X/Y" missense field,
    or None when it is absent or not a single-residue substitution."""
    if not amino_acids or "/" not in amino_acids:
        return None
    ref, _, alt = amino_acids.partition("/")
    ref, alt = ref.strip(), alt.strip()
    if len(ref) == 1 and len(alt) == 1:
        return ref, alt
    return None


class PS1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Per-gene PS1 spec: some VCEPs restrict the splice extension to
        # non-canonical splice nucleotides (InSiGHT MMR genes).
        from acmg_classifier.criteria.ps1_genes import PS1Spec
        self._spec = PS1Spec(cfg.disease_prevalence_tsv)
        # SCN paralogue amino-acid alignment (PS1 analogous-residue route).
        from acmg_classifier.criteria.ps1_paralog import PS1ParalogMap
        self._paralog_map = PS1ParalogMap(cfg.ps1_paralog_map_tsv)

    @property
    def _min_stars(self) -> int:
        v = getattr(self._cfg, "ps1_min_stars", 2)
        return v if isinstance(v, int) and v >= 0 else 2

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # PS1 has two routes:
        #  * Missense -> the amino-acid rule: a DIFFERENT nucleotide change
        #    producing the SAME amino-acid substitution as a known pathogenic
        #    variant (matched on hgvs_p).
        #  * Splice / intronic -> the ClinGen SVI splicing extension: these have
        #    no protein change, so the amino-acid rule can never fire. Instead a
        #    DIFFERENT nucleotide change at the SAME splice-site position is
        #    recognised as having the same predicted splicing effect (matched on
        #    genomic position).
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.PS1, "No primary consequence")

        # VCEP gate: a VCEP that declined PS1 for the gene withholds it (CDH1).
        if self._spec.is_not_applicable(pc.gene_symbol):
            return CriteriaResult.not_met(
                ACMGCriterion.PS1, f"{pc.gene_symbol}: VCEP designates PS1 not applicable"
            )

        if pc.consequence == ConsequenceType.MISSENSE:
            return self._evaluate_missense(variant, pc)
        if pc.consequence in (
            ConsequenceType.SPLICE_ACCEPTOR,
            ConsequenceType.SPLICE_DONOR,
            ConsequenceType.SPLICE_REGION,
            ConsequenceType.INTRON,
            ConsequenceType.SYNONYMOUS,
        ):
            return self._evaluate_splice(variant, pc)
        return CriteriaResult.not_met(ACMGCriterion.PS1, "Not a missense or splice variant")

    def _evaluate_missense(self, variant: VariantRecord, pc) -> CriteriaResult:
        # Same amino-acid change at the same codon via a DIFFERENT nucleotide.
        # Exclude the variant itself by genomic coordinate so it cannot "match
        # itself". Config.ps1_min_stars (default 2): the prior assertion must be a multi-submitter or expert-reviewed record.
        from acmg_classifier.local_db.clinvar_sqlite import query_same_aa_change
        hits = query_same_aa_change(
            self._cfg.clinvar_sqlite,
            pc.gene_symbol,
            pc.hgvs_p,
            exclude_chrom=variant.chrom,
            exclude_pos=variant.pos,
            exclude_ref=variant.ref,
            exclude_alt=variant.alt,
            min_stars=self._min_stars,
        )
        from acmg_classifier.criteria.vcep_comparators import filter_comparators
        hits, vcep_note = filter_comparators(
            pc.gene_symbol, "PS1", hits, variant.assembly.value,
        )
        if hits:
            strength = _cap(_ps1_strength(hits), self._spec.max_strength(pc.gene_symbol))
            evidence = (
                f"ClinVar same AA (different nucleotide), comparator {strength.value}: "
                f"{', '.join(h.variation_id or '' for h in hits[:3])}"
                + (f" [{vcep_note}]" if vcep_note else "")
            )
            return CriteriaResult.met(ACMGCriterion.PS1, strength, evidence)

        # Paralogue (analogous-residue) route: a VCEP may grant PS1 from the same
        # amino-acid change at the ANALOGOUS (same-numbered) residue of a sibling
        # gene — RASopathy groups (HRAS/NRAS/KRAS, MAP2K1/MAP2K2, SOS1/SOS2) and
        # HBA2/HBA1. Only consulted when the same-gene rule did not fire.
        siblings = self._spec.paralog_group(pc.gene_symbol)
        if siblings:
            # No self-exclusion / codon-proximity guard for paralogues: the hit is
            # in a DIFFERENT gene (different chromosome), so it can never be the
            # variant itself and the same-codon proximity rule does not apply.
            para = []
            for sib in siblings:
                para.extend(query_same_aa_change(
                    self._cfg.clinvar_sqlite, sib, pc.hgvs_p, min_stars=self._min_stars,
                ))
            if para:
                fixed = self._spec.paralog_strength(pc.gene_symbol)
                strength = fixed if fixed is not None else _cap(
                    _ps1_strength(para), self._spec.max_strength(pc.gene_symbol)
                )
                evidence = (
                    f"ClinVar same AA in paralogue gene ({', '.join(siblings)}), "
                    f"{strength.value}: {', '.join(h.variation_id or '' for h in para[:3])}"
                )
                return CriteriaResult.met(ACMGCriterion.PS1, strength, evidence)

        # Paralogue alignment route: the same amino-acid change at the ANALOGOUS
        # residue of a paralogue, translated via the alignment table.
        #  * SCN1A/2A/3A/8A (GN067-070): >=2 paralogue Pathogenic → PS1_Strong;
        #    a single P/LP → PS1_Supporting.
        #  * KCNQ1 (GN112): a Pathogenic KCNQ2 same-AA variant → PS1_Moderate
        #    (fixed; LP-only does not qualify). The fixed strength is taken from
        #    ps1_paralog_strength.
        aa = _aa_ref_alt(pc.amino_acids)
        if aa and self._paralog_map.has_gene(pc.gene_symbol):
            from acmg_classifier.local_db.clinvar_sqlite import query_paralog_aa_change
            ref, alt = aa
            analogs = self._paralog_map.analogs(pc.gene_symbol, pc.protein_position)
            path_genes: set[str] = set()
            any_plp: list = []
            for sib, spos in analogs.items():
                hits = query_paralog_aa_change(
                    self._cfg.clinvar_sqlite, sib, f"{ref}{spos}{alt}", min_stars=self._min_stars,
                )
                for h in hits:
                    any_plp.append(h)
                    sig = (h.clinical_significance or "").lower()
                    if "pathogenic" in sig and sig != "likely pathogenic":
                        path_genes.add(sib)
            fixed = self._spec.paralog_strength(pc.gene_symbol)
            strength = None
            if fixed is not None:
                # KCNQ1: requires a Pathogenic paralogue hit (not LP-only).
                if path_genes:
                    strength = fixed
            elif any_plp:
                strength = (CriterionStrength.STRONG if len(path_genes) >= 2
                            else CriterionStrength.SUPPORTING)
            if strength is not None:
                strength = _cap(strength, self._spec.max_strength(pc.gene_symbol))
                ids = ", ".join(h.variation_id or "" for h in any_plp[:3])
                return CriteriaResult.met(
                    ACMGCriterion.PS1, strength,
                    f"ClinVar same AA at analogous paralogue residue "
                    f"({len(path_genes)} paralogue P / {len(any_plp)} P-LP), "
                    f"{strength.value}: {ids}",
                )

        return CriteriaResult.not_met(
            ACMGCriterion.PS1,
            f"No ClinVar >={self._min_stars} star same-AA hit (excluding self)"
            + (f" [{vcep_note}]" if vcep_note else ""),
        )

    def _evaluate_splice(self, variant: VariantRecord, pc) -> CriteriaResult:
        # PS1's splice extension is opt-in per VCEP. Genes whose PS1 is the
        # original missense-only ACMG rule (no splice extension — e.g. GAA, the
        # HCM genes) must NOT receive PS1 for a splice/intronic variant.
        mode = self._spec.splice_mode(pc.gene_symbol)
        if not mode:
            return CriteriaResult.not_met(
                ACMGCriterion.PS1,
                f"{pc.gene_symbol}: PS1 is missense-only (no splice extension)",
            )
        # A "noncanonical" extension excludes canonical ±1/±2 sites
        # (SPLICE_DONOR / SPLICE_ACCEPTOR) — those are PVS1 territory.
        if mode == "noncanonical" and pc.consequence in (
            ConsequenceType.SPLICE_ACCEPTOR,
            ConsequenceType.SPLICE_DONOR,
        ):
            return CriteriaResult.not_met(
                ACMGCriterion.PS1,
                f"{pc.gene_symbol}: PS1 splice limited to non-canonical sites "
                "(canonical ±1/±2 is PVS1)",
            )

        # Comparators (ClinGen SVI splicing, Walker et al. 2023, Table 2):
        #  * a different nucleotide change at the SAME position as a P/LP variant;
        #  * for a variant at a canonical ±1,2 position, also a P/LP variant at the
        #    other base of the same ±1,2 dinucleotide.
        # Comparators elsewhere in the splice motif are not used. The strength
        # depends on PVS1 / PP3 of the variant and is set by the registry
        # (``apply_ps1_splice_table2``) once those are final; here a provisional
        # Strong (Pathogenic comparator) / Moderate (LP only) is returned with a
        # machine-readable tag.
        from acmg_classifier.local_db.clinvar_sqlite import (
            query_same_splice_site, query_splice_dinucleotide_partner,
        )
        hits = query_same_splice_site(
            self._cfg.clinvar_sqlite,
            pc.gene_symbol,
            variant.chrom, variant.pos, variant.ref, variant.alt,
            min_stars=self._min_stars,
        )
        if pc.consequence == ConsequenceType.SYNONYMOUS:
            # Exonic variant: only comparators that act through splicing (no
            # amino-acid change) predict the same splicing event.
            hits = [h for h in hits if _no_protein_change(h)]
        canon = canonical_offset(pc)
        if canon is not None and len(variant.ref) == 1 and len(variant.alt) == 1:
            seen = {h.variation_id for h in hits}
            hits += [h for h in query_splice_dinucleotide_partner(
                self._cfg.clinvar_sqlite, pc.gene_symbol, variant.chrom, variant.pos,
                canon[0], canon[1], min_stars=self._min_stars,
            ) if h.variation_id not in seen]
        from acmg_classifier.criteria.vcep_comparators import filter_comparators
        hits, vcep_note = filter_comparators(
            pc.gene_symbol, "PS1", hits, variant.assembly.value,
        )
        if not hits:
            return CriteriaResult.not_met(
                ACMGCriterion.PS1,
                f"No ClinVar >={self._min_stars} star same-splice-site P/LP hit"
                + (f" [{vcep_note}]" if vcep_note else ""),
            )
        provisional = _ps1_strength(hits)
        comp = "P" if provisional == CriterionStrength.STRONG else "LP"
        is_canon = canon is not None or pc.consequence in (
            ConsequenceType.SPLICE_DONOR, ConsequenceType.SPLICE_ACCEPTOR)
        site = "canonical" if is_canon else "noncanonical"
        strength = _cap(provisional, self._spec.max_strength(pc.gene_symbol))
        evidence = (
            f"ClinVar same splice-site position (different nucleotide"
            f"{'; same ±1,2 dinucleotide' if canon is not None else ''}), "
            f"comparator {comp}: " + ", ".join(h.variation_id or "" for h in hits[:3])
            + f" {SPLICE_TAG}{site}:{comp}]"
            + (f" [{vcep_note}]" if vcep_note else "")
        )
        return CriteriaResult.met(ACMGCriterion.PS1, strength, evidence)
