"""BS1 -- allele frequency greater than expected for disorder."""
from __future__ import annotations
import re

from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.criteria.allele_frequency import DiseaseThresholds, variant_listed

# Bare protein change from an HGVS p. string, tolerant of a transcript prefix
# and ClinVar/VEP parentheses; "*" is normalised to "Ter" so p.Gln368* and
# p.Gln368Ter compare equal. Lower-cased for case-insensitive matching.
_PCHANGE_RE = re.compile(r"p\.\(?([A-Za-z]{3}\d+(?:[A-Za-z]{3}|Ter|\*|=))\)?")


def _norm_pchange(hgvs_p: str | None) -> str:
    if not hgvs_p:
        return ""
    m = _PCHANGE_RE.search(hgvs_p)
    return m.group(1).replace("*", "Ter").lower() if m else ""


# Genes whose BS1_Supporting tier is raised to the full BS1 strength when gnomAD
# shows at least this many homozygotes (HHT VCEP, ACVRL1 / ENG).
_BS1_SUPPORTING_HOM_UPGRADE: dict[str, int] = {"ACVRL1": 2, "ENG": 2}


class BS1Evaluator(CriterionEvaluator):
    """
    BS1 threshold is disorder-specific, derived from disease_prevalence.tsv
    (Whiffin/Ware maximum credible AF, floored at 0.05%). Falls back to a
    conservative flat default if the gene is absent or lacks parameters.
    """

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Thresholds are loaded once at instantiation rather than per-variant
        # to avoid re-parsing the TSV across thousands of variant evaluations.
        self._thresholds = DiseaseThresholds(cfg.disease_prevalence_tsv)
        # gnomAD coverage for the ENIGMA BRCA1/2 read-depth condition; optional.
        from acmg_classifier.local_db.coverage_db import CoverageDB
        try:
            self._coverage = CoverageDB(getattr(cfg, "gnomad_coverage_db", None))
        except Exception:  # noqa: BLE001 - no coverage DB -> gate skipped
            self._coverage = None

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BS1, "No valid gnomAD record")

        # Prefer FAF95 (the conservative 95% CI lower bound). Use it as-is when
        # present — including a legitimate 0.0 — and only fall back to raw
        # popmax/global AF when FAF95 is genuinely absent (None). This mirrors
        # BA1; the previous `faf95 or popmax or af` chain collapsed a real
        # FAF95 of 0.0 to the raw AF, over-firing BS1 on wide-CI rare variants.
        faf = gd.faf95_popmax
        metric = "gnomAD FAF95"
        if faf is None:
            faf = gd.popmax_af if gd.popmax_af is not None else (gd.af or 0.0)
            metric = "gnomAD popmax AF (FAF95 unavailable)"

        # Per-gene threshold lookup: if the gene appears in the prevalence
        # table, use its disease-specific cutoff; otherwise the default applies.
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else ""
        gt = self._thresholds.get(gene)
        threshold = gt.bs1

        # Variant-level BS1 exclusion: a recurrent disease allele the VCEP bars
        # from BS1 regardless of its population frequency (e.g. MYOC p.Gln368Ter,
        # whose ~2.6% allelic contribution is disease-driven, not benign).
        # Entries may be protein ("p.") or MANE transcript ("c.") changes — the
        # LGMD VCEP "BS1 exception" variants include intronic changes.
        # The VCEP lists are written in MANE Select numbering (rule_consequence).
        rc = annotation.rule_consequence
        hit = variant_listed(gt.bs1_exclude, rc.hgvs_c if rc else None, rc.hgvs_p if rc else None)
        if hit:
            return CriteriaResult.not_met(
                ACMGCriterion.BS1,
                f"{gene}: VCEP excludes BS1 for {hit} "
                "(recurrent / common pathogenic or hypomorphic allele)",
            )

        # Non-cancer subset (ENIGMA BRCA1/2): the VCEP judges BS1 on gnomAD's
        # non-cancer subset, so use its popmax FAF95 (recomputed at build time
        # from the per-group non-cancer AC/AN). Fall back to the overall FAF95
        # (noting it) when the companion non-cancer value is unavailable — e.g. a
        # companion DB predating the column, or a variant absent from the subset.
        if gt.af_subset == "non_cancer":
            if gd.faf95_non_cancer is not None:
                faf = gd.faf95_non_cancer
                metric = "gnomAD FAF95 (non-cancer)"
            else:
                metric += " [non-cancer subset unavailable → overall]"

        # X-linked "in males" genes (RPGR, RS1): the VCEP cutoff is on the male
        # (XY) allele frequency. Fall back to the overall FAF when AF_XY is
        # unavailable (gnomAD DB predating the af_xy column).
        if gt.af_basis == "males" and gd.af_xy is not None:
            faf = gd.af_xy
            metric = "gnomAD AF_XY (males)"

        # Point-estimate basis (mirrors BA1): VCEPs that define the cutoff on the
        # grpmax/popmax POINT AF rather than FAF95 (af_basis="popmax"), unless
        # globally disabled (Config.popmax_af_basis=False).
        if (gt.af_basis == "popmax" and self._cfg.popmax_af_basis
                and gd.popmax_af is not None):
            faf = gd.popmax_af
            metric = "gnomAD popmax AF (point)"

        # "99.99% CI" (CDH1, F8 / F9 VCEPs): the filtering AF at 99.99% confidence,
        # from the GrpMax population's allele count / number. Without those counts
        # the gnomAD FAF95 (or the point AF) stays in use, and the metric says so.
        if gt.af_basis == "faf9999":
            from acmg_classifier.criteria.allele_frequency import filtering_af
            f9999 = filtering_af(gd.ac_grpmax, gd.an_grpmax)
            if f9999 is not None:
                faf = f9999
                metric = "gnomAD GrpMax filtering AF (99.99% CI)"
            else:
                metric += " [99.99% CI not computable: GrpMax AC/AN unavailable]"

        # No FAF95 (gnomAD computes it only when the grpmax group has more than one
        # observation): the ClinGen/VCEP frequency rules are defined on FAF95 in
        # populations with >1 observation, so a single-observation point AF must
        # not establish BS1. Only the explicit point-AF / male-AF bases above
        # replace FAF95.
        self.note(af_metric=metric, af_value=f"{faf:.3g}", bs1_cutoff=f"{threshold:.3g}")
        if metric.startswith("gnomAD popmax AF (FAF95 unavailable)"):
            return CriteriaResult.not_met(
                ACMGCriterion.BS1,
                f"gnomAD FAF95 unavailable (single observation; popmax AF={faf:.3g}); "
                f"BS1 not applied",
            )

        # Format with 3 significant figures (not fixed 4 decimals): disease-
        # specific cutoffs are often ~1e-4, where ".4f" rounds both the FAF and
        # the threshold to misleadingly coarse values (0.000358 → "0.0004",
        # 0.000316 → "0.0003"). The comparison itself uses full float precision.
        if faf >= threshold:
            from acmg_classifier.criteria.allele_frequency import min_alleles_block
            few = min_alleles_block(gene, gd)
            if not few:
                from acmg_classifier.criteria.allele_frequency import depth_block
                few = depth_block(self._coverage, gene, variant)
            if few:
                return CriteriaResult.not_met(
                    ACMGCriterion.BS1,
                    f"{metric}={faf:.3g} >= BS1 threshold {threshold:.3g} but {few}",
                )
            return CriteriaResult.met(
                ACMGCriterion.BS1,
                strength=gt.bs1_strength,
                evidence=f"{metric}={faf:.3g} >= BS1 threshold {threshold:.3g} "
                         f"({gt.bs1_strength.value}) for {gene or 'default'}",
            )
        sup = gt.bs1_supporting
        if sup is not None and faf >= sup:
            from acmg_classifier.criteria.allele_frequency import min_alleles_block
            few = min_alleles_block(gene, gd)
            if not few:
                from acmg_classifier.criteria.allele_frequency import depth_block
                few = depth_block(self._coverage, gene, variant)
            if few:
                return CriteriaResult.not_met(
                    ACMGCriterion.BS1,
                    f"{metric}={faf:.3g} >= BS1_Supporting threshold {sup:.3g} but {few}",
                )
            # HHT VCEP (ACVRL1 / ENG): BS1 (Strong) also "if variant meets
            # BS1_Supporting and has >= 2 homozygotes".
            need_hom = _BS1_SUPPORTING_HOM_UPGRADE.get(gene)
            if need_hom is not None and (gd.nhomalt or 0) >= need_hom:
                return CriteriaResult.met(
                    ACMGCriterion.BS1, strength=gt.bs1_strength,
                    evidence=(f"{metric}={faf:.3g} >= BS1_Supporting threshold {sup:.3g} and "
                              f"{gd.nhomalt} homozygotes (>= {need_hom}) "
                              f"({gt.bs1_strength.value}) for {gene}"),
                )
            return CriteriaResult.met(
                ACMGCriterion.BS1, strength=CriterionStrength.SUPPORTING,
                evidence=(f"{metric}={faf:.3g} >= BS1_Supporting threshold {sup:.3g} "
                          f"(< BS1 threshold {threshold:.3g}) for {gene}"),
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BS1,
            f"{metric}={faf:.3g} < BS1 threshold {threshold:.3g}",
        )
