"""BS1 -- allele frequency greater than expected for disorder."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.criteria.allele_frequency import DiseaseThresholds


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
        if faf is None:
            faf = gd.popmax_af or gd.af or 0.0

        # Per-gene threshold lookup: if the gene appears in the prevalence
        # table, use its disease-specific cutoff; otherwise the default applies.
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else ""
        threshold = self._thresholds.get(gene).bs1

        # Format with 3 significant figures (not fixed 4 decimals): disease-
        # specific cutoffs are often ~1e-4, where ".4f" rounds both the FAF and
        # the threshold to misleadingly coarse values (0.000358 → "0.0004",
        # 0.000316 → "0.0003"). The comparison itself uses full float precision.
        if faf >= threshold:
            return CriteriaResult.met(
                ACMGCriterion.BS1,
                evidence=f"gnomAD FAF95={faf:.3g} >= BS1 threshold {threshold:.3g} for {gene or 'default'}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BS1,
            f"gnomAD FAF95={faf:.3g} < BS1 threshold {threshold:.3g}",
        )
