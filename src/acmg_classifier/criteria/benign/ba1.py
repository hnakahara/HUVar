"""BA1 -- allele frequency >5% in gnomAD (stand-alone benign)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.criteria.allele_frequency import DiseaseThresholds


class BA1Evaluator(CriterionEvaluator):
    """BA1 cutoff is disease-specific (Whiffin/Ware: min(0.05, 10 x maxAF)),
    falling back to the ACMG 2015 stand-alone benign 5% when no per-gene
    parameters are available."""

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._thresholds = DiseaseThresholds(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # Absent or filter-failed gnomAD records cannot establish BA1 — we
        # need a *high* observed frequency, which requires a valid record.
        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BA1, "No valid gnomAD record")

        # Same prefer-FAF95 logic as PM2/PS4: FAF gives the most conservative
        # estimate; only fall back to raw AF if FAF is unavailable.
        faf = gd.faf95_popmax
        if faf is None:
            faf = gd.popmax_af or gd.af or 0.0

        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else ""
        threshold = self._thresholds.get(gene).ba1

        if faf >= threshold:
            return CriteriaResult.met(
                ACMGCriterion.BA1,
                evidence=f"gnomAD FAF95_popmax={faf:.4f} >= BA1 threshold {threshold:.4f} for {gene or 'default'}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BA1,
            f"gnomAD FAF95_popmax={faf:.4f} < BA1 threshold {threshold:.4f}",
        )
