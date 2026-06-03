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

        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else ""
        gt = self._thresholds.get(gene)
        threshold = gt.ba1

        # Same prefer-FAF95 logic as PM2/PS4: FAF gives the most conservative
        # estimate; only fall back to raw AF if FAF is unavailable.
        faf = gd.faf95_popmax
        if faf is None:
            faf = gd.popmax_af or gd.af or 0.0

        # X-linked genes whose VCEP defines the cutoff "in males" (RPGR, RS1,
        # ABCD1, SLC6A8, OTC): compare against the male (XY) allele frequency.
        # Fall back to the overall FAF when AF_XY is unavailable (e.g. a gnomAD
        # DB built before the af_xy column).
        metric = "gnomAD FAF95_popmax"
        if gt.af_basis == "males" and gd.af_xy is not None:
            faf = gd.af_xy
            metric = "gnomAD AF_XY (males)"

        # 3 significant figures rather than 4 fixed decimals — disease-specific
        # BA1 cutoffs can be ~1e-3/1e-5, which ".4f" rounds misleadingly. The
        # comparison uses full float precision regardless of this display.
        if faf >= threshold:
            return CriteriaResult.met(
                ACMGCriterion.BA1,
                evidence=f"{metric}={faf:.3g} >= BA1 threshold {threshold:.3g} for {gene or 'default'}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BA1,
            f"{metric}={faf:.3g} < BA1 threshold {threshold:.3g}",
        )
