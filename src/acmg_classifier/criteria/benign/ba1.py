"""BA1 -- allele frequency >5% in gnomAD (stand-alone benign)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# The 5% threshold is the original ACMG 2015 stand-alone benign cutoff.
# It is intentionally NOT lowered for recessive genes (unlike PM2): if a
# variant is at 5% population frequency it cannot be a Mendelian disease
# allele under any inheritance model.
_BA1_FAF95_THRESHOLD = 0.05


class BA1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

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

        if faf >= _BA1_FAF95_THRESHOLD:
            return CriteriaResult.met(
                ACMGCriterion.BA1,
                evidence=f"gnomAD FAF95_popmax={faf:.4f} >= {_BA1_FAF95_THRESHOLD}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BA1,
            f"gnomAD FAF95_popmax={faf:.4f} < {_BA1_FAF95_THRESHOLD}",
        )
