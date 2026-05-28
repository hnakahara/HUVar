"""BS2 -- observed in healthy adult individual (gnomAD homozygote/hemizygote counts)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_MIN_HOMALT = 5   # at least 5 homozygotes in gnomAD to trigger BS2
_MIN_HEMI = 5     # for X-linked


class BS2Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BS2, "No valid gnomAD record")

        nhomalt = gd.nhomalt or 0
        nhemi = gd.nhemi or 0

        if nhomalt >= _MIN_HOMALT:
            return CriteriaResult.met(
                ACMGCriterion.BS2,
                evidence=f"gnomAD nhomalt={nhomalt} >= {_MIN_HOMALT}",
            )
        if nhemi >= _MIN_HEMI:
            return CriteriaResult.met(
                ACMGCriterion.BS2,
                evidence=f"gnomAD nhemi={nhemi} >= {_MIN_HEMI}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BS2,
            f"gnomAD nhomalt={nhomalt}, nhemi={nhemi} (below threshold)",
        )
