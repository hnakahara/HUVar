"""PM4 -- protein length change due to in-frame indel or stop-loss."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class PM4Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.PM4, "No consequence")

        if pc.consequence in (
            ConsequenceType.INFRAME_INSERTION,
            ConsequenceType.INFRAME_DELETION,
            ConsequenceType.STOP_LOST,
        ):
            if annotation.repeat and annotation.repeat.in_repeat:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM4,
                    f"In-frame indel in repeat ({annotation.repeat.repeat_class}); use BP3",
                )
            return CriteriaResult.met(
                ACMGCriterion.PM4,
                evidence=f"{pc.consequence.value} outside repeat region",
            )
        return CriteriaResult.not_met(ACMGCriterion.PM4, "Not an in-frame indel or stop-loss")
