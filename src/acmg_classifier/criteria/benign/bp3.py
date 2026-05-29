"""BP3 -- in-frame indel in repetitive region without known function."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class BP3Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # BP3 is the benign-side counterpart to PM4: in-frame indels in
        # repetitive regions are likely tolerated because length variation
        # in repeats is biologically common. PM4 explicitly hands these
        # variants to BP3 (see pm4.py).
        pc = annotation.primary_consequence
        if pc is None or pc.consequence not in (
            ConsequenceType.INFRAME_INSERTION,
            ConsequenceType.INFRAME_DELETION,
        ):
            return CriteriaResult.not_met(ACMGCriterion.BP3, "Not an in-frame indel")

        # Repeat overlap is determined by Dfam/RepeatMasker (see
        # local_db.repeatmasker_db). The "without known function" qualifier
        # from the ACMG text is not enforced here — we assume the repeat
        # tracks already exclude domains with known function.
        rep = annotation.repeat
        if rep and rep.in_repeat:
            return CriteriaResult.met(
                ACMGCriterion.BP3,
                evidence=f"In-frame indel in repeat element ({rep.repeat_class}: {rep.repeat_name})",
            )
        return CriteriaResult.not_met(ACMGCriterion.BP3, "Not in a repeat region (Dfam)")
