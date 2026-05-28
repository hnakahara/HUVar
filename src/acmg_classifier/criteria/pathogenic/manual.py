"""Manual pathogenic criteria: PS2, PM3, PM6, PP4.

PS3 (functional) and PP1 (cosegregation) have dedicated evaluators that mine ClinVar
SCV comments (with manual supplement override), so they are no longer handled here.
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_MANUAL_CRITERIA = (
    ACMGCriterion.PS2,
    ACMGCriterion.PM3,
    ACMGCriterion.PM6,
    ACMGCriterion.PP4,
)


class ManualPathogenicEvaluator(CriterionEvaluator):
    """Reads manual evidence from supplement TSV for criteria requiring human curation."""

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> list[CriteriaResult]:  # type: ignore[override]
        results = []
        sup = supplement or []
        for criterion in _MANUAL_CRITERIA:
            entries = [e for e in sup if e.criterion == criterion]
            if entries:
                entry = entries[0]
                results.append(CriteriaResult.met(criterion, entry.strength, entry.evidence))
            else:
                results.append(CriteriaResult.not_met(criterion, "No manual evidence provided"))
        return results
