"""Manual benign criteria: BS3, BS4, BP2, BP5."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Benign criteria that cannot be automated from public databases:
#   BS3 = functional studies show no damaging effect (requires literature)
#   BS4 = lack of segregation in affected family members
#   BP2 = observed in trans with pathogenic in dominant gene / in cis
#   BP5 = found in case with alternate molecular basis
# Same supplement-based curator workflow as the pathogenic manual.py.
_MANUAL_CRITERIA = (
    ACMGCriterion.BS3,
    ACMGCriterion.BS4,
    ACMGCriterion.BP2,
    ACMGCriterion.BP5,
)


class ManualBenignEvaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> list[CriteriaResult]:  # type: ignore[override]
        # Returns a list (one record per criterion) so the registry sees a
        # not_met entry for criteria the curator did not provide — keeps the
        # audit trail complete.
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
