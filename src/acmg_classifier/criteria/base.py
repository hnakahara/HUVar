from __future__ import annotations
from abc import ABC, abstractmethod

from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class CriterionEvaluator(ABC):
    """Abstract base for all ACMG criterion evaluators."""

    @abstractmethod
    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult: ...
