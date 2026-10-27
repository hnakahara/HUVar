from __future__ import annotations
from abc import ABC, abstractmethod

from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class CriterionEvaluator(ABC):
    """Abstract base for all ACMG criterion evaluators."""

    # True for criteria that are entirely ClinVar-codon / VCEP-position based
    # (PM1, PS1, PM5). When the criteria are evaluated on a non-default
    # transcript (MANE Plus Clinical), the registry hands these evaluators
    # ``annotation.rule_view()`` so they read MANE Select numbering — the
    # numbering of the VCEP tables and of ClinVar records.
    rule_numbered: bool = False

    @abstractmethod
    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult: ...
