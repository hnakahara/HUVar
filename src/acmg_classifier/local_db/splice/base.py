"""Abstract base class for splice predictors."""
from __future__ import annotations
from abc import ABC, abstractmethod

from acmg_classifier.models.annotation import SpliceScore
from acmg_classifier.models.variant import VariantRecord


class SplicePredictor(ABC):
    @abstractmethod
    def predict(self, variant: VariantRecord) -> SpliceScore: ...

    @abstractmethod
    def is_available(self) -> bool: ...
