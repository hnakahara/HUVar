"""Abstract base for splice predictors.

The orchestrator selects between SQUIRLS and SpliceAI at runtime based on
config; both must expose the same `predict()` signature so PVS1/PP3/BP4/BP7
can consume splice scores without caring which tool produced them. SpliceAI
requires a commercial Illumina licence — that is why SQUIRLS exists as the
default open-source fallback rather than as an always-on alternative.
"""
from __future__ import annotations
from abc import ABC, abstractmethod

from acmg_classifier.models.annotation import SpliceScore
from acmg_classifier.models.variant import VariantRecord


class SplicePredictor(ABC):
    @abstractmethod
    def predict(self, variant: VariantRecord) -> SpliceScore: ...

    @abstractmethod
    def is_available(self) -> bool: ...
