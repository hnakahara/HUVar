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

    def precompute(self, variants: list[VariantRecord]) -> None:
        """Optionally pre-compute scores for the whole batch before predict().

        No-op by default: tabix-backed predictors (SpliceAI, SQUIRLS) look up
        precomputed scores per variant, so they need no batch step. Runtime
        predictors like MMSplice override this to run one expensive model pass
        over the entire batch and cache the results, which predict() then reads.
        """
        return None
