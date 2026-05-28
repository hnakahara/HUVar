"""BP4 -- computational evidence suggesting no impact (Bergquist 2024 thresholds)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


def _alphamissense_bp4(score: float) -> CriterionStrength | None:
    """AlphaMissense BP4 strength per Bergquist 2024 Table 2.

    No Strong (-4) category exists for AlphaMissense; the strongest BP4 the
    table assigns is ThreePoint at score ≤ 0.070.
    """
    if score <= 0.070:
        return CriterionStrength.THREE_POINT
    if score <= 0.099:
        return CriterionStrength.MODERATE
    if score <= 0.169:
        return CriterionStrength.SUPPORTING
    return None


def _spliceai_bp4(max_delta: float) -> CriterionStrength | None:
    """SpliceAI BP4 per Walker 2023."""
    if max_delta <= 0.10:
        return CriterionStrength.SUPPORTING
    return None


def _squirls_bp4(score: float) -> CriterionStrength | None:
    """SQUIRLS BP4 (approximate, not Walker 2023 calibrated)."""
    if score < 0.20:
        return CriterionStrength.SUPPORTING
    return None


class BP4Evaluator(CriterionEvaluator):
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
            return CriteriaResult.not_met(ACMGCriterion.BP4, "No consequence")

        if pc.consequence == ConsequenceType.MISSENSE:
            sp = annotation.splice
            if sp and sp.is_available and sp.tool == "spliceai" and sp.max_delta is not None:
                if sp.max_delta >= 0.20:
                    return CriteriaResult.not_met(
                        ACMGCriterion.BP4,
                        f"SpliceAI max_delta={sp.max_delta:.3f} — predicted splice impact, BP4 not applicable",
                    )
            am = annotation.alphamissense
            if am and am.score is not None:
                strength = _alphamissense_bp4(am.score)
                if strength:
                    return CriteriaResult.met(
                        ACMGCriterion.BP4, strength,
                        f"AlphaMissense={am.score:.3f} ({strength.value})",
                    )
                return CriteriaResult.not_met(
                    ACMGCriterion.BP4,
                    f"AlphaMissense={am.score:.3f} (not in benign range)",
                )

        if pc.consequence in (
            ConsequenceType.SPLICE_REGION,
            ConsequenceType.INTRON,
            ConsequenceType.SYNONYMOUS,
        ):
            sp = annotation.splice
            if sp and sp.is_available:
                if sp.tool == "spliceai" and sp.max_delta is not None:
                    strength = _spliceai_bp4(sp.max_delta)
                    score_str = f"SpliceAI max_delta={sp.max_delta:.3f}"
                elif sp.tool == "squirls" and sp.raw_score is not None:
                    strength = _squirls_bp4(sp.raw_score)
                    score_str = f"SQUIRLS={sp.raw_score:.3f} (approximate)"
                else:
                    return CriteriaResult.not_met(ACMGCriterion.BP4, "Splice score unavailable")

                if strength:
                    return CriteriaResult.met(
                        ACMGCriterion.BP4, strength,
                        f"{score_str} ({strength.value})",
                    )
            return CriteriaResult.not_met(ACMGCriterion.BP4, "Splice score not in benign range")

        return CriteriaResult.not_met(ACMGCriterion.BP4, "Consequence not applicable for BP4")
