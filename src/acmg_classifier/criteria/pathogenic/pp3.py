"""PP3 -- computational evidence of deleterious effect (Bergquist 2024 thresholds)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion, ConsequenceType, CriterionStrength, InSilicoTool,
)
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


def _alphamissense_pp3(score: float) -> CriterionStrength | None:
    if score >= 0.990:
        return CriterionStrength.STRONG
    if score >= 0.972:
        return CriterionStrength.THREE_POINT
    if score >= 0.906:
        return CriterionStrength.MODERATE
    if score >= 0.792:
        return CriterionStrength.SUPPORTING
    return None


def _esm1b_pp3(llr: float) -> CriterionStrength | None:
    """Bergquist 2024 Table 2 ESM1b PP3 thresholds (lower LLR ⇒ more pathogenic)."""
    if llr <= -24.0:
        return CriterionStrength.STRONG
    if llr <= -14.0:
        return CriterionStrength.THREE_POINT
    if llr <= -12.2:
        return CriterionStrength.MODERATE
    if llr <= -10.7:
        return CriterionStrength.SUPPORTING
    return None


def _squirls_pp3(score: float) -> CriterionStrength | None:
    if score >= 0.50:
        return CriterionStrength.MODERATE
    if score >= 0.20:
        return CriterionStrength.SUPPORTING
    return None


def _spliceai_pp3(max_delta: float) -> CriterionStrength | None:
    if max_delta >= 0.20:
        return CriterionStrength.MODERATE
    return None


class PP3Evaluator(CriterionEvaluator):
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
            return CriteriaResult.not_met(ACMGCriterion.PP3, "No consequence")

        if pc.consequence == ConsequenceType.MISSENSE:
            sp = annotation.splice
            if sp and sp.is_available and sp.tool == "spliceai" and sp.max_delta is not None:
                if sp.max_delta >= 0.20:
                    return CriteriaResult.met(
                        ACMGCriterion.PP3, CriterionStrength.MODERATE,
                        f"SpliceAI max_delta={sp.max_delta:.3f} (Moderate) — missense with predicted splice impact",
                    )
            if self._cfg.insilico_tool == InSilicoTool.ESM1B:
                es = annotation.esm1b
                if es and es.llr is not None:
                    strength = _esm1b_pp3(es.llr)
                    if strength:
                        return CriteriaResult.met(
                            ACMGCriterion.PP3, strength,
                            f"ESM1b LLR={es.llr:.3f} ({strength.value})",
                        )
                    return CriteriaResult.not_met(
                        ACMGCriterion.PP3,
                        f"ESM1b LLR={es.llr:.3f} (indeterminate or benign)",
                    )
                return CriteriaResult.not_met(ACMGCriterion.PP3, "No in-silico score available")

            am = annotation.alphamissense
            if am and am.score is not None:
                strength = _alphamissense_pp3(am.score)
                if strength:
                    return CriteriaResult.met(
                        ACMGCriterion.PP3, strength,
                        f"AlphaMissense={am.score:.3f} ({strength.value})",
                    )
                return CriteriaResult.not_met(
                    ACMGCriterion.PP3,
                    f"AlphaMissense={am.score:.3f} (indeterminate or benign)",
                )
            return CriteriaResult.not_met(ACMGCriterion.PP3, "No in-silico score available")

        if pc.consequence in (
            ConsequenceType.SPLICE_REGION,
            ConsequenceType.INTRON,
            ConsequenceType.SYNONYMOUS,
        ):
            sp = annotation.splice
            if sp and sp.is_available:
                if sp.tool == "spliceai" and sp.max_delta is not None:
                    strength = _spliceai_pp3(sp.max_delta)
                    score_str = f"SpliceAI max_delta={sp.max_delta:.3f}"
                elif sp.tool == "squirls" and sp.raw_score is not None:
                    strength = _squirls_pp3(sp.raw_score)
                    score_str = f"SQUIRLS={sp.raw_score:.3f} (thresholds approximate)"
                else:
                    return CriteriaResult.not_met(ACMGCriterion.PP3, "Splice score unavailable")

                if strength:
                    return CriteriaResult.met(
                        ACMGCriterion.PP3, strength,
                        f"{score_str} ({strength.value})",
                    )
            return CriteriaResult.not_met(ACMGCriterion.PP3, "Splice score not pathogenic")

        return CriteriaResult.not_met(ACMGCriterion.PP3, "Consequence not applicable for PP3")
