"""BP7 -- synonymous/deep-intronic variant with no predicted splice impact (Walker 2023 expansion)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_DEEP_INTRONIC_DONOR_MIN = 7   # >= +7 from donor
_DEEP_INTRONIC_ACCEPTOR_MAX = -21  # <= -21 from acceptor


def _is_deep_intronic(pc) -> bool:
    """Return True if variant is far enough from canonical splice site (Walker 2023)."""
    dist = pc.intron_distance_from_splice
    if dist is None:
        return False
    return dist >= _DEEP_INTRONIC_DONOR_MIN or dist <= _DEEP_INTRONIC_ACCEPTOR_MAX


def _splice_benign(annotation: AnnotationData) -> bool:
    """Return True if splice tool predicts no impact."""
    sp = annotation.splice
    if sp is None or not sp.is_available:
        return False
    if sp.tool == "spliceai" and sp.max_delta is not None:
        return sp.max_delta <= 0.10  # Walker 2023
    if sp.tool == "squirls" and sp.raw_score is not None:
        return sp.raw_score < 0.20  # approximate
    return False


class BP7Evaluator(CriterionEvaluator):
    """
    Walker 2023 expanded BP7 applies to:
    1. Synonymous variants with no predicted splice impact
    2. Intronic variants >= +7 or <= -21 bp from canonical splice site
    """

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
            return CriteriaResult.not_met(ACMGCriterion.BP7, "No consequence")

        if pc.consequence == ConsequenceType.SYNONYMOUS:
            if _splice_benign(annotation):
                sp = annotation.splice
                return CriteriaResult.met(
                    ACMGCriterion.BP7,
                    evidence=f"Synonymous + {sp.tool} score benign",
                )
            return CriteriaResult.not_met(
                ACMGCriterion.BP7,
                "Synonymous but splice impact not ruled out",
            )

        if pc.consequence == ConsequenceType.INTRON:
            if _is_deep_intronic(pc):
                if _splice_benign(annotation):
                    return CriteriaResult.met(
                        ACMGCriterion.BP7,
                        evidence=(
                            f"Deep intronic (dist={pc.intron_distance_from_splice}) "
                            "and splice score benign (Walker 2023)"
                        ),
                    )
                return CriteriaResult.met(
                    ACMGCriterion.BP7,
                    evidence=f"Deep intronic (dist={pc.intron_distance_from_splice}; Walker 2023)",
                )

        return CriteriaResult.not_met(
            ACMGCriterion.BP7,
            "Not a synonymous or deep-intronic variant",
        )
