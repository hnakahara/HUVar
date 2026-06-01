"""BP7 -- synonymous/deep-intronic variant with no predicted splice impact (Walker 2023 expansion)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Walker 2023 splice-impact safe-distance thresholds. Variants closer to
# the splice site are inside the canonical donor/acceptor consensus and
# cannot satisfy BP7 even with a low splice-tool score. The asymmetric
# cutoffs (+7 vs -21) reflect biology: branch points and polypyrimidine
# tracts extend much further upstream of the acceptor than downstream of
# the donor.
_DEEP_INTRONIC_DONOR_MIN = 7       # >= +7  (downstream of donor)
_DEEP_INTRONIC_ACCEPTOR_MAX = -21  # <= -21 (upstream of acceptor)


def _is_deep_intronic(pc) -> bool:
    """Variant is far enough from any canonical splice site to be plausibly
    benign under Walker 2023. Returns False when distance is unknown — we
    refuse to assume safety without evidence."""
    dist = pc.intron_distance_from_splice
    if dist is None:
        return False
    return dist >= _DEEP_INTRONIC_DONOR_MIN or dist <= _DEEP_INTRONIC_ACCEPTOR_MAX


def _splice_benign(annotation: AnnotationData) -> bool:
    """Splice predictor agrees the variant has no impact.

    Walker 2023 calibrates SpliceAI ≤ 0.10 as "no impact" (Supporting BP4).
    SQUIRLS lacks an equivalent calibration so we use < 0.20 as a practical
    approximation matching the lower BP4 bound — see README caveat."""
    sp = annotation.splice
    if sp is None or not sp.is_available:
        return False
    if sp.tool == "spliceai" and sp.max_delta is not None:
        return sp.max_delta <= 0.10
    if sp.tool == "squirls" and sp.raw_score is not None:
        return sp.raw_score < 0.20
    # MMSplice DISABLED — retained, commented out, for later:
    # if sp.tool == "mmsplice" and sp.raw_score is not None:
    #     # |delta_logit_psi| < 0.5 → minimal predicted splice effect (matches BP4).
    #     return abs(sp.raw_score) < 0.5
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

        # Synonymous branch: requires a splice-tool agreement because a
        # silent codon change can still create/destroy a splice site. We do
        # NOT auto-fire on "synonymous alone" — that was the pre-Walker
        # behaviour and is now known to mis-classify exonic splice variants.
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

        # Intronic branch: distance alone is enough to fire if the variant
        # is well outside the splice consensus. If a splice predictor *also*
        # agrees the variant is benign we still fire BP7 (just with stronger
        # evidence in the message) — distance is independently sufficient.
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
