"""Bayesian point-based classification (Tavtigian 2020 + Bergquist 2024 extension)."""
from __future__ import annotations
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import Pathogenicity


class ClassifierBayesian:
    """
    Score thresholds (Tavtigian 2020, prior P=0.10 — asymmetric: the benign side
    needs less evidence because the prior already favours benign):
        ≥ 10   → Pathogenic
        6–9    → Likely Pathogenic
        0–5    → VUS
        -1–-6  → Likely Benign
        ≤ -7   → Benign
    """

    def classify(
        self, results: list[CriteriaResult]
    ) -> tuple[int, Pathogenicity]:
        """Sum signed strength points and threshold into a category.

        The Bayesian model collapses the legacy 2015 combination rules into
        a single signed-integer score (negative = benign-leaning, positive
        = pathogenic-leaning). Thresholds derive from the OddsPath
        calibration in Tavtigian 2020 — they are not arbitrary."""
        # BA1 is NOT an override here: it contributes its points (-8, like a
        # Very Strong benign criterion) and the class is decided by the final
        # score like every other criterion. BA1 alone (-8) is Benign; BA1 with
        # strong pathogenic evidence can therefore end as VUS (conflicting
        # evidence). The 2015 Table 5 classifier keeps BA1 as stand-alone Benign.
        # CriteriaResult.points already encodes sign (pathogenic positive,
        # benign negative) and suppression (returns 0). A simple sum is all
        # the Bayesian framework needs.
        score = sum(r.points for r in results)

        if score >= 10:
            classification = Pathogenicity.PATHOGENIC
        elif score >= 6:
            classification = Pathogenicity.LIKELY_PATHOGENIC
        elif score >= 0:
            classification = Pathogenicity.VUS
        elif score >= -6:
            classification = Pathogenicity.LIKELY_BENIGN
        else:
            classification = Pathogenicity.BENIGN

        return score, classification
