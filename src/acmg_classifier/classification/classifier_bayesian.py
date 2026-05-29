"""Bayesian point-based classification (Tavtigian 2020 + Bergquist 2024 extension)."""
from __future__ import annotations
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import Pathogenicity


class ClassifierBayesian:
    """
    Score thresholds:
        ≥ 10  → Pathogenic
        6–9   → Likely Pathogenic
        -5–5  → VUS
        -6–-9 → Likely Benign
        ≤ -10 → Benign
    """

    def classify(
        self, results: list[CriteriaResult]
    ) -> tuple[int, Pathogenicity]:
        """Sum signed strength points and threshold into a category.

        The Bayesian model collapses the legacy 2015 combination rules into
        a single signed-integer score (negative = benign-leaning, positive
        = pathogenic-leaning). Thresholds derive from the OddsPath
        calibration in Tavtigian 2020 — they are not arbitrary."""
        # BA1 is stand-alone benign per Tavtigian 2020 — it sits OUTSIDE the
        # point-sum framework. Even when other pathogenic evidence is
        # present, BA1's "population frequency too high for any Mendelian
        # disease" interpretation overrides. We still return the raw score
        # so reporting can show what the points would have said.
        from acmg_classifier.models.enums import ACMGCriterion
        if any(r.criterion == ACMGCriterion.BA1 and r.triggered and not r.suppressed for r in results):
            score = sum(r.points for r in results)
            return score, Pathogenicity.BENIGN

        # CriteriaResult.points already encodes sign (pathogenic positive,
        # benign negative) and suppression (returns 0). A simple sum is all
        # the Bayesian framework needs.
        score = sum(r.points for r in results)

        if score >= 10:
            classification = Pathogenicity.PATHOGENIC
        elif score >= 6:
            classification = Pathogenicity.LIKELY_PATHOGENIC
        elif score >= -5:
            classification = Pathogenicity.VUS
        elif score >= -9:
            classification = Pathogenicity.LIKELY_BENIGN
        else:
            classification = Pathogenicity.BENIGN

        return score, classification
