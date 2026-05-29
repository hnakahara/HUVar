from __future__ import annotations
from typing import Optional
from pydantic import BaseModel

from acmg_classifier.models.enums import ACMGCriterion, CriterionDirection, CriterionStrength

# Bayesian point values per strength level (Tavtigian 2020 + Bergquist 2024).
# The Bayesian formulation maps each strength level to an integer "odds-of-
# pathogenicity" exponent so that evidence can be summed instead of resolved
# through the legacy 2015 rule combinatorics. Sums map to categories via the
# thresholds in classifier_bayesian.py. Benign criteria contribute negative
# points so a single sum captures both directions.
STRENGTH_POINTS: dict[tuple[CriterionStrength, CriterionDirection], int] = {
    (CriterionStrength.VERY_STRONG, CriterionDirection.PATHOGENIC): 8,
    (CriterionStrength.STRONG, CriterionDirection.PATHOGENIC): 4,
    (CriterionStrength.THREE_POINT, CriterionDirection.PATHOGENIC): 3,
    (CriterionStrength.MODERATE, CriterionDirection.PATHOGENIC): 2,
    (CriterionStrength.SUPPORTING, CriterionDirection.PATHOGENIC): 1,
    (CriterionStrength.VERY_STRONG, CriterionDirection.BENIGN): -8,
    (CriterionStrength.STRONG, CriterionDirection.BENIGN): -4,
    (CriterionStrength.THREE_POINT, CriterionDirection.BENIGN): -3,
    (CriterionStrength.MODERATE, CriterionDirection.BENIGN): -2,
    (CriterionStrength.SUPPORTING, CriterionDirection.BENIGN): -1,
    (CriterionStrength.NOT_MET, CriterionDirection.PATHOGENIC): 0,
    (CriterionStrength.NOT_MET, CriterionDirection.BENIGN): 0,
    (CriterionStrength.INDETERMINATE, CriterionDirection.PATHOGENIC): 0,
    (CriterionStrength.INDETERMINATE, CriterionDirection.BENIGN): 0,
}

# Default strength for each criterion per ACMG 2015 + ClinGen SVI
DEFAULT_STRENGTH: dict[ACMGCriterion, tuple[CriterionStrength, CriterionDirection]] = {
    ACMGCriterion.PVS1: (CriterionStrength.VERY_STRONG, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PS1: (CriterionStrength.STRONG, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PS2: (CriterionStrength.STRONG, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PS3: (CriterionStrength.STRONG, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PS4: (CriterionStrength.STRONG, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PM1: (CriterionStrength.MODERATE, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PM2: (CriterionStrength.SUPPORTING, CriterionDirection.PATHOGENIC),  # SVI update: Supporting
    ACMGCriterion.PM3: (CriterionStrength.MODERATE, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PM4: (CriterionStrength.MODERATE, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PM5: (CriterionStrength.MODERATE, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PM6: (CriterionStrength.MODERATE, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PP1: (CriterionStrength.SUPPORTING, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PP2: (CriterionStrength.SUPPORTING, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PP3: (CriterionStrength.SUPPORTING, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PP4: (CriterionStrength.SUPPORTING, CriterionDirection.PATHOGENIC),
    ACMGCriterion.PP5: (CriterionStrength.SUPPORTING, CriterionDirection.PATHOGENIC),
    ACMGCriterion.BA1: (CriterionStrength.VERY_STRONG, CriterionDirection.BENIGN),
    ACMGCriterion.BS1: (CriterionStrength.STRONG, CriterionDirection.BENIGN),
    ACMGCriterion.BS2: (CriterionStrength.STRONG, CriterionDirection.BENIGN),
    ACMGCriterion.BS3: (CriterionStrength.STRONG, CriterionDirection.BENIGN),
    ACMGCriterion.BS4: (CriterionStrength.STRONG, CriterionDirection.BENIGN),
    ACMGCriterion.BP1: (CriterionStrength.SUPPORTING, CriterionDirection.BENIGN),
    ACMGCriterion.BP2: (CriterionStrength.SUPPORTING, CriterionDirection.BENIGN),
    ACMGCriterion.BP3: (CriterionStrength.SUPPORTING, CriterionDirection.BENIGN),
    ACMGCriterion.BP4: (CriterionStrength.SUPPORTING, CriterionDirection.BENIGN),
    ACMGCriterion.BP5: (CriterionStrength.SUPPORTING, CriterionDirection.BENIGN),
    ACMGCriterion.BP6: (CriterionStrength.SUPPORTING, CriterionDirection.BENIGN),
    ACMGCriterion.BP7: (CriterionStrength.SUPPORTING, CriterionDirection.BENIGN),
}


class CriteriaResult(BaseModel):
    """Result of evaluating a single ACMG criterion for one variant."""

    criterion: ACMGCriterion
    triggered: bool
    strength: CriterionStrength
    direction: CriterionDirection
    evidence: str = ""             # human-readable explanation
    suppressed: bool = False       # e.g. PP3 suppressed when PVS1 triggered

    @property
    def points(self) -> int:
        """Bayesian contribution for this criterion (0 when not triggered or
        suppressed). Suppression is used when a higher-tier criterion already
        encodes the same evidence — e.g. PP3 is suppressed when PVS1 fires so
        in-silico evidence is not double-counted."""
        if not self.triggered or self.suppressed:
            return 0
        return STRENGTH_POINTS.get((self.strength, self.direction), 0)

    @classmethod
    def not_met(
        cls,
        criterion: ACMGCriterion,
        evidence: str = "",
    ) -> "CriteriaResult":
        """Build a "not triggered" result that still records the criterion's
        direction from DEFAULT_STRENGTH. We always emit a CriteriaResult per
        criterion (rather than omitting it) so downstream code can report the
        full evidence trail, including which criteria were explicitly checked
        and rejected."""
        default_strength, direction = DEFAULT_STRENGTH[criterion]
        return cls(
            criterion=criterion,
            triggered=False,
            strength=CriterionStrength.NOT_MET,
            direction=direction,
            evidence=evidence,
        )

    @classmethod
    def met(
        cls,
        criterion: ACMGCriterion,
        strength: Optional[CriterionStrength] = None,
        evidence: str = "",
    ) -> "CriteriaResult":
        """Build a "triggered" result. `strength` is optional because most
        criteria fire at their ACMG default level; only criteria with formal
        strength modifiers (e.g. PVS1_Strong, PS1_Moderate per ClinGen SVI)
        need to pass an explicit override."""
        default_strength, direction = DEFAULT_STRENGTH[criterion]
        return cls(
            criterion=criterion,
            triggered=True,
            strength=strength or default_strength,
            direction=direction,
            evidence=evidence,
        )
