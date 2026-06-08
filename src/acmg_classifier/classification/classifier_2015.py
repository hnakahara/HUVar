"""Rule-based ACMG 2015 classification."""
from __future__ import annotations
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion, CriterionDirection, CriterionStrength, Pathogenicity,
)


def _count(
    results: list[CriteriaResult],
    criterion: ACMGCriterion,
    strength: CriterionStrength,
) -> int:
    """Count results for one (criterion, strength) pair.

    NOTE: currently unused by the public Classifier2015.classify path —
    retained for tests / future strength-aware rules. See
    docs/cleanup-candidates.md."""
    return sum(
        1 for r in results
        if r.criterion == criterion
        and r.triggered
        and not r.suppressed
        and r.strength == strength
    )


def _has(results: list[CriteriaResult], criterion: ACMGCriterion) -> bool:
    """Convenience predicate. Currently unused — see _count note."""
    return any(
        r.criterion == criterion and r.triggered and not r.suppressed
        for r in results
    )


def _bucket(direction: CriterionDirection, strength: CriterionStrength) -> str | None:
    """Map a (direction, strength) pair to an ACMG 2015 evidence bucket.

    The 2015 combination table (Richards et al. 2015, Table 5) only knows
    four native strength tiers per direction:
        Pathogenic : Very Strong (pvs) > Strong (ps) > Moderate (pm) > Supporting (pp)
        Benign     : Stand-alone (ba)  > Strong (bs) >       —       > Supporting (bp)
    so we bucket by the *actual* strength a criterion fired at — not its name.
    This is what makes e.g. PM2 (SVI default: Supporting) count as Supporting
    rather than Moderate, and lets PP3/BP4 fired at Moderate/Strong (Bergquist
    2024 in-silico tiers) count at their true level.

    THREE_POINT (a Bergquist 2024 extension, +3) has no native 2015 tier, and
    the benign side has no Moderate tier at all. Both are resolved by rounding
    DOWN to the nearest weaker native tier — the conservative choice that never
    inflates evidence beyond what the original framework can express:
        pathogenic ThreePoint   -> Moderate (pm)
        benign ThreePoint/Moderate -> Supporting (bp)

    Returns None for NOT_MET / INDETERMINATE (no contribution)."""
    if direction == CriterionDirection.PATHOGENIC:
        if strength == CriterionStrength.VERY_STRONG:
            return "pvs"
        if strength == CriterionStrength.STRONG:
            return "ps"
        if strength in (CriterionStrength.THREE_POINT, CriterionStrength.MODERATE):
            return "pm"
        if strength == CriterionStrength.SUPPORTING:
            return "pp"
        return None
    # Benign
    if strength == CriterionStrength.VERY_STRONG:
        return "ba"
    if strength == CriterionStrength.STRONG:
        return "bs"
    if strength in (
        CriterionStrength.THREE_POINT,
        CriterionStrength.MODERATE,
        CriterionStrength.SUPPORTING,
    ):
        return "bp"
    return None


def _triggered(results: list[CriteriaResult]) -> dict[str, list[str]]:
    """Group triggered criteria into ACMG 2015 evidence buckets.

    The buckets (pvs/ps/pm/pp/ba/bs/bp) are what the rule table in
    Classifier2015.classify counts against. Every criterion is bucketed by the
    *strength it actually fired at* (see _bucket), not by its name prefix — so a
    strength-modified criterion (PVS1 capped to Moderate, PM2 at Supporting,
    PP3 promoted to Strong, etc.) is counted at its true tier."""
    out: dict[str, list[str]] = {
        "pvs": [], "ps": [], "pm": [], "pp": [],
        "ba": [], "bs": [], "bp": [],
    }
    for r in results:
        if not r.triggered or r.suppressed:
            continue
        bucket = _bucket(r.direction, r.strength)
        if bucket is not None:
            out[bucket].append(r.criterion.value)
    return out


class Classifier2015:
    """Implements the ACMG 2015 combination rules (Table 5)."""

    def classify(
        self, results: list[CriteriaResult]
    ) -> tuple[Pathogenicity, str]:
        """Apply ACMG 2015 Table 5 combination rules.

        Returns (final pathogenicity, human-readable rule string). The rule
        string lists every criterion that contributed, preserving the
        evidence trail so a clinician can audit the call without re-running
        the classifier."""
        g = _triggered(results)
        pvs = len(g["pvs"])
        ps = len(g["ps"])
        pm = len(g["pm"])
        pp = len(g["pp"])
        ba = len(g["ba"])
        bs = len(g["bs"])
        bp = len(g["bp"])

        triggered_names = (
            g["pvs"] + g["ps"] + g["pm"] + g["pp"]
            + g["ba"] + g["bs"] + g["bp"]
        )
        rule_str = " + ".join(triggered_names) if triggered_names else "none"

        # BA1 is "stand-alone" benign — a single trigger ends classification
        # before any pathogenic combination is considered. Per ACMG 2015 BA1
        # always wins over otherwise-pathogenic evidence (it represents
        # population-level frequency that excludes Mendelian disease).
        if ba:
            return Pathogenicity.BENIGN, rule_str

        # --- Pathogenic ---
        if (
            (pvs >= 1 and (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2))
            or (ps >= 2)
            or (ps >= 1 and (pm >= 3 or (pm >= 2 and pp >= 2) or (pm >= 1 and pp >= 4)))
        ):
            return Pathogenicity.PATHOGENIC, rule_str

        # --- Likely Pathogenic ---
        if (
            (pvs >= 1 and pm >= 1)
            or (ps >= 1 and pm >= 1)
            or (ps >= 1 and pp >= 2)
            or (pm >= 3)
            or (pm >= 2 and pp >= 2)
            or (pm >= 1 and pp >= 4)
        ):
            return Pathogenicity.LIKELY_PATHOGENIC, rule_str

        # --- Benign ---
        if bs >= 2 or (bs >= 1 and bp >= 1):
            # not quite Benign unless 2×Strong benign or 1×Strong + supporting
            # ACMG table says: ≥2 Strong OR (1 Strong + 1 Supporting)
            return Pathogenicity.BENIGN, rule_str

        # --- Likely Benign ---
        if bs >= 1 or bp >= 2:
            return Pathogenicity.LIKELY_BENIGN, rule_str

        return Pathogenicity.VUS, rule_str
