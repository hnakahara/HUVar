"""Rule-based ACMG 2015 classification."""
from __future__ import annotations
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength, Pathogenicity


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


def _triggered(results: list[CriteriaResult]) -> dict[str, list[str]]:
    """Group triggered criteria into ACMG 2015 evidence buckets.

    The buckets (pvs/ps/pm/pp/ba/bs/bp) are what the rule table in
    Classifier2015.classify counts against. PVS1 with a *downgraded*
    strength (e.g. PVS1 capped to Moderate per the ClinGen SVI cap) is
    re-bucketed to the lower tier so the 2015 combination rules treat
    it correctly — otherwise a capped PVS1 would still count as 'pvs'
    and trigger Pathogenic with only one extra Moderate evidence."""
    out: dict[str, list[str]] = {
        "pvs": [], "ps": [], "pm": [], "pp": [],
        "ba": [], "bs": [], "bp": [],
    }
    for r in results:
        if not r.triggered or r.suppressed:
            continue
        c = r.criterion.value
        if c == "PVS1":
            # PVS1 is the only criterion with a sliding strength (per the
            # ClinGen 2019 decision tree). Re-bucket based on the actual
            # strength that fired, not the criterion name.
            if r.strength in (CriterionStrength.VERY_STRONG,):
                out["pvs"].append(c)
            elif r.strength == CriterionStrength.STRONG:
                out["ps"].append(c)
            elif r.strength == CriterionStrength.MODERATE:
                out["pm"].append(c)
            else:
                out["pp"].append(c)
        elif c.startswith("PS"):
            out["ps"].append(c)
        elif c.startswith("PM"):
            out["pm"].append(c)
        elif c.startswith("PP"):
            out["pp"].append(c)
        elif c == "BA1":
            out["ba"].append(c)
        elif c.startswith("BS"):
            out["bs"].append(c)
        elif c.startswith("BP"):
            out["bp"].append(c)
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
