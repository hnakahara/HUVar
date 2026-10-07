"""Shared helpers for criteria that come from curated sources rather than from
automated evaluation (Genome Medicine revision).

Priority for PS3 / PS4 / PP1: (1) curator-supplied supplement (user curation or
the precompiled eRepo supplement), (2) the criterion as applied by a ClinGen
expert panel to the same variant in ClinVar (>=3 stars), (3) automated evidence
where still allowed (PS3 text mining, capped at Supporting). The evidence text
is prefixed with the source label so users can see where each call came from.
"""
from __future__ import annotations

from typing import Optional

from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.supplement import SupplementEntry

SRC_SUPPLEMENT = "[curated supplement]"
SRC_EXPERT_PANEL = "[ClinVar expert panel]"
SRC_TEXT_MINING = "[ClinVar text mining]"
SRC_INFO = "[info only — not applied]"

_STRENGTH = {
    "VeryStrong": CriterionStrength.VERY_STRONG,
    "Strong": CriterionStrength.STRONG,
    "Moderate": CriterionStrength.MODERATE,
    "Supporting": CriterionStrength.SUPPORTING,
}


def from_supplement(criterion: ACMGCriterion,
                    supplement: Optional[list[SupplementEntry]]) -> Optional[CriteriaResult]:
    for e in supplement or []:
        if e.criterion == criterion:
            return CriteriaResult.met(criterion, e.strength, f"{SRC_SUPPLEMENT} {e.evidence}")
    return None


def from_expert_panel(criterion: ACMGCriterion, cfg, variant) -> Optional[CriteriaResult]:
    if getattr(cfg, "exclude_self_expert_panel", False) is True:
        return None
    from acmg_classifier.local_db.clinvar_sqlite import query_expert_panel_criteria
    ep = query_expert_panel_criteria(cfg.clinvar_sqlite, variant.chrom, variant.pos,
                                     variant.ref, variant.alt)
    st = ep.get(criterion.value)
    if not st:
        return None
    return CriteriaResult.met(
        criterion, _STRENGTH.get(st, CriterionStrength.STRONG),
        f"{SRC_EXPERT_PANEL} {criterion.value}_{st} applied by the ClinGen expert panel "
        "for this variant (ClinVar review status >=3 stars)",
    )
