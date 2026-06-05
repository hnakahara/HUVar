"""Unit tests for manual-evidence override across ALL criteria (both modes)."""
from pathlib import Path

import pytest

from acmg_classifier.config import Config
from acmg_classifier.criteria.registry import CriteriaRegistry
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength, SupplementMode
from acmg_classifier.models.supplement import SupplementEntry


def _tool_results() -> list[CriteriaResult]:
    """Automated calls a tool might emit: PVS1_strong, PS1, BP1, PM2_supporting."""
    return [
        CriteriaResult.met(ACMGCriterion.PVS1, CriterionStrength.STRONG, "tool PVS1"),
        CriteriaResult.met(ACMGCriterion.PS1, CriterionStrength.STRONG, "tool PS1"),
        CriteriaResult.met(ACMGCriterion.BP1, CriterionStrength.SUPPORTING, "tool BP1"),
        CriteriaResult.met(ACMGCriterion.PM2, CriterionStrength.SUPPORTING, "tool PM2"),
    ]


def _supplement() -> list[SupplementEntry]:
    return [
        SupplementEntry(variant_id="x", criterion=ACMGCriterion.PVS1,
                        strength=CriterionStrength.MODERATE, evidence="manual PVS1"),
        SupplementEntry(variant_id="x", criterion=ACMGCriterion.PS1,
                        strength=CriterionStrength.MODERATE, evidence="manual PS1"),
        SupplementEntry(variant_id="x", criterion=ACMGCriterion.PM2,
                        strength=CriterionStrength.SUPPORTING, evidence="manual PM2"),
    ]


def _registry(mode: SupplementMode) -> CriteriaRegistry:
    return CriteriaRegistry(Config(data_dir=Path("./data"), supplement_mode=mode))


def _by_crit(results):
    return {r.criterion: r for r in results}


def test_merge_overrides_strength_and_keeps_tool_only_criteria():
    reg = _registry(SupplementMode.MERGE)
    results = _tool_results()
    reg._apply_supplement_override(results, _supplement())
    by = _by_crit(results)

    # Curator strength overrides the tool's call on named criteria.
    assert by[ACMGCriterion.PVS1].strength == CriterionStrength.MODERATE
    assert by[ACMGCriterion.PS1].strength == CriterionStrength.MODERATE
    assert "manual override" in by[ACMGCriterion.PVS1].evidence
    # Tool-only criterion (no manual entry) is retained in merge mode.
    assert by[ACMGCriterion.BP1].triggered
    assert by[ACMGCriterion.BP1].strength == CriterionStrength.SUPPORTING


def test_manual_only_discards_tool_only_criteria():
    reg = _registry(SupplementMode.MANUAL_ONLY)
    results = _tool_results()
    reg._apply_supplement_override(results, _supplement())
    by = _by_crit(results)

    # Named criteria keep the curator's strength.
    assert by[ACMGCriterion.PVS1].strength == CriterionStrength.MODERATE
    assert by[ACMGCriterion.PM2].triggered
    # BP1 had no manual entry → dropped in manual-only mode.
    assert not by[ACMGCriterion.BP1].triggered
    assert by[ACMGCriterion.BP1].strength == CriterionStrength.NOT_MET


def test_merge_adds_criterion_without_evaluator():
    """A curator may name a criterion that has no automated evaluator (e.g. PP5)."""
    reg = _registry(SupplementMode.MERGE)
    results = _tool_results()
    reg._apply_supplement_override(
        results,
        [SupplementEntry(variant_id="x", criterion=ACMGCriterion.PP5,
                         strength=CriterionStrength.SUPPORTING, evidence="manual PP5")],
    )
    by = _by_crit(results)
    assert ACMGCriterion.PP5 in by
    assert by[ACMGCriterion.PP5].triggered


def test_merge_with_no_supplement_is_noop():
    reg = _registry(SupplementMode.MERGE)
    results = _tool_results()
    before = [(r.criterion, r.strength, r.triggered) for r in results]
    reg._apply_supplement_override(results, None)
    after = [(r.criterion, r.strength, r.triggered) for r in results]
    assert before == after


def test_manual_only_with_no_supplement_falls_back_to_tool():
    """A variant absent from the supplement keeps the tool's automated calls,
    even in manual-only mode (per-variant fallback)."""
    reg = _registry(SupplementMode.MANUAL_ONLY)
    results = _tool_results()
    before = [(r.criterion, r.strength, r.triggered) for r in results]
    reg._apply_supplement_override(results, None)
    after = [(r.criterion, r.strength, r.triggered) for r in results]
    assert before == after
    # Tool calls are preserved.
    assert all(r.triggered for r in results)
