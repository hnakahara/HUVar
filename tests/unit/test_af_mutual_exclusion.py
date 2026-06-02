"""Unit tests for BA1 / BS1 / PM2 mutual exclusivity (ClinGen SVI).

A variant may carry at most one allele-frequency criterion. The priority
ladder is BA1 > BS1 > PM2; lower-priority triggered criteria are suppressed
(kept in the audit trail, 0 points) rather than dropped.
"""
from acmg_classifier.criteria.registry import _apply_af_mutual_exclusion
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion


def _met(crit):
    return CriteriaResult.met(crit, evidence="x")


def _not_met(crit):
    return CriteriaResult.not_met(crit, evidence="x")


def _by(results, crit):
    return next(r for r in results if r.criterion == crit)


def test_ba1_suppresses_bs1_and_pm2():
    results = [_met(ACMGCriterion.BA1), _met(ACMGCriterion.BS1), _met(ACMGCriterion.PM2)]
    _apply_af_mutual_exclusion(results)
    assert not _by(results, ACMGCriterion.BA1).suppressed
    assert _by(results, ACMGCriterion.BS1).suppressed
    assert _by(results, ACMGCriterion.PM2).suppressed
    assert "[suppressed: BA1 active]" in _by(results, ACMGCriterion.BS1).evidence


def test_bs1_suppresses_pm2_when_no_ba1():
    results = [_not_met(ACMGCriterion.BA1), _met(ACMGCriterion.BS1), _met(ACMGCriterion.PM2)]
    _apply_af_mutual_exclusion(results)
    assert not _by(results, ACMGCriterion.BS1).suppressed
    assert _by(results, ACMGCriterion.PM2).suppressed
    assert "[suppressed: BS1 active]" in _by(results, ACMGCriterion.PM2).evidence


def test_pm2_alone_is_kept():
    results = [_not_met(ACMGCriterion.BA1), _not_met(ACMGCriterion.BS1), _met(ACMGCriterion.PM2)]
    _apply_af_mutual_exclusion(results)
    assert not _by(results, ACMGCriterion.PM2).suppressed


def test_no_af_criteria_triggered_is_noop():
    results = [_not_met(ACMGCriterion.BA1), _not_met(ACMGCriterion.BS1), _not_met(ACMGCriterion.PM2)]
    _apply_af_mutual_exclusion(results)
    assert all(not r.suppressed for r in results)


def test_ba1_suppresses_pm2_when_bs1_absent():
    # BS1 not triggered must not break selecting PM2 as a loser under BA1.
    results = [_met(ACMGCriterion.BA1), _not_met(ACMGCriterion.BS1), _met(ACMGCriterion.PM2)]
    _apply_af_mutual_exclusion(results)
    assert not _by(results, ACMGCriterion.BA1).suppressed
    assert _by(results, ACMGCriterion.PM2).suppressed


def test_already_suppressed_winner_is_skipped():
    # If BA1 was already suppressed upstream, BS1 should win over PM2.
    ba1 = _met(ACMGCriterion.BA1)
    ba1.suppressed = True
    results = [ba1, _met(ACMGCriterion.BS1), _met(ACMGCriterion.PM2)]
    _apply_af_mutual_exclusion(results)
    assert not _by(results, ACMGCriterion.BS1).suppressed
    assert _by(results, ACMGCriterion.PM2).suppressed


def test_suppressed_criteria_contribute_zero_points():
    results = [_met(ACMGCriterion.BA1), _met(ACMGCriterion.BS1), _met(ACMGCriterion.PM2)]
    _apply_af_mutual_exclusion(results)
    # BS1 (Strong benign, -4) and PM2 (Supporting path, +1) must zero out.
    assert _by(results, ACMGCriterion.BS1).points == 0
    assert _by(results, ACMGCriterion.PM2).points == 0
    assert _by(results, ACMGCriterion.BA1).points != 0
