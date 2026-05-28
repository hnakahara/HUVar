"""Unit tests for ACMG 2015 and Bayesian classifiers."""
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength, Pathogenicity
from acmg_classifier.classification.classifier_2015 import Classifier2015
from acmg_classifier.classification.classifier_bayesian import ClassifierBayesian


def _pvs1():
    return CriteriaResult.met(ACMGCriterion.PVS1, CriterionStrength.VERY_STRONG)

def _ps1():
    return CriteriaResult.met(ACMGCriterion.PS1)

def _pm2():
    return CriteriaResult.met(ACMGCriterion.PM2, CriterionStrength.SUPPORTING)

def _ba1():
    return CriteriaResult.met(ACMGCriterion.BA1)

def _bs1():
    return CriteriaResult.met(ACMGCriterion.BS1)

def _bp4():
    return CriteriaResult.met(ACMGCriterion.BP4)


def _not_met(criterion):
    return CriteriaResult.not_met(criterion)


clf_2015 = Classifier2015()
clf_bay = ClassifierBayesian()


class TestClassifier2015:
    def test_ba1_benign(self):
        results = [_ba1()]
        path, _ = clf_2015.classify(results)
        assert path == Pathogenicity.BENIGN

    def test_pvs1_pm2_likely_pathogenic(self):
        results = [_pvs1(), _pm2()]
        path, _ = clf_2015.classify(results)
        assert path == Pathogenicity.LIKELY_PATHOGENIC

    def test_pvs1_ps1_pathogenic(self):
        results = [_pvs1(), _ps1()]
        path, _ = clf_2015.classify(results)
        assert path == Pathogenicity.PATHOGENIC

    def test_two_strong_pathogenic(self):
        ps1 = CriteriaResult.met(ACMGCriterion.PS1)
        ps2 = CriteriaResult.met(ACMGCriterion.PS2)
        path, _ = clf_2015.classify([ps1, ps2])
        assert path == Pathogenicity.PATHOGENIC

    def test_no_criteria_vus(self):
        results = [_not_met(c) for c in ACMGCriterion]
        path, _ = clf_2015.classify(results)
        assert path == Pathogenicity.VUS

    def test_bs1_bp4_likely_benign(self):
        results = [_bs1(), _bp4()]
        path, _ = clf_2015.classify(results)
        assert path == Pathogenicity.BENIGN  # 1 Strong + 1 Supporting = Benign (ACMG 2015 Table 5)


class TestClassifierBayesian:
    def test_pvs1_score_is_8(self):
        score, _ = clf_bay.classify([_pvs1()])
        assert score == 8

    def test_ba1_score_is_minus8(self):
        score, _ = clf_bay.classify([_ba1()])
        assert score == -8

    def test_pvs1_pm2_likely_pathogenic(self):
        score, path = clf_bay.classify([_pvs1(), _pm2()])
        assert score == 9
        assert path == Pathogenicity.LIKELY_PATHOGENIC

    def test_pvs1_ps1_pathogenic(self):
        score, path = clf_bay.classify([_pvs1(), _ps1()])
        assert score == 12
        assert path == Pathogenicity.PATHOGENIC

    def test_zero_vus(self):
        _, path = clf_bay.classify([])
        assert path == Pathogenicity.VUS

    def test_ba1_benign(self):
        score, path = clf_bay.classify([_ba1()])
        assert score == -8
        assert path == Pathogenicity.BENIGN
