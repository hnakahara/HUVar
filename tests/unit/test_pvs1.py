"""Unit tests for PVS1 decision tree."""
import pytest
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pvs1.nmd_predictor import predicts_nmd, is_last_exon, is_penultimate_exon
from acmg_classifier.pvs1.transcript_evaluator import gene_has_lof_mechanism


def _consequence(ctype, exon="5/24"):
    return ConsequenceInfo(
        transcript_id="NM_007294.4",
        gene_id="ENSG00000012048",
        gene_symbol="BRCA1",
        consequence=ctype,
        biotype="protein_coding",
        is_mane_select=True,
        exon=exon,
    )


class TestNMDPredictor:
    def test_nmd_predicted_early_exon(self):
        c = _consequence(ConsequenceType.STOP_GAINED, exon="5/24")
        assert predicts_nmd(c) is True

    def test_nmd_not_last_exon(self):
        c = _consequence(ConsequenceType.STOP_GAINED, exon="24/24")
        assert predicts_nmd(c) is False

    def test_single_exon(self):
        c = _consequence(ConsequenceType.STOP_GAINED, exon="1/1")
        assert predicts_nmd(c) is False

    def test_is_last_exon(self):
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="12/12")
        assert is_last_exon(c) is True

    def test_is_penultimate_exon(self):
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="11/12")
        assert is_penultimate_exon(c) is True


class TestGeneLoFMechanism:
    def test_lof_intolerant_low_loeuf(self):
        assert gene_has_lof_mechanism(None, gnomad_loeuf=0.10) is True

    def test_lof_tolerant_high_loeuf(self):
        assert gene_has_lof_mechanism(None, gnomad_loeuf=0.80) is False

    def test_no_loeuf_defaults_to_true(self):
        assert gene_has_lof_mechanism(None, gnomad_loeuf=None) is True


class TestPVS1DecisionTree:
    def setup_method(self):
        from unittest.mock import MagicMock
        self.cfg = MagicMock()

    def test_frameshift_nmd_no_rescue_very_strong(self):
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="5/24")
        gd = GnomADData(loeuf=0.10)
        ann = AnnotationData(consequences=[c], gnomad=gd)
        v = VariantRecord(chrom="chr17", pos=100, ref="G", alt="GA", assembly=Assembly.GRCH38)
        strength, evidence = evaluate_pvs1(v, ann, self.cfg)
        assert strength == CriterionStrength.VERY_STRONG

    def test_start_loss_moderate(self):
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        c = _consequence(ConsequenceType.START_LOST, exon="1/24")
        ann = AnnotationData(consequences=[c])
        v = VariantRecord(chrom="chr17", pos=100, ref="G", alt="A", assembly=Assembly.GRCH38)
        strength, evidence = evaluate_pvs1(v, ann, self.cfg)
        assert strength == CriterionStrength.MODERATE

    def test_last_exon_frameshift_no_domain_moderate(self):
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        c = _consequence(ConsequenceType.FRAMESHIFT, exon="24/24")
        gd = GnomADData(loeuf=0.10)
        ann = AnnotationData(consequences=[c], gnomad=gd)
        v = VariantRecord(chrom="chr17", pos=100, ref="GA", alt="G", assembly=Assembly.GRCH38)
        strength, evidence = evaluate_pvs1(v, ann, self.cfg)
        assert strength == CriterionStrength.MODERATE
