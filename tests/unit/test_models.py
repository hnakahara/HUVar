"""Unit tests for Pydantic models."""
import pytest
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.enums import Assembly, VariantType, CriterionStrength, CriterionDirection
from acmg_classifier.models.criteria import CriteriaResult, STRENGTH_POINTS
from acmg_classifier.models.enums import ACMGCriterion


def test_variant_chrom_normalisation():
    v = VariantRecord(chrom="17", pos=100, ref="G", alt="A", assembly=Assembly.GRCH38)
    assert v.chrom == "chr17"


def test_variant_key():
    v = VariantRecord(chrom="chr17", pos=43044295, ref="G", alt="A", assembly=Assembly.GRCH38)
    assert v.key == "chr17:43044295:G:A"


def test_variant_type_snv():
    v = VariantRecord(chrom="chr1", pos=100, ref="G", alt="A", assembly=Assembly.GRCH38)
    assert v.variant_type == VariantType.SNV
    assert v.is_snv


def test_variant_type_indel():
    v = VariantRecord(chrom="chr1", pos=100, ref="GT", alt="G", assembly=Assembly.GRCH38)
    assert v.variant_type == VariantType.INDEL
    assert v.is_indel


def test_criteria_result_points_pathogenic():
    r = CriteriaResult.met(ACMGCriterion.PVS1)
    assert r.points == 8


def test_criteria_result_points_benign():
    r = CriteriaResult.met(ACMGCriterion.BA1)
    assert r.points == -8


def test_criteria_result_not_met_zero_points():
    r = CriteriaResult.not_met(ACMGCriterion.PM2)
    assert r.points == 0
    assert not r.triggered


def test_criteria_result_suppressed_zero_points():
    r = CriteriaResult.met(ACMGCriterion.PP3)
    r.suppressed = True
    assert r.points == 0


def test_strength_points_completeness():
    for criterion in (ACMGCriterion.PVS1, ACMGCriterion.PS1, ACMGCriterion.PM2,
                      ACMGCriterion.BA1, ACMGCriterion.BS1, ACMGCriterion.BP4):
        result = CriteriaResult.met(criterion)
        assert result.points != 0
