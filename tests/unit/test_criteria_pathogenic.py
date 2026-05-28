"""Unit tests for pathogenic criteria evaluators."""
import pytest
from unittest.mock import MagicMock, patch
from acmg_classifier.models.annotation import (
    AnnotationData, GnomADData, AlphaMissenseData, ClinVarRecord, ConsequenceInfo
)
from acmg_classifier.models.enums import Assembly, ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord


def _snv():
    return VariantRecord(chrom="chr17", pos=43044295, ref="G", alt="A", assembly=Assembly.GRCH38)


def _consequence(ctype=ConsequenceType.MISSENSE, gene="BRCA1", exon="5/24",
                 protein_pos=1699, aa_change="R1699Q", codon_pos=1699):
    return ConsequenceInfo(
        transcript_id="NM_007294.4",
        gene_id="ENSG00000012048",
        gene_symbol=gene,
        consequence=ctype,
        biotype="protein_coding",
        is_mane_select=True,
        exon=exon,
        protein_position=protein_pos,
        amino_acid_change=aa_change,
        codon_position=codon_pos,
    )


class TestPM2:
    def setup_method(self):
        cfg = MagicMock()
        from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
        self.evaluator = PM2Evaluator(cfg)

    def test_pm2_absent_no_record(self):
        ann = AnnotationData(gnomad=None, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered
        assert r.strength == CriterionStrength.SUPPORTING

    def test_pm2_absent_ac_zero(self):
        gd = GnomADData(ac=0, af=0.0, faf95_popmax=0.0, filter_pass=True)
        ann = AnnotationData(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_pm2_common_not_met(self):
        gd = GnomADData(ac=50, af=0.001, faf95_popmax=0.002, filter_pass=True)
        ann = AnnotationData(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered


class TestPP3:
    def setup_method(self):
        cfg = MagicMock()
        from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
        self.evaluator = PP3Evaluator(cfg)

    def test_pp3_strong_alphamissense(self):
        am = AlphaMissenseData(score=0.995)
        ann = AnnotationData(
            alphamissense=am,
            consequences=[_consequence(ConsequenceType.MISSENSE)],
        )
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered
        assert r.strength == CriterionStrength.STRONG

    def test_pp3_moderate_alphamissense(self):
        am = AlphaMissenseData(score=0.950)
        ann = AnnotationData(
            alphamissense=am,
            consequences=[_consequence(ConsequenceType.MISSENSE)],
        )
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered
        assert r.strength == CriterionStrength.MODERATE

    def test_pp3_indeterminate_not_triggered(self):
        am = AlphaMissenseData(score=0.500)
        ann = AnnotationData(
            alphamissense=am,
            consequences=[_consequence(ConsequenceType.MISSENSE)],
        )
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def test_pp3_benign_score_not_triggered(self):
        am = AlphaMissenseData(score=0.050)
        ann = AnnotationData(
            alphamissense=am,
            consequences=[_consequence(ConsequenceType.MISSENSE)],
        )
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered


class TestPP5:
    def setup_method(self):
        cfg = MagicMock()
        from acmg_classifier.criteria.pathogenic.pp5 import PP5Evaluator
        self.evaluator = PP5Evaluator(cfg)

    def test_pp5_triggered_expert_panel(self):
        cv = ClinVarRecord(
            variation_id="123456",
            clinical_significance="Pathogenic",
            review_status="reviewed by expert panel",
            star_rating=3,
        )
        ann = AnnotationData(clinvar_vcf=[cv])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_pp5_not_triggered_one_star(self):
        cv = ClinVarRecord(
            variation_id="789",
            clinical_significance="Pathogenic",
            review_status="criteria provided, single submitter",
            star_rating=1,
        )
        ann = AnnotationData(clinvar_vcf=[cv])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def test_pp5_not_triggered_vus(self):
        cv = ClinVarRecord(
            variation_id="999",
            clinical_significance="Uncertain significance",
            review_status="reviewed by expert panel",
            star_rating=3,
        )
        ann = AnnotationData(clinvar_vcf=[cv])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered
