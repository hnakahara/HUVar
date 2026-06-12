"""Unit tests for pathogenic criteria evaluators."""
from unittest.mock import MagicMock
from acmg_classifier.models.annotation import (
    AnnotationData, GnomADData, AlphaMissenseData, ESM1bData, ClinVarRecord, ConsequenceInfo
)
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, CriterionStrength, InSilicoTool,
)
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

    def test_pm2_uses_raw_af_not_faf95(self):
        """PM2 judges on the RAW grpmax AF (ClinGen Hearing Loss VCEP), not the
        FAF95 lower bound. Raw AF above the dominant cutoff → NOT met, even
        though the conservative FAF95 is tiny."""
        gd = GnomADData(ac=10, af=0.0005, popmax_af=0.0005,
                        faf95_popmax=0.00001, filter_pass=True)
        ann = AnnotationData(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def test_pm2_faf95_zero_not_absent_when_raw_af_present(self):
        """A FAF95 of 0.0 must NOT be read as "absent" when the variant is
        actually observed (raw AF 0.0005). PM2 uses the raw AF → NOT met."""
        gd = GnomADData(ac=10, af=0.0005, popmax_af=0.0005,
                        faf95_popmax=0.0, filter_pass=True)
        ann = AnnotationData(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def test_pm2_raw_af_below_threshold_met(self):
        gd = GnomADData(ac=5, af=0.00005, popmax_af=0.00005, filter_pass=True)
        ann = AnnotationData(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered
        assert r.strength == CriterionStrength.SUPPORTING

    def test_pm2_popmax_af_none_falls_back_to_global_af(self):
        gd = GnomADData(ac=10, af=0.0005, popmax_af=None, filter_pass=True)
        ann = AnnotationData(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered  # raw global AF 0.0005 >= dominant cutoff

    def test_pm2_filter_failed_ac0_still_met(self):
        # eRepo: a filter-failed (AC0) record is effectively absent → PM2 applies
        # (was a blanket false negative when filter failure blocked PM2).
        gd = GnomADData(ac=0, af=0.0, popmax_af=0.0, filter_pass=False)
        ann = AnnotationData(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered
        assert "filter" in r.evidence.lower()

    def test_pm2_filter_failed_rare_still_met(self):
        # Filter-failed but extremely rare (AF 1e-6 < dominant 1e-4) → PM2 applies.
        gd = GnomADData(ac=1, af=1e-6, popmax_af=1e-6, filter_pass=False)
        ann = AnnotationData(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_pm2_filter_failed_common_not_met(self):
        # A genuinely common filter-failed call (e.g. inflated segdup AF) still
        # fails the threshold → no false positive introduced by the fix.
        gd = GnomADData(ac=500, af=0.01, popmax_af=0.01, filter_pass=False)
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


class TestPP3ESM1b:
    def setup_method(self):
        cfg = MagicMock()
        cfg.insilico_tool = InSilicoTool.ESM1B
        from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
        self.evaluator = PP3Evaluator(cfg)

    def _ann(self, llr):
        return AnnotationData(
            esm1b=ESM1bData(llr=llr),
            consequences=[_consequence(ConsequenceType.MISSENSE)],
        )

    def test_pp3_strong_at_minus24(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-24.0))
        assert r.triggered
        assert r.strength == CriterionStrength.STRONG

    def test_pp3_three_point_at_minus14(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-14.0))
        assert r.triggered
        assert r.strength == CriterionStrength.THREE_POINT

    def test_pp3_moderate_at_minus12_2(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-12.2))
        assert r.triggered
        assert r.strength == CriterionStrength.MODERATE

    def test_pp3_supporting_at_minus10_7(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-10.7))
        assert r.triggered
        assert r.strength == CriterionStrength.SUPPORTING

    def test_pp3_indeterminate_at_minus5(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-5.0))
        assert not r.triggered

    def test_pp3_benign_range_not_triggered(self):
        r = self.evaluator.evaluate(_snv(), self._ann(10.0))
        assert not r.triggered

    def test_pp3_falls_through_when_no_esm1b(self):
        # Tool is ESM1B but score missing → "No in-silico score available"
        ann = AnnotationData(
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
