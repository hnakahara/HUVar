"""Unit tests for benign criteria evaluators using mock annotations."""
from acmg_classifier.models.annotation import AnnotationData, GnomADData, RepeatMaskerRegion
from acmg_classifier.models.annotation import ConsequenceInfo, SpliceScore
from acmg_classifier.models.enums import Assembly, ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord


def _snv(chrom="chr17", pos=100, ref="G", alt="A"):
    return VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=Assembly.GRCH38)


def _consequence(ctype=ConsequenceType.MISSENSE, exon="5/12", gene="BRCA1"):
    return ConsequenceInfo(
        transcript_id="NM_007294.4",
        gene_id="ENSG00000012048",
        gene_symbol=gene,
        consequence=ctype,
        biotype="protein_coding",
        is_mane_select=True,
        exon=exon,
    )


def _annotation(**kwargs):
    return AnnotationData(**kwargs)


class TestBA1:
    def setup_method(self):
        from unittest.mock import MagicMock
        cfg = MagicMock()
        cfg.clinvar_sqlite = MagicMock()
        cfg.clinvar_sqlite.exists.return_value = False
        self.cfg = cfg
        from acmg_classifier.criteria.benign.ba1 import BA1Evaluator
        self.evaluator = BA1Evaluator(cfg)

    def test_ba1_triggered_high_faf(self):
        gd = GnomADData(faf95_popmax=0.06, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered
        assert r.criterion == ACMGCriterion.BA1

    def test_ba1_not_triggered_rare(self):
        gd = GnomADData(faf95_popmax=0.001, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def test_ba1_no_gnomad(self):
        ann = _annotation(consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered


class TestBS2:
    def setup_method(self):
        from unittest.mock import MagicMock
        self.cfg = MagicMock()
        from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
        self.evaluator = BS2Evaluator(self.cfg)

    def test_bs2_triggered_homozygotes(self):
        gd = GnomADData(nhomalt=10, filter_pass=True)
        ann = _annotation(gnomad=gd)
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_bs2_not_triggered_low_nhomalt(self):
        gd = GnomADData(nhomalt=2, filter_pass=True)
        ann = _annotation(gnomad=gd)
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered


class TestBP3:
    def setup_method(self):
        from unittest.mock import MagicMock
        self.cfg = MagicMock()
        from acmg_classifier.criteria.benign.bp3 import BP3Evaluator
        self.evaluator = BP3Evaluator(self.cfg)

    def test_bp3_in_frame_in_repeat(self):
        rep = RepeatMaskerRegion(in_repeat=True, repeat_class="SINE", repeat_name="Alu")
        ann = _annotation(
            consequences=[_consequence(ConsequenceType.INFRAME_DELETION)],
            repeat=rep,
        )
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_bp3_not_in_repeat(self):
        rep = RepeatMaskerRegion(in_repeat=False)
        ann = _annotation(
            consequences=[_consequence(ConsequenceType.INFRAME_DELETION)],
            repeat=rep,
        )
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered


class TestBP7:
    def setup_method(self):
        from unittest.mock import MagicMock
        self.cfg = MagicMock()
        from acmg_classifier.criteria.benign.bp7 import BP7Evaluator
        self.evaluator = BP7Evaluator(self.cfg)

    def test_bp7_synonymous_with_benign_splice(self):
        sp = SpliceScore(tool="squirls", is_available=True, raw_score=0.05)
        ann = _annotation(
            consequences=[_consequence(ConsequenceType.SYNONYMOUS)],
            splice=sp,
        )
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_bp7_deep_intronic(self):
        c = _consequence(ConsequenceType.INTRON)
        c = c.model_copy(update={"intron_distance_from_splice": 15})
        sp = SpliceScore(tool="squirls", is_available=True, raw_score=0.05)
        ann = _annotation(consequences=[c], splice=sp)
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_bp7_missense_not_triggered(self):
        ann = _annotation(consequences=[_consequence(ConsequenceType.MISSENSE)])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered
