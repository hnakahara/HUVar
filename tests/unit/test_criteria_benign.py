"""Unit tests for benign criteria evaluators using mock annotations."""
from acmg_classifier.models.annotation import (
    AnnotationData, GnomADData, RepeatMaskerRegion,
    ConsequenceInfo, SpliceScore, AlphaMissenseData, ESM1bData,
)
from acmg_classifier.models.enums import (
    Assembly, ACMGCriterion, ConsequenceType, CriterionStrength, InSilicoTool,
)
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


class TestBS1:
    def setup_method(self):
        from unittest.mock import MagicMock
        cfg = MagicMock()
        # No per-gene overrides → the conservative default (0.5%) applies.
        cfg.disease_prevalence_tsv.exists.return_value = False
        self.cfg = cfg
        from acmg_classifier.criteria.benign.bs1 import BS1Evaluator
        self.evaluator = BS1Evaluator(cfg)

    def test_bs1_triggered_high_faf(self):
        gd = GnomADData(faf95_popmax=0.01, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered
        assert r.criterion == ACMGCriterion.BS1

    def test_bs1_not_triggered_rare(self):
        gd = GnomADData(faf95_popmax=0.001, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def test_bs1_faf95_zero_does_not_fall_back_to_raw_af(self):
        """FAF95=0.0 (a real value, not "missing") must be used as-is.

        Regression: the old `faf95_popmax or popmax_af or af` chain treated a
        FAF95 lower bound of 0.0 as falsy and silently fell back to the raw
        grpmax/global AF, over-firing BS1 on wide-CI variants. The conservative
        FAF95 of 0.0 is below threshold, so BS1 must NOT trigger here."""
        gd = GnomADData(faf95_popmax=0.0, popmax_af=0.006, af=0.006, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def test_bs1_faf95_missing_falls_back_to_popmax(self):
        """When FAF95 is genuinely absent (None), fall back to raw AF — same
        None-aware behaviour as BA1, so very-rare records still get evaluated."""
        gd = GnomADData(faf95_popmax=None, popmax_af=0.01, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_bs1_no_gnomad(self):
        ann = _annotation(consequences=[_consequence()])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered


class TestMalesAlleleFrequency:
    """X-linked genes whose VCEP states the cutoff "in males" compare against
    gnomAD AF_XY (RPGR/RS1/ABCD1/SLC6A8/OTC), falling back to overall FAF when
    AF_XY is unavailable (gnomAD DB predating the af_xy column)."""

    def _cfg(self, tmp_path):
        from unittest.mock import MagicMock
        cfg = MagicMock()
        tsv = tmp_path / "disease_prevalence.tsv"
        tsv.write_text(
            "gene_symbol\tbs1_threshold\tba1_threshold\taf_basis\n"
            "RPGR\t0.000083\t0.05\tmales\n",
            encoding="utf-8",
        )
        cfg.disease_prevalence_tsv = tsv
        return cfg

    def test_ba1_uses_male_af_over_overall_faf(self, tmp_path):
        from acmg_classifier.criteria.benign.ba1 import BA1Evaluator
        ev = BA1Evaluator(self._cfg(tmp_path))
        # Overall FAF (0.04) < 5%, but male AF (0.06) ≥ 5% → BA1 fires on males.
        gd = GnomADData(faf95_popmax=0.04, af_xy=0.06, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence(gene="RPGR")])
        r = ev.evaluate(_snv(), ann)
        assert r.triggered
        assert "AF_XY (males)" in r.evidence

    def test_ba1_falls_back_to_overall_when_af_xy_missing(self, tmp_path):
        from acmg_classifier.criteria.benign.ba1 import BA1Evaluator
        ev = BA1Evaluator(self._cfg(tmp_path))
        # af_xy None (old DB) → use overall FAF 0.06 ≥ 5% → still fires, but the
        # evidence reflects the overall-FAF fallback, not AF_XY.
        gd = GnomADData(faf95_popmax=0.06, af_xy=None, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence(gene="RPGR")])
        r = ev.evaluate(_snv(), ann)
        assert r.triggered
        assert "AF_XY" not in r.evidence

    def test_bs1_uses_male_af(self, tmp_path):
        from acmg_classifier.criteria.benign.bs1 import BS1Evaluator
        ev = BS1Evaluator(self._cfg(tmp_path))
        # Overall FAF below BS1 (8.3e-5) but male AF above it → BS1 fires.
        gd = GnomADData(faf95_popmax=0.00001, af_xy=0.0002, filter_pass=True)
        ann = _annotation(gnomad=gd, consequences=[_consequence(gene="RPGR")])
        r = ev.evaluate(_snv(), ann)
        assert r.triggered
        assert "AF_XY (males)" in r.evidence


_BS2_TSV = (
    "gene_symbol\tbs2\tinheritance\n"
    "PTEN\tapplicable\tAD\n"
    "MYH7\tnot_applicable\tAD\n"
    "ABCA4\tapplicable\tAR\n"
    "RPGR\tapplicable\tXL\n"
)


class TestBS2:
    def _cfg(self, tmp_path, write_tsv=False):
        from unittest.mock import MagicMock
        cfg = MagicMock()
        p = tmp_path / "disease_prevalence.tsv"
        if write_tsv:
            p.write_text(_BS2_TSV, encoding="utf-8")
        cfg.disease_prevalence_tsv = p   # absent -> empty loader (mode-agnostic)
        cfg.bs2_min_homalt = 5
        cfg.bs2_min_hemi = 5
        cfg.bs2_min_het = 5
        return cfg

    def _ev(self, cfg):
        from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
        return BS2Evaluator(cfg)

    def _ann(self, gene=None, **gd):
        kw = {"gnomad": GnomADData(filter_pass=True, **gd)}
        if gene:
            kw["consequences"] = [_consequence(gene=gene)]
        return _annotation(**kw)

    # --- mode-agnostic fallback (no VCEP row) ---
    def test_homozygotes_fallback_triggers(self, tmp_path):
        r = self._ev(self._cfg(tmp_path)).evaluate(_snv(), self._ann(nhomalt=10))
        assert r.triggered

    def test_low_nhomalt_not_met(self, tmp_path):
        r = self._ev(self._cfg(tmp_path)).evaluate(_snv(), self._ann(nhomalt=2))
        assert not r.triggered

    # --- dominant rule (heterozygous carriers) ---
    def test_dominant_carriers_trigger(self, tmp_path):
        ev = self._ev(self._cfg(tmp_path, write_tsv=True))
        r = ev.evaluate(_snv(), self._ann(gene="PTEN", ac=6, nhomalt=0))  # 6 carriers
        assert r.triggered
        assert "dominant" in r.evidence

    def test_dominant_below_threshold(self, tmp_path):
        ev = self._ev(self._cfg(tmp_path, write_tsv=True))
        r = ev.evaluate(_snv(), self._ann(gene="PTEN", ac=3, nhomalt=0))
        assert not r.triggered

    # --- VCEP applicability gate ---
    def test_not_applicable_blocks(self, tmp_path):
        # A VCEP that bars gnomAD population data blocks the gnomAD-count BS2
        # path; with no ClinVar expert-panel BS2 for the variant, BS2 stays
        # not-met (the fallback only fires on a >=3-star BS2 citation).
        ev = self._ev(self._cfg(tmp_path, write_tsv=True))
        r = ev.evaluate(_snv(), self._ann(gene="MYH7", ac=50, nhomalt=20))
        assert not r.triggered
        assert "bars gnomAD-based BS2" in r.evidence

    # --- recessive ignores heterozygotes ---
    def test_recessive_uses_homozygotes_not_carriers(self, tmp_path):
        ev = self._ev(self._cfg(tmp_path, write_tsv=True))
        # Many heterozygous carriers but few homozygotes -> AR must NOT fire.
        r = ev.evaluate(_snv(), self._ann(gene="ABCA4", ac=100, nhomalt=2))
        assert not r.triggered
        r2 = ev.evaluate(_snv(), self._ann(gene="ABCA4", ac=20, nhomalt=8))
        assert r2.triggered

    def test_xlinked_uses_hemizygotes(self, tmp_path):
        ev = self._ev(self._cfg(tmp_path, write_tsv=True))
        r = ev.evaluate(_snv(), self._ann(gene="RPGR", nhemi=7, ac=7))
        assert r.triggered and "X-linked" in r.evidence


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


class TestBP4AlphaMissense:
    def setup_method(self):
        from unittest.mock import MagicMock
        self.cfg = MagicMock()
        from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
        self.evaluator = BP4Evaluator(self.cfg)

    def _ann(self, score):
        return _annotation(
            alphamissense=AlphaMissenseData(score=score),
            consequences=[_consequence(ConsequenceType.MISSENSE)],
        )

    def test_three_point_at_0_07(self):
        # Bergquist 2024: AlphaMissense BP4 ≤ 0.070 is ThreePoint, NOT Strong.
        r = self.evaluator.evaluate(_snv(), self._ann(0.070))
        assert r.triggered
        assert r.strength == CriterionStrength.THREE_POINT

    def test_moderate_at_0_099(self):
        r = self.evaluator.evaluate(_snv(), self._ann(0.099))
        assert r.triggered
        assert r.strength == CriterionStrength.MODERATE

    def test_supporting_at_0_169(self):
        r = self.evaluator.evaluate(_snv(), self._ann(0.169))
        assert r.triggered
        assert r.strength == CriterionStrength.SUPPORTING

    def test_indeterminate_at_0_5(self):
        r = self.evaluator.evaluate(_snv(), self._ann(0.5))
        assert not r.triggered


class TestBP4ESM1b:
    def setup_method(self):
        from unittest.mock import MagicMock
        cfg = MagicMock()
        cfg.insilico_tool = InSilicoTool.ESM1B
        from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
        self.evaluator = BP4Evaluator(cfg)

    def _ann(self, llr):
        return _annotation(
            esm1b=ESM1bData(llr=llr),
            consequences=[_consequence(ConsequenceType.MISSENSE)],
        )

    def test_three_point_at_8_8(self):
        # ESM1b BP4 ≥ 8.8 is ThreePoint (no Strong category).
        r = self.evaluator.evaluate(_snv(), self._ann(8.8))
        assert r.triggered
        assert r.strength == CriterionStrength.THREE_POINT

    def test_moderate_at_minus3_2(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-3.2))
        assert r.triggered
        assert r.strength == CriterionStrength.MODERATE

    def test_supporting_at_minus6_3(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-6.3))
        assert r.triggered
        assert r.strength == CriterionStrength.SUPPORTING

    def test_indeterminate_at_minus8(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-8.0))
        assert not r.triggered

    def test_pathogenic_range_not_triggered(self):
        r = self.evaluator.evaluate(_snv(), self._ann(-15.0))
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

    def test_bp7_deep_intronic_without_splice_not_met(self):
        # Walker 2023: deep-intronic BP7 still requires a splice prediction of
        # no impact. With no splice predictor, distance alone must NOT fire.
        c = _consequence(ConsequenceType.INTRON)
        c = c.model_copy(update={"intron_distance_from_splice": 15})
        ann = _annotation(consequences=[c], splice=None)
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def test_bp7_deep_intronic_with_splice_impact_not_met(self):
        # Deep-intronic but the splice tool predicts impact → BP7 withheld.
        c = _consequence(ConsequenceType.INTRON)
        c = c.model_copy(update={"intron_distance_from_splice": 15})
        sp = SpliceScore(tool="openspliceai", is_available=True, max_delta=0.8)
        ann = _annotation(consequences=[c], splice=sp)
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered

    def _phylop_stub(self, score):
        class _Stub:
            def is_available(self_inner):
                return True

            def value(self_inner, chrom, pos):
                return score
        return _Stub()

    def test_bp7_synonymous_highly_conserved_blocked(self):
        # Splice benign but the nucleotide is highly conserved → BP7 withheld.
        self.cfg.bp7_phylop_max = 2.0
        self.evaluator._phylop = self._phylop_stub(7.5)
        sp = SpliceScore(tool="openspliceai", is_available=True, max_delta=0.02)
        ann = _annotation(consequences=[_consequence(ConsequenceType.SYNONYMOUS)], splice=sp)
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered
        assert "conserved" in r.evidence.lower()

    def test_bp7_synonymous_not_conserved_fires(self):
        # Splice benign and NOT highly conserved → BP7 fires.
        self.cfg.bp7_phylop_max = 2.0
        self.evaluator._phylop = self._phylop_stub(0.3)
        sp = SpliceScore(tool="openspliceai", is_available=True, max_delta=0.02)
        ann = _annotation(consequences=[_consequence(ConsequenceType.SYNONYMOUS)], splice=sp)
        r = self.evaluator.evaluate(_snv(), ann)
        assert r.triggered

    def test_bp7_missense_not_triggered(self):
        ann = _annotation(consequences=[_consequence(ConsequenceType.MISSENSE)])
        r = self.evaluator.evaluate(_snv(), ann)
        assert not r.triggered
