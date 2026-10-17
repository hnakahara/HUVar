"""VCEP rule gaps found by the eRepo audit (Genome Medicine revision):
PP3 splice strength / gene cutoffs, MDEP PM2 copy rule, PM5 conditions
(splice exclusion, PP3 / REVEL requirements, MDEP Grantham downgrade, F8/F9 two
LP comparators, CDH1 canonical splice), MYOC PM4 and the HNF1A promoter PM1."""
from unittest.mock import MagicMock

from acmg_classifier.criteria import pm5_conditions as cond
from acmg_classifier.criteria.pathogenic.pm1 import _hnf1a_promoter_site
from acmg_classifier.criteria.pathogenic.pm2 import _mdep_incomputable_faf
from acmg_classifier.criteria.pathogenic.pm4 import _myoc_pm4
from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
from acmg_classifier.criteria.pathogenic.pp3 import (
    PP3Evaluator, _spliceai_pp3, _splice_cutoff,
)
from acmg_classifier.criteria.pm5_ptc import evaluate_pm5_ptc
from acmg_classifier.criteria.registry import _apply_pm5_requires_pp3
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, GnomADData, SpliceScore,
)
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion, Assembly, ConsequenceType as C, CriterionStrength as S,
)
from acmg_classifier.models.variant import VariantRecord

from tests.unit.test_pm5_grantham import _db, _row


def _pc(gene, cons=C.MISSENSE, pos=175, aa="R/W", hgvs_c=None, hgvs_p=None, tx="NM_1"):
    return ConsequenceInfo(
        transcript_id=tx, gene_id="ENSG", gene_symbol=gene, consequence=cons,
        biotype="protein_coding", is_mane_select=True, protein_position=pos,
        codon_position=pos, amino_acids=aa, hgvs_c=hgvs_c,
        hgvs_p=hgvs_p or "NM:p.Arg175Trp",
    )


def _snv(ref="C", alt="T"):
    return VariantRecord(chrom="chr1", pos=100, ref=ref, alt=alt, assembly=Assembly.GRCH38)


# ------------------------------------------------------------------ PP3 (G1/G2)
class TestPp3Splice:
    def test_splice_is_supporting(self):
        assert _spliceai_pp3(0.95) == S.SUPPORTING
        assert _spliceai_pp3(0.19) is None

    def test_gene_cutoffs(self):
        assert _splice_cutoff("PAH") == (0.5, True)
        assert _splice_cutoff("RUNX1") == (0.38, False)
        assert _splice_cutoff("BRCA1") == (0.2, False)
        assert _spliceai_pp3(0.45, *_splice_cutoff("DYSF")) is None
        assert _spliceai_pp3(0.50, *_splice_cutoff("DYSF")) == S.SUPPORTING
        assert _spliceai_pp3(0.50, *_splice_cutoff("PAH")) is None  # "> 0.5"

    def _ev(self):
        cfg = MagicMock()
        ev = PP3Evaluator.__new__(PP3Evaluator)
        ev._cfg = cfg
        return ev

    def test_missense_keeps_stronger_protein_call(self):
        ev = self._ev()
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_pc("BRCA1")],
                             splice=SpliceScore(tool="spliceai", max_delta=0.6, is_available=True))
        prot = CriteriaResult.met(ACMGCriterion.PP3, S.STRONG, "REVEL=0.95 (Strong)")
        assert ev._with_missense_splice(ann.consequences[0], ann, prot).strength == S.STRONG

    def test_missense_splice_rescues_unmet_protein(self):
        ev = self._ev()
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_pc("BRCA1")],
                             splice=SpliceScore(tool="spliceai", max_delta=0.6, is_available=True))
        prot = CriteriaResult.not_met(ACMGCriterion.PP3, "REVEL=0.30")
        r = ev._with_missense_splice(ann.consequences[0], ann, prot)
        assert r.triggered and r.strength == S.SUPPORTING

    def test_missense_splice_below_gene_cutoff(self):
        ev = self._ev()
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_pc("PAH")],
                             splice=SpliceScore(tool="spliceai", max_delta=0.4, is_available=True))
        prot = CriteriaResult.not_met(ACMGCriterion.PP3, "REVEL=0.30")
        assert not ev._with_missense_splice(ann.consequences[0], ann, prot).triggered


# ------------------------------------------------------------------ PM2 (G3)
class TestMdepPm2:
    def test_single_copy(self):
        ok, _ = _mdep_incomputable_faf(GnomADData(ac=1, ac_grpmax=1, popmax_pop="sas"))
        assert ok

    def test_two_copies_nfe(self):
        ok, _ = _mdep_incomputable_faf(GnomADData(ac=2, ac_grpmax=2, popmax_pop="nfe"))
        assert ok

    def test_two_copies_other_group(self):
        ok, _ = _mdep_incomputable_faf(GnomADData(ac=2, ac_grpmax=2, popmax_pop="afr"))
        assert not ok


# ------------------------------------------------------------------ PM5 (G6/G7/G8)
def _pm5(tmp_path, rows, tsv_rows=""):
    tsv = tmp_path / "disease_prevalence.tsv"
    tsv.write_text("gene_symbol\tpm5_grantham\tpm5_excludes\tpm5_max\tpm5_lp\n" + tsv_rows,
                   encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = _db(tmp_path, rows)
    cfg.pm5_min_stars = 1
    return PM5Evaluator(cfg)


class TestPm5Conditions:
    def test_mdep_grantham_downgrade(self, tmp_path):
        # R->H (29) vs comparator R->W (101): candidate smaller -> Supporting
        ev = _pm5(tmp_path, [_row("1", "GCK", "NM:p.Arg175Trp", "R175W", 175, "Pathogenic")],
                  "GCK\tge\t\tStrong\t\n")
        ann = AnnotationData(gnomad=GnomADData(),
                             consequences=[_pc("GCK", aa="R/H", hgvs_p="NM:p.Arg175His")])
        r = ev.evaluate(_snv(), ann)
        assert r.triggered and r.strength == S.SUPPORTING

    def test_candidate_splice_excluded(self, tmp_path):
        ev = _pm5(tmp_path, [_row("1", "CTLA4", "NM:p.Arg175Gln", "R175Q", 175, "Pathogenic")])
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_pc("CTLA4")],
                             splice=SpliceScore(tool="spliceai", max_delta=0.3, is_available=True))
        r = ev.evaluate(_snv(), ann)
        assert not r.triggered and "splicing" in r.evidence

    def test_dysf_revel_required(self, tmp_path):
        assert cond.REVEL_MIN["DYSF"] == 0.7

    def test_f8_two_lp_moderate(self, tmp_path):
        ev = _pm5(tmp_path, [
            _row("1", "F8", "NM:p.Arg175Gln", "R175Q", 175, "Likely pathogenic"),
            _row("2", "F8", "NM:p.Arg175Gly", "R175G", 175, "Likely pathogenic"),
        ])
        r = ev.evaluate(_snv(), AnnotationData(gnomad=GnomADData(), consequences=[_pc("F8")]))
        assert r.triggered and r.strength == S.MODERATE

    def test_requires_pp3(self):
        pm5 = CriteriaResult.met(ACMGCriterion.PM5, S.MODERATE, "ClinVar same codon")
        pp3 = CriteriaResult.not_met(ACMGCriterion.PP3, "REVEL=0.3")
        ann = AnnotationData(gnomad=GnomADData(), consequences=[_pc("MLH1")])
        _apply_pm5_requires_pp3([pm5, pp3], ann)
        assert pm5.suppressed

    def test_comparator_modes(self):
        assert cond.splice_ok("GUCY2D", 0.15, "comparator", 0.12) is True
        assert cond.splice_ok("GUCY2D", 0.25, "comparator", 0.12) is False
        assert cond.splice_ok("RPE65", 0.50, "comparator", 0.30) is False
        assert cond.splice_ok("PIK3R1", 0.10, "comparator") is True
        assert cond.splice_ok("RPGR", 0.30, "comparator") is None  # candidate-only

    def test_cdh1_canonical_splice(self):
        pc = _pc("CDH1", cons=C.SPLICE_DONOR, pos=None, aa=None,
                 hgvs_c="c.1137+1G>A", tx="NM_004360.5")
        s, ev = evaluate_pm5_ptc(pc, AnnotationData(gnomad=GnomADData(), consequences=[pc]))
        assert s == S.SUPPORTING and ev.startswith("PM5_PTC")


# ------------------------------------------------------------------ PM4 (G9)
class TestMyocPm4:
    def test_truncating_large_moderate(self):
        s, _ = _myoc_pm4(_pc("MYOC", cons=C.STOP_GAINED, pos=272), _snv())
        assert s == S.MODERATE

    def test_truncating_small_supporting(self):
        s, _ = _myoc_pm4(_pc("MYOC", cons=C.FRAMESHIFT, pos=489), _snv())
        assert s == S.SUPPORTING

    def test_outside_domain(self):
        assert _myoc_pm4(_pc("MYOC", cons=C.STOP_GAINED, pos=200), _snv()) is None


# ------------------------------------------------------------------ PM1 (G10)
def test_hnf1a_promoter():
    assert _hnf1a_promoter_site(_pc("HNF1A", cons=C.UPSTREAM, hgvs_c="NM_000545.8:c.-218T>C"))
    assert _hnf1a_promoter_site(_pc("HNF1A", cons=C.UPSTREAM, hgvs_c="c.-197C>T")) is None
