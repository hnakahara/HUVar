"""PVS1 additions (Genome Medicine revision):

* generic SVI tree, NMD-escape branch: role of region unknown -> >10% of the
  protein altered = PVS1_Strong, <=10% = PVS1_Moderate (was N/A);
* ENIGMA BRCA1/BRCA2 Table 4 PTC rules; ITGA2B transmembrane critical region;
  GCK stop-loss; read-through frameshifts (MECP2/TCF4/CDKL5/UBE3A, SLC9A6).
"""
import pytest

from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import ConsequenceType as C, CriterionStrength as S
from acmg_classifier.pvs1.decision_tree import _nmd_branch
from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1


def _pc(gene, tx, cons, pos, hgvs_p, exon=None):
    return ConsequenceInfo(
        transcript_id=tx, gene_id="ENSG", gene_symbol=gene, consequence=cons,
        biotype="protein_coding", is_mane_select=True, protein_position=pos,
        hgvs_p=f"NP_0:{hgvs_p}", exon=exon,
    )


def _nmd(pc):
    return _nmd_branch(AnnotationData(gnomad=GnomADData(), consequences=[pc]), True, False, pc)


class TestSviTenPercent:
    def test_serpinc1_escape_over_10pct_strong(self):
        s, ev = _nmd(_pc("SERPINC1", "NM_000488.4", C.STOP_GAINED, 391, "p.Arg391Ter"))
        assert s == S.STRONG and ">10%" in ev

    def test_ldlr_escape_under_10pct_moderate(self):
        s, ev = _nmd(_pc("LDLR", "NM_000527.5", C.STOP_GAINED, 849, "p.Ser849Ter"))
        assert s == S.MODERATE and "<=10%" in ev

    def test_frameshift_counts_from_first_altered_residue(self):
        # SLC6A8 p.Phe564HisfsTer21: 72/635 aa altered (11.3%) -> Strong
        s, _ = _nmd(_pc("SLC6A8", "NM_005629.4", C.FRAMESHIFT, 564, "p.Phe564HisfsTer21"))
        assert s == S.STRONG

    def test_nmd_predicted_unchanged(self):
        s, _ = _nmd(_pc("LDLR", "NM_000527.5", C.STOP_GAINED, 100, "p.Gln100Ter"))
        assert s == S.VERY_STRONG


class TestBrcaTable4:
    @pytest.mark.parametrize("gene,tx,ptc,exp", [
        ("BRCA2", "NM_000059.4", 3247, S.VERY_STRONG),   # last exon, < p.T3310
        ("BRCA2", "NM_000059.4", 3309, S.VERY_STRONG),
        ("BRCA2", "NM_000059.4", 3326, S.NOT_MET),       # > p.E3309
        ("BRCA1", "NM_007294.4", 190, S.NOT_MET),        # exon 8(9)
        ("BRCA1", "NM_007294.4", 220, S.NOT_MET),        # exon 9(10)
        ("BRCA1", "NM_007294.4", 300, S.VERY_STRONG),
        ("BRCA1", "NM_007294.4", 1855, S.NOT_MET),
    ])
    def test_nonsense(self, gene, tx, ptc, exp):
        assert evaluate_vcep_pvs1(_pc(gene, tx, C.STOP_GAINED, ptc, f"p.Gln{ptc}Ter"))[0] == exp

    def test_frameshift_uses_ptc(self):
        # BRCA2 p.Cys3233TrpfsTer15 -> PTC 3247 -> PVS1
        pc = _pc("BRCA2", "NM_000059.4", C.FRAMESHIFT, 3233, "p.Cys3233TrpfsTer15")
        assert evaluate_vcep_pvs1(pc)[0] == S.VERY_STRONG


class TestItga2b:
    def test_alters_tm_domain_strong(self):
        pc = _pc("ITGA2B", "NM_000419.5", C.FRAMESHIFT, 992, "p.Glu992GlyfsTer42")
        s, ev = evaluate_vcep_pvs1(pc)
        assert s == S.STRONG and "critical region" in ev

    def test_after_tm_domain_10pct_moderate(self):
        pc = _pc("ITGA2B", "NM_000419.5", C.FRAMESHIFT, 1031, "p.Leu1031ArgfsTer10")
        assert evaluate_vcep_pvs1(pc)[0] == S.MODERATE


class TestStopLossAndReadthrough:
    def test_gck_stop_lost(self):
        pc = _pc("GCK", "NM_000162.5", C.STOP_LOST, 466, "p.Ter466ArgextTer144")
        assert evaluate_vcep_pvs1(pc)[0] == S.VERY_STRONG

    def test_other_gene_stop_lost_defers(self):
        pc = _pc("TCF4", "NM_001083962.2", C.STOP_LOST, 672, "p.Ter672GlnextTer20")
        assert evaluate_vcep_pvs1(pc) is None

    @pytest.mark.parametrize("gene,tx,pos,hp,exp", [
        ("MECP2", "NM_001110792.2", 490, "p.Pro490ArgfsTer40", S.VERY_STRONG),  # ptc 529 > 499
        ("MECP2", "NM_001110792.2", 490, "p.Pro490ArgfsTer5", S.MODERATE),      # ptc 494: band
        ("TCF4", "NM_001083962.2", 660, "p.Ser660AlafsTer?", S.VERY_STRONG),
        ("SLC9A6", "NM_001379110.1", 670, "p.Leu670SerfsTer30", S.MODERATE),
    ])
    def test_readthrough(self, gene, tx, pos, hp, exp):
        assert evaluate_vcep_pvs1(_pc(gene, tx, C.FRAMESHIFT, pos, hp))[0] == exp
