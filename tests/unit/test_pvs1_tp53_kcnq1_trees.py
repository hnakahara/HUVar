"""TP53 (PVS1 flowchart, Figure 1) and KCNQ1 (PVS1 decision tree v2025) PVS1."""
import pytest

from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType as CT, CriterionStrength as S
from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1


def _pc(gene, tx, csq, hgvs_c="", hgvs_p="", pos=None, intron=None):
    return ConsequenceInfo(transcript_id=tx, gene_id="g", gene_symbol=gene, consequence=csq,
                           biotype="protein_coding", is_mane_select=True, hgvs_c=hgvs_c,
                           hgvs_p=hgvs_p, protein_position=pos, intron=intron)


def _tp53(csq, **kw):
    return evaluate_vcep_pvs1(_pc("TP53", "NM_000546.6", csq, **kw))


@pytest.mark.parametrize("hgvs_p,pos,want", [
    ("NP_000537.3:p.Arg342Ter", 342, S.VERY_STRONG),
    ("NP_000537.3:p.Lys351Ter", 351, S.STRONG),
    ("NP_000537.3:p.Glu358Ter", 358, S.MODERATE),
])
def test_tp53_nonsense(hgvs_p, pos, want):
    assert _tp53(CT.STOP_GAINED, hgvs_p=hgvs_p, pos=pos)[0] == want


def test_tp53_frameshift_by_ptc_codon():
    # First altered codon 345 but the PTC lands at 360 (downstream of p.Leu350,
    # no NMD) -> PVS1_Moderate by the flowchart's "frameshift induced PTC" rule.
    r = _tp53(CT.FRAMESHIFT, hgvs_p="NP_000537.3:p.Ser345ArgfsTer16", pos=345)
    assert r[0] == S.MODERATE and "PTC at codon 360" in r[1]
    assert _tp53(CT.FRAMESHIFT, hgvs_p="NP_000537.3:p.Pro300GlnfsTer20",
                 pos=300)[0] == S.VERY_STRONG
    # Read-through of the natural stop codon.
    assert _tp53(CT.FRAMESHIFT, hgvs_p="NP_000537.3:p.Asp391GlufsTer26",
                 pos=391)[0] == S.MODERATE


@pytest.mark.parametrize("csq,hgvs_c,want", [
    (CT.SPLICE_DONOR, "NM_000546.6:c.-29+1G>A", S.NOT_MET),       # E1 donor
    (CT.SPLICE_ACCEPTOR, "NM_000546.6:c.-28-2A>G", S.NOT_MET),    # E2 acceptor
    (CT.SPLICE_DONOR, "NM_000546.6:c.74+1G>A", S.VERY_STRONG),    # E2 donor
    (CT.SPLICE_ACCEPTOR, "NM_000546.6:c.783-1G>A", S.NOT_MET),    # E8 acceptor
    (CT.SPLICE_DONOR, "NM_000546.6:c.1100+1G>A", S.MODERATE),     # E10 donor
    (CT.SPLICE_ACCEPTOR, "NM_000546.6:c.1101-2A>G", S.MODERATE),  # E11 acceptor
    (CT.SPLICE_ACCEPTOR, "NM_000546.6:c.376-1G>A", S.VERY_STRONG),
])
def test_tp53_splice(csq, hgvs_c, want):
    assert _tp53(csq, hgvs_c=hgvs_c)[0] == want


@pytest.mark.parametrize("csq,hgvs_c,want", [
    (CT.SPLICE_DONOR, "NM_000218.3:c.386+1G>A", S.VERY_STRONG),   # exon 1, NMD
    (CT.SPLICE_ACCEPTOR, "NM_000218.3:c.781-2A>G", S.STRONG),     # exon 6 in-frame
    (CT.SPLICE_DONOR, "NM_000218.3:c.1251+1G>A", S.STRONG),       # exon 9 in-frame
    (CT.SPLICE_DONOR, "NM_000218.3:c.1393+1G>A", S.VERY_STRONG),  # exon 10, NMD
    (CT.SPLICE_ACCEPTOR, "NM_000218.3:c.1686-1G>A", S.MODERATE),  # exon 14, no NMD
    (CT.SPLICE_ACCEPTOR, "NM_000218.3:c.1795-2A>G", S.MODERATE),  # exon 16 (last)
])
def test_kcnq1_splice(csq, hgvs_c, want):
    assert evaluate_vcep_pvs1(_pc("KCNQ1", "NM_000218.3", csq, hgvs_c=hgvs_c))[0] == want


@pytest.mark.parametrize("pos,want", [(500, S.VERY_STRONG), (600, S.MODERATE),
                                      (650, S.SUPPORTING)])
def test_kcnq1_truncation(pos, want):
    pc = _pc("KCNQ1", "NM_000218.3", CT.STOP_GAINED, hgvs_p=f"NP_000209.2:p.Arg{pos}Ter",
             pos=pos)
    assert evaluate_vcep_pvs1(pc)[0] == want


# --- GATM / IDUA / RPGR decision trees ------------------------------------------

def _ev(gene, tx, csq, **kw):
    return evaluate_vcep_pvs1(_pc(gene, tx, csq, **kw))


def test_gatm_cys407_and_splice():
    # NMD escape (exon 9) before Cys407 -> Strong even when <10% is removed.
    pc = _pc("GATM", "NM_001482.3", CT.STOP_GAINED, hgvs_p="NP_001473.1:p.Arg400Ter",
             pos=400)
    pc.exon = "9/9"
    r = evaluate_vcep_pvs1(pc)
    assert r[0] == S.STRONG and "critical region" in r[1]
    pc = _pc("GATM", "NM_001482.3", CT.STOP_GAINED, hgvs_p="NP_001473.1:p.Gln410Ter",
             pos=410)
    pc.exon = "9/9"
    assert evaluate_vcep_pvs1(pc)[0] == S.MODERATE
    for hc, want in (("c.69+1G>A", S.MODERATE), ("c.676-2A>G", S.STRONG),
                     ("c.978+1G>T", S.STRONG), ("c.1159+2T>C", S.MODERATE),
                     ("c.1160-1G>A", S.STRONG), ("c.288+1G>A", S.VERY_STRONG)):
        csq = CT.SPLICE_DONOR if "+" in hc else CT.SPLICE_ACCEPTOR
        assert _ev("GATM", "NM_001482.3", csq, hgvs_c=f"NM_001482.3:{hc}")[0] == want, hc


def test_idua_ptc_frameshift_and_start():
    tx = "NM_000203.5"
    assert _ev("IDUA", tx, CT.STOP_GAINED, hgvs_p="NP_000194.2:p.Trp402Ter",
               pos=402)[0] == S.VERY_STRONG
    assert _ev("IDUA", tx, CT.STOP_GAINED, hgvs_p="NP_000194.2:p.Gln600Ter",
               pos=600)[0] == S.MODERATE
    # Frameshift from codon 520 with its stop at 610 (escapes NMD): 91 codons altered.
    r = _ev("IDUA", tx, CT.FRAMESHIFT, hgvs_p="NP_000194.2:p.Leu520ArgfsTer91", pos=520)
    assert r[0] == S.STRONG and ">10%" in r[1]
    assert _ev("IDUA", tx, CT.FRAMESHIFT, hgvs_p="NP_000194.2:p.Leu620ArgfsTer5",
               pos=620)[0] == S.MODERATE
    assert _ev("IDUA", tx, CT.START_LOST)[0] == S.VERY_STRONG
    assert _ev("IDUA", tx, CT.SPLICE_ACCEPTOR,
               hgvs_c=f"{tx}:c.1829-2A>G")[0] == S.MODERATE


def test_rpgr_orf15_tail():
    tx = "NM_001034853.2"
    assert _ev("RPGR", tx, CT.STOP_GAINED, hgvs_p="NP_001030025.1:p.Glu900Ter",
               pos=900)[0] == S.VERY_STRONG
    assert _ev("RPGR", tx, CT.STOP_GAINED, hgvs_p="NP_001030025.1:p.Glu1140Ter",
               pos=1140)[0] == S.MODERATE


# --- InSiGHT MMR decision tree --------------------------------------------------

_MMR = {"MLH1": "NM_000249.4", "MSH2": "NM_000251.3", "MSH6": "NM_000179.3",
        "PMS2": "NM_000535.7"}


@pytest.mark.parametrize("gene,hgvs_p,pos,want", [
    ("MSH2", "p.Lys870ArgfsTer30", 870, S.MODERATE),     # PTC 899 > 891
    ("MSH2", "p.Lys850ArgfsTer10", 850, S.VERY_STRONG),  # PTC 859
    ("MLH1", "p.Arg755Ter", 755, S.MODERATE),
    ("PMS2", "p.Glu800Ter", 800, S.MODERATE),
    ("MSH6", "p.Thr1350ArgfsTer?", 1350, S.MODERATE),    # read-through, first codon
    ("MLH1", "p.Leu740SerfsTer?", 740, S.VERY_STRONG),
])
def test_mmr_truncation_by_ptc(gene, hgvs_p, pos, want):
    csq = CT.FRAMESHIFT if "fs" in hgvs_p else CT.STOP_GAINED
    assert _ev(gene, _MMR[gene], csq, hgvs_p=f"NP_0:{hgvs_p}", pos=pos)[0] == want


@pytest.mark.parametrize("gene,hgvs_c,want", [
    ("MLH1", "c.116+1G>A", S.VERY_STRONG),   # exon 1: frameshift + NMD
    ("MLH1", "c.208-2A>G", S.STRONG),        # exon 3 in-frame
    ("MSH2", "c.942+1G>A", S.STRONG),        # exon 5 in-frame
    ("MSH2", "c.2459-1G>A", S.MODERATE),     # exon 15: frameshift, no NMD
    ("MSH6", "c.4001+2T>C", S.MODERATE),     # exon 9: frameshift, no NMD
    ("PMS2", "c.989-1G>T", S.STRONG),        # exon 10 in-frame
    ("PMS2", "c.705+1G>A", S.STRONG),        # exon 6 in-frame
    ("PMS2", "c.803+1G>A", S.VERY_STRONG),   # exon 7: frameshift + NMD
])
def test_mmr_splice(gene, hgvs_c, want):
    csq = CT.SPLICE_DONOR if "+" in hgvs_c else CT.SPLICE_ACCEPTOR
    assert _ev(gene, _MMR[gene], csq, hgvs_c=f"{_MMR[gene]}:{hgvs_c}")[0] == want


@pytest.mark.parametrize("gene,want", [("MLH1", S.VERY_STRONG), ("MSH2", S.NOT_MET),
                                       ("MSH6", S.STRONG), ("PMS2", S.STRONG)])
def test_mmr_initiation(gene, want):
    assert _ev(gene, _MMR[gene], CT.START_LOST)[0] == want
