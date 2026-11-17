"""Rules implemented from the files attached to the CSpec documents (VCEP PVS1
trees / splice tables, frequency-exclusion lists, PM1 tables, corrections
documents), audited in November 2026."""
from pathlib import Path

import pytest

from acmg_classifier.criteria.allele_frequency import DiseaseThresholds, variant_listed
from acmg_classifier.criteria.bs2_genes import BS2Applicability
from acmg_classifier.criteria.pm1_hotspots import PM1_CO_REQUIREMENTS, PM1Hotspots
from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType as CT, CriterionStrength as S
from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1
from acmg_classifier.pvs1.vcep_pvs1_exons import SpliceExonOverrides

_SHARED = Path(__file__).resolve().parents[2] / "resources" / "shared"
_OVERRIDES = SpliceExonOverrides(_SHARED / "vcep_pvs1_splice_exons.tsv")


def _pc(gene, cq, tx="NM_x", pos=None, hgvs_c=None, hgvs_p=None, intron=None, exon=None):
    return ConsequenceInfo(
        transcript_id=tx, gene_id="ENSG", gene_symbol=gene, consequence=cq,
        biotype="protein_coding", protein_position=pos, hgvs_c=hgvs_c, hgvs_p=hgvs_p,
        intron=intron, exon=exon,
    )


# --- per-allele splice tables (ATM v1.6, PALB2 v1.3) -------------------------

@pytest.mark.parametrize("cq, hgvs, expected", [
    (CT.SPLICE_ACCEPTOR, "c.332-1G>A", S.STRONG),        # list B (HEAT in-frame)
    (CT.SPLICE_DONOR, "c.5496+1G>T", S.STRONG),          # list B
    (CT.SPLICE_DONOR, "c.5496+2T>A", S.VERY_STRONG),     # list A
    (CT.SPLICE_ACCEPTOR, "c.2467-2A>C", S.SUPPORTING),   # list C
    (CT.SPLICE_ACCEPTOR, "c.2467-2A>G", S.VERY_STRONG),  # list A
    (CT.SPLICE_ACCEPTOR, "c.8269-1G>T", S.SUPPORTING),   # list F
    (CT.SPLICE_ACCEPTOR, "c.8269-1G>A", S.VERY_STRONG),  # list E
    (CT.SPLICE_DONOR, "c.-31+1G>A", S.NOT_MET),          # 5'UTR
    (CT.SPLICE_DONOR, "c.6347+2T>C", S.NOT_MET),         # creates a functional GC donor
    (CT.SPLICE_DONOR, "c.6347+2T>A", S.VERY_STRONG),     # list D
    (CT.SPLICE_DONOR, "c.7515+2C>T", S.NOT_MET),         # GC donor improved
])
def test_atm_splice_table(cq, hgvs, expected):
    pc = _pc("ATM", cq, tx="NM_000051.4", hgvs_c=f"NM_000051.4:{hgvs}")
    s, ev = evaluate_vcep_pvs1(pc)
    assert s == expected
    if expected != S.VERY_STRONG:
        assert "ATM PVS1 guide v1.6" in ev


@pytest.mark.parametrize("cq, hgvs, expected", [
    (CT.SPLICE_ACCEPTOR, "c.49-1G>A", S.MODERATE),
    (CT.SPLICE_ACCEPTOR, "c.49-2A>T", S.MODERATE),
    (CT.SPLICE_ACCEPTOR, "c.49-2A>G", S.VERY_STRONG),   # list A
    (CT.SPLICE_DONOR, "c.211+2T>G", S.SUPPORTING),
    (CT.SPLICE_DONOR, "c.108+2T>C", S.NOT_MET),
    (CT.SPLICE_DONOR, "c.108+2T>A", S.VERY_STRONG),     # list C
    (CT.SPLICE_DONOR, "c.3350+2C>T", S.NOT_MET),
    (CT.SPLICE_ACCEPTOR, "c.212-1G>C", S.STRONG),       # list D
    (CT.SPLICE_ACCEPTOR, "c.1685-1G>C", S.STRONG),      # list D
    (CT.SPLICE_ACCEPTOR, "c.1685-1G>A", S.VERY_STRONG), # list A
])
def test_palb2_splice_table(cq, hgvs, expected):
    pc = _pc("PALB2", cq, tx="NM_024675.4", hgvs_c=f"NM_024675.4:{hgvs}")
    assert evaluate_vcep_pvs1(pc)[0] == expected


def test_site_table_needs_the_mane_transcript():
    pc = _pc("ATM", CT.SPLICE_ACCEPTOR, tx="NM_001351834.2",
             hgvs_c="NM_001351834.2:c.332-1G>A")
    assert evaluate_vcep_pvs1(pc)[0] == S.VERY_STRONG


# --- initiation codon -------------------------------------------------------

@pytest.mark.parametrize("gene", ["ATM", "PALB2", "OTC", "ABCA4"])
def test_initiation_codon_pvs1(gene):
    assert evaluate_vcep_pvs1(_pc(gene, CT.START_LOST))[0] == S.VERY_STRONG


def test_abca4_truncation_and_splice_defer_to_generic_tree():
    assert evaluate_vcep_pvs1(_pc("ABCA4", CT.STOP_GAINED, pos=500)) is None
    assert evaluate_vcep_pvs1(_pc("ABCA4", CT.SPLICE_DONOR, intron="5/49")) is None


# --- LGMD genes: re-initiation window, NMD, in-frame splice -----------------

@pytest.mark.parametrize("gene", ["DYSF", "SGCA", "SGCB", "SGCG", "SGCD", "CAPN3", "ANO5"])
def test_lgmd_first_100bp_moderate(gene):
    s, ev = evaluate_vcep_pvs1(_pc(gene, CT.STOP_GAINED, pos=20, hgvs_p="p.Gln20Ter"))
    assert s == S.MODERATE and "re-initiation" in ev
    assert evaluate_vcep_pvs1(_pc(gene, CT.FRAMESHIFT, pos=34))[0] == S.MODERATE


def test_lgmd_nmd_after_codon_34():
    pc = _pc("CAPN3", CT.STOP_GAINED, tx="NM_000070.3", pos=200, hgvs_p="p.Arg200Ter")
    assert evaluate_vcep_pvs1(pc)[0] == S.VERY_STRONG


def test_lgmd_inframe_splice_10pct_rule():
    # CAPN3 intron-5 acceptor skips exon 6 (in-frame, <10% of the protein).
    pc = _pc("CAPN3", CT.SPLICE_ACCEPTOR, tx="NM_000070.3", intron="5/23",
             hgvs_c="NM_000070.3:c.946-1G>A")
    s, ev = evaluate_vcep_pvs1(pc)
    assert s == S.MODERATE and "in-frame skip of exon 6/24" in ev


# --- splice_inframe on existing VCEP genes ------------------------------------

def test_mybpc3_inframe_exon_moderate_and_frameshift_exon_pvs1():
    inframe = _pc("MYBPC3", CT.SPLICE_DONOR, tx="NM_000256.3", intron="2/34")
    assert evaluate_vcep_pvs1(inframe)[0] == S.MODERATE
    frameshift = _pc("MYBPC3", CT.SPLICE_DONOR, tx="NM_000256.3", intron="5/34")
    assert evaluate_vcep_pvs1(frameshift)[0] == S.VERY_STRONG


def test_scid_inframe_skip_strong_with_review_note():
    # ADA intron-2 acceptor skips the in-frame exon 3.
    pc = _pc("ADA", CT.SPLICE_ACCEPTOR, tx="NM_000022.4", intron="2/11")
    s, ev = evaluate_vcep_pvs1(pc)
    assert s == S.STRONG and "downstream" in ev


def test_scn2a_inframe_exon_moderate():
    pc = _pc("SCN2A", CT.SPLICE_DONOR, tx="NM_001040142.2", intron="4/26")
    assert evaluate_vcep_pvs1(pc)[0] == S.MODERATE


# --- exon rows in vcep_pvs1_splice_exons.tsv -----------------------------------

@pytest.mark.parametrize("gene, cq, intron, expected", [
    ("PTEN", CT.SPLICE_ACCEPTOR, "2/8", S.STRONG),       # exon 3
    ("PTEN", CT.SPLICE_DONOR, "8/8", S.STRONG),          # exon 8
    ("ACADVL", CT.SPLICE_ACCEPTOR, "12/19", S.STRONG),   # exon 13
    ("ACADVL", CT.SPLICE_DONOR, "3/19", S.MODERATE),     # exon 3
    ("F9", CT.SPLICE_DONOR, "1/7", S.STRONG),            # exon 1
    ("F8", CT.SPLICE_ACCEPTOR, "25/25", S.STRONG),       # exon 26
    ("NEB", CT.SPLICE_ACCEPTOR, "54/182", S.VERY_STRONG),  # critical exon 55
    ("MECP2", CT.SPLICE_ACCEPTOR, "2/2", S.STRONG),      # e1 exon 3 = e2 exon 4
    ("HNF1A", CT.SPLICE_ACCEPTOR, "9/9", S.STRONG),
    ("UBE3A", CT.SPLICE_DONOR, "9/12", S.STRONG),
    ("SLC9A6", CT.SPLICE_DONOR, "4/17", S.STRONG),
    ("SLC6A8", CT.SPLICE_ACCEPTOR, "1/12", S.STRONG),    # gene without a PVS1 spec
])
def test_exon_rows(gene, cq, intron, expected):
    s, ev = evaluate_vcep_pvs1(_pc(gene, cq, intron=intron), _OVERRIDES)
    assert s == expected and "exon-aware" in ev


def test_neb_noncritical_inframe_exon_moderate():
    pc = _pc("NEB", CT.SPLICE_DONOR, tx="NM_001164508.2", intron="54/182")
    assert evaluate_vcep_pvs1(pc, _OVERRIDES)[0] == S.MODERATE


def test_gene_without_spec_and_row_defers():
    assert evaluate_vcep_pvs1(_pc("SLC6A8", CT.SPLICE_ACCEPTOR, intron="2/12"),
                              _OVERRIDES) is None


# --- frequency exclusions / BS2 / PM1 ----------------------------------------

def _thresholds():
    return DiseaseThresholds(_SHARED / "disease_prevalence.tsv")


@pytest.mark.parametrize("gene, hgvs_c, ba1_too", [
    ("ABCA4", "NM_000350.3:c.5603A>T", False),
    ("ABCA4", "NM_000350.3:c.4253+43G>A", False),
    ("HBB", "NM_000518.5:c.20A>T", True),        # Hb S
    ("HBB", "NM_000518.5:c.126_129delCTTT", True),
    ("HBA2", "NM_000517.6:c.95+2_95+6del", True),
])
def test_bs1_ba1_exclusions(gene, hgvs_c, ba1_too):
    gt = _thresholds().get(gene)
    assert variant_listed(gt.bs1_exclude, hgvs_c, None)
    assert bool(variant_listed(gt.ba1_exclude, hgvs_c, None)) == ba1_too


def test_dclre1c_bs2_homozygotes_only():
    bs2 = BS2Applicability(_SHARED / "disease_prevalence.tsv")
    assert bs2.hom_only("DCLRE1C")
    assert bs2.count("DCLRE1C") == 1


def test_acadvl_pm1_460_466():
    pm1 = PM1Hotspots(_SHARED / "pm1_hotspots.tsv")
    assert pm1.lookup("ACADVL", 463) == S.MODERATE
    assert pm1.lookup("ACADVL", 470) is None


@pytest.mark.parametrize("gene", ["FOXN1", "JAK3"])
def test_pm1_caveat_not_bs1_bs2_ba1(gene):
    assert PM1_CO_REQUIREMENTS[gene] == ((), ("BA1", "BS1", "BS2"))
