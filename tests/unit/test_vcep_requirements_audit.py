"""Rules corrected after the VCEP requirements audit of November 2026 (each gene's
CSpec criteria text against the tables / evaluators)."""
import csv
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria import pm5_conditions as cond
from acmg_classifier.criteria.bp_genes import BPApplicability
from acmg_classifier.criteria.comparator_tiers import PM5_RULES, tier_strength
from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
from acmg_classifier.criteria.pathogenic.pp3 import _splice_cutoff
from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator
from acmg_classifier.criteria.pm4_regions import PM4Regions
from acmg_classifier.criteria.ps1_paralog import PS1ParalogMap
from acmg_classifier.criteria.registry import _apply_pp3_pm1_cap
from acmg_classifier.criteria.revel_genes import RevelSpec
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion as C, Assembly, ConsequenceType as CT
from acmg_classifier.models.enums import CriterionStrength as S
from acmg_classifier.models.variant import VariantRecord

_SHARED = Path(__file__).resolve().parents[2] / "resources" / "shared"
_DP = _SHARED / "disease_prevalence.tsv"
_V = VariantRecord(chrom="chr1", pos=100, ref="C", alt="A", assembly=Assembly.GRCH38)


def _ann(gene, pos=10, aa="G/V", gnomad=None, cq=CT.MISSENSE):
    p3 = {"G": "Gly", "V": "Val", "Y": "Tyr", "C": "Cys", "M": "Met", "R": "Arg",
          "D": "Asp", "N": "Asn", "H": "His"}
    ref, alt = aa.split("/")
    return AnnotationData(gnomad=gnomad, consequences=[ConsequenceInfo(
        transcript_id="NM_1", gene_id="E", gene_symbol=gene, consequence=cq,
        biotype="protein_coding", is_mane_select=True, protein_position=pos,
        hgvs_p=f"NP_0:p.{p3[ref]}{pos}{p3[alt]}", amino_acids=aa)])


# --- PS1: MAP2K1/MAP2K2 and SOS1/SOS2 analogous residues -----------------------

def test_paralog_alignment_offsets():
    m = PS1ParalogMap(_SHARED / "ps1_paralog_map.tsv")
    assert m.analogs("MAP2K1", 130) == {"MAP2K2": 134}     # p.Tyr130 <-> p.Tyr134
    assert m.analogs("MAP2K2", 134) == {"MAP2K1": 130}
    assert m.analogs("SOS1", 269) == {"SOS2": 267}         # p.Met269 <-> p.Met267
    assert m.analogs("KCNQ1", 144) == {"KCNQ2": 114}       # unchanged
    assert m.analogs("SCN1A", 1) == {"SCN2A": 1, "SCN3A": 1, "SCN8A": 1}


_SCHEMA = """CREATE TABLE variants (variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT,
 alt TEXT, gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
 codon_position INTEGER, clinical_significance TEXT, review_status TEXT, star_rating INTEGER);"""


def _ps1(tmp_path, rows):
    db = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(db)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, gene_symbol, hgvs_p, amino_acid_change, "
        "codon_position, clinical_significance, review_status, star_rating) "
        "VALUES (?,?,?,?,?,?,?,?)", rows)
    con.commit()
    con.close()
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.clinvar_sqlite = db
    cfg.ps1_paralog_map_tsv = _SHARED / "ps1_paralog_map.tsv"
    return PS1Evaluator(cfg)


def test_map2k1_uses_the_aligned_map2k2_residue(tmp_path):
    rows = [("1", "MAP2K2", "NP_2:p.Tyr134Cys", "Y134C", 134, "Pathogenic", "x", 2)]
    r = _ps1(tmp_path, rows).evaluate(_V, _ann("MAP2K1", 130, "Y/C"))
    assert r.triggered and r.strength == S.STRONG and "paralogue" in r.evidence


def test_map2k1_same_number_in_map2k2_is_not_analogous(tmp_path):
    rows = [("1", "MAP2K2", "NP_2:p.Tyr130Cys", "Y130C", 130, "Pathogenic", "x", 2)]
    assert not _ps1(tmp_path, rows).evaluate(_V, _ann("MAP2K1", 130, "Y/C")).triggered


def test_sos2_uses_the_aligned_sos1_residue_lp_is_moderate(tmp_path):
    rows = [("1", "SOS1", "NP_1:p.Met269Arg", "M269R", 269, "Likely pathogenic", "x", 2)]
    r = _ps1(tmp_path, rows).evaluate(_V, _ann("SOS2", 267, "M/R"))
    assert r.triggered and r.strength == S.MODERATE


# --- PM4 ------------------------------------------------------------------------

def test_pm4_rows_mecp2_ube3a_in_mane_numbering():
    r = PM4Regions(_SHARED / "pm4_regions.tsv")
    assert r.in_tier("MECP2", 102) and r.in_tier("MECP2", 174) and not r.in_tier("MECP2", 95)
    assert r.in_tier("MECP2", 316)
    assert r.indel_strength("MECP2", 400) == "not_met"        # Pro-rich 393-417
    assert r.indel_strength("MECP2", 385) == S.MODERATE
    assert r.in_tier("UBE3A", 840)
    for g in ("MECP2", "CDKL5", "FOXG1", "TCF4", "UBE3A"):
        assert r.small_outside(g) == 2


def test_mecp2_ube3a_stoploss_strong():
    from acmg_classifier.criteria.pathogenic.pm4 import PM4Evaluator
    r = PM4Regions(_SHARED / "pm4_regions.tsv")
    assert r.stoploss_strength("MECP2") == S.STRONG
    assert r.stoploss_strength("UBE3A") == S.STRONG
    for g in ("CDKL5", "FOXG1", "TCF4"):
        assert r.stoploss_strength(g) is None          # default Moderate
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.pm4_regions_tsv = _SHARED / "pm4_regions.tsv"
    res = PM4Evaluator(cfg).evaluate(_V, _ann("MECP2", 499, cq=CT.STOP_LOST))
    assert res.triggered and res.strength == S.STRONG


def test_pm4_requires_and_conservation_rows():
    r = PM4Regions(_SHARED / "pm4_regions.tsv")
    assert r.requires("LDLR") == ("PM2",) and r.requires("RS1") == ("PM2",)
    assert r.requires("MECP2") == ()
    for g in ("AIPL1", "GUCY2D", "RPE65"):
        assert r.conserved_phylop(g) == 2.0


def _registry_stub():
    from acmg_classifier.criteria.registry import CriteriaRegistry
    reg = CriteriaRegistry.__new__(CriteriaRegistry)
    reg._pm4_regions = PM4Regions(_SHARED / "pm4_regions.tsv")
    return reg


@pytest.mark.parametrize("pm2_met, suppressed", [(True, False), (False, True)])
def test_ldlr_pm4_needs_pm2(pm2_met, suppressed):
    pm4 = CriteriaResult.met(C.PM4, S.MODERATE, "in-frame indel")
    pm2 = (CriteriaResult.met(C.PM2, S.SUPPORTING, "rare") if pm2_met
           else CriteriaResult.not_met(C.PM2, "common"))
    _registry_stub()._apply_pm4_exclusions(
        [pm4, pm2], _ann("LDLR", cq=CT.INFRAME_DELETION))
    assert pm4.suppressed is suppressed


# --- PM2: ITGA2B / ITGB3 ----------------------------------------------------------

def _pm2(tmp_path):
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.gene_inheritance_tsv = tmp_path / "missing.tsv"
    return PM2Evaluator(cfg)


@pytest.mark.parametrize("gene", ["ITGA2B", "ITGB3"])
def test_platelet_pm2_moderate_only_when_absent(tmp_path, gene):
    ev = _pm2(tmp_path)
    absent = ev.evaluate(_V, _ann(gene, gnomad=GnomADData(filter_pass=True, ac=0, af=0.0,
                                                          popmax_af=0.0)))
    assert absent.triggered and absent.strength == S.MODERATE
    rare = ev.evaluate(_V, _ann(gene, gnomad=GnomADData(filter_pass=True, ac=3, af=2e-6,
                                                        popmax_af=2e-5)))
    assert rare.triggered and rare.strength == S.SUPPORTING
    common = ev.evaluate(_V, _ann(gene, gnomad=GnomADData(filter_pass=True, ac=300, af=2e-4,
                                                          popmax_af=5e-4)))
    assert not common.triggered


# --- PP3 / BP4 --------------------------------------------------------------------

@pytest.mark.parametrize("gene", ["SCN1A", "SCN1B", "SCN2A", "SCN3A", "SCN8A"])
def test_scn_revel_pp3_capped_at_moderate(gene):
    assert RevelSpec(_DP).get(gene).pp3 == {S.SUPPORTING: 0.644, S.MODERATE: 0.773}


def _pair(pp3, pm1):
    return [CriteriaResult.met(C.PP3, pp3, "REVEL"), CriteriaResult.met(C.PM1, pm1, "hotspot")]


def test_pah_pp3_strong_withholds_pm1():
    res = _pair(S.STRONG, S.MODERATE)
    _apply_pp3_pm1_cap(res, _ann("PAH"))
    assert res[0].strength == S.STRONG and res[1].suppressed


def test_pah_pp3_moderate_keeps_pm1_moderate():
    res = _pair(S.MODERATE, S.MODERATE)
    _apply_pp3_pm1_cap(res, _ann("PAH"))
    assert not res[1].suppressed and res[1].strength == S.MODERATE


def test_pah_three_point_pp3_leaves_one_point_for_pm1():
    res = _pair(S.THREE_POINT, S.MODERATE)
    _apply_pp3_pm1_cap(res, _ann("PAH"))
    assert not res[1].suppressed and res[1].strength == S.SUPPORTING


def test_scn_pp3_capped_whatever_the_predictor():
    res = _pair(S.STRONG, S.MODERATE)
    _apply_pp3_pm1_cap(res, _ann("SCN2A"))
    assert res[0].strength == S.MODERATE and "caps PP3" in res[0].evidence
    assert not res[1].suppressed and res[1].strength == S.MODERATE


def test_other_genes_untouched():
    res = _pair(S.STRONG, S.MODERATE)
    _apply_pp3_pm1_cap(res, _ann("BRCA1"))
    assert res[0].strength == S.STRONG and not res[1].suppressed


def test_vhl_bp4_splice_only_and_pp3_splice_cutoff():
    assert BPApplicability(_DP).bp4_scope("VHL") == "splice_only"
    assert _splice_cutoff("VHL") == (0.50, True)


def test_gp1ba_gp1bb_bp4_need_zero_spliceai():
    rows = {r["gene_symbol"]: r for r in csv.DictReader(_DP.open(encoding="utf-8"),
                                                       delimiter="\t")}
    for g in ("GP1BA", "GP1BB", "GP9"):
        assert rows[g]["bp4_splice_cutoff"] == "0"
        assert rows[g]["bp4_missense_splice_max"] == "0"


# --- PM5 ----------------------------------------------------------------------------

@pytest.mark.parametrize("gene", ["DYSF", "CAPN3", "ANO5", "SGCA", "SGCB", "SGCG", "SGCD"])
def test_lgmd_pm5_conditions(gene):
    assert cond.SPLICE[gene] == (0.50, "both", "max")
    assert cond.REVEL_MIN[gene] == 0.7
    assert cond.splice_ok(gene, 0.6, "candidate") is False
    assert cond.splice_ok(gene, 0.1, "comparator", 0.1) is True


@pytest.mark.parametrize("n_p, n_lp, expected", [
    (2, 0, S.MODERATE), (1, 0, S.SUPPORTING), (0, 1, S.SUPPORTING),
    (0, 5, S.SUPPORTING), (1, 3, S.SUPPORTING), (0, 0, None),
])
def test_rpgr_pm5_moderate_needs_two_pathogenic(n_p, n_lp, expected):
    assert tier_strength(PM5_RULES["RPGR"], n_p, n_lp)[0] == expected


# =========================== second batch (conservative-side gaps) ===============

from acmg_classifier.criteria.allele_frequency import DiseaseThresholds  # noqa: E402
from acmg_classifier.criteria.benign.bp7 import BP7Evaluator  # noqa: E402
from acmg_classifier.criteria.benign.bs1 import BS1Evaluator  # noqa: E402
from acmg_classifier.criteria.benign.bs2 import BS2Evaluator  # noqa: E402
from acmg_classifier.criteria.bs2_genes import BS2Applicability  # noqa: E402
from acmg_classifier.criteria.pathogenic.pp3 import _splice_tier  # noqa: E402
from acmg_classifier.criteria.registry import _apply_bs1_bs2_exclusion  # noqa: E402


def _gd(**kw):
    return GnomADData(filter_pass=True, **kw)


def _rows():
    return {r["gene_symbol"]: r for r in csv.DictReader(_DP.open(encoding="utf-8"),
                                                       delimiter="\t")}


# --- PS1: RAS family and BRAF / RAF1 ------------------------------------------------

def test_ras_family_and_raf_alignment():
    m = PS1ParalogMap(_SHARED / "ps1_paralog_map.tsv")
    assert m.analogs("HRAS", 12) == {"KRAS": 12, "NRAS": 12, "MRAS": 22, "RIT1": 30,
                                     "RRAS2": 23}
    assert m.analogs("RIT1", 79)["HRAS"] == 61 and m.analogs("RRAS2", 72)["KRAS"] == 61
    assert m.analogs("BRAF", 594) == {"RAF1": 486}
    assert m.analogs("RAF1", 259) == {"BRAF": 365}


def test_mras_uses_the_aligned_kras_residue(tmp_path):
    rows = [("1", "KRAS", "NP_k:p.Gly12Val", "G12V", 12, "Pathogenic", "x", 2)]
    ev = _ps1(tmp_path, rows)
    r = ev.evaluate(_V, _ann("MRAS", 22, "G/V"))
    assert r.triggered and r.strength == S.STRONG
    assert not ev.evaluate(_V, _ann("MRAS", 12, "G/V")).triggered


def test_kras_still_finds_same_numbered_hras(tmp_path):
    rows = [("1", "HRAS", "NP_h:p.Gly12Val", "G12V", 12, "Pathogenic", "x", 2)]
    assert _ps1(tmp_path, rows).evaluate(_V, _ann("KRAS", 12, "G/V")).triggered


# --- BS1_Supporting -----------------------------------------------------------------

def _bs1(gene, **gd):
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.popmax_af_basis = True
    return BS1Evaluator(cfg).evaluate(_V, _ann(gene, gnomad=_gd(**gd)))


def test_bs1_supporting_thresholds_loaded():
    t = DiseaseThresholds(_DP)
    assert t.get("PTEN").bs1_supporting == 0.0000043
    assert t.get("BRCA1").bs1_supporting == 0.00002
    assert t.get("ABCA4").bs1_supporting == 0.00163
    assert t.get("GJB2").bs1_supporting == 0.0007 and t.get("MYO15A").bs1_supporting == 0.0007
    assert t.get("KCNQ4").bs1_supporting is None and t.get("COCH").bs1_supporting is None
    assert t.get("TP53").bs1_supporting is None


def test_pten_bs1_supporting_tier():
    r = _bs1("PTEN", faf95_popmax=1e-5, popmax_af=2e-5, af=1e-5, ac=20)
    assert r.triggered and r.strength == S.SUPPORTING
    assert _bs1("PTEN", faf95_popmax=5e-5, popmax_af=6e-5, af=5e-5, ac=90).strength == S.STRONG
    assert not _bs1("PTEN", faf95_popmax=1e-6, popmax_af=3e-6, af=1e-6, ac=2).triggered


def test_hht_bs1_supporting_with_two_homozygotes_is_strong():
    r = _bs1("ENG", faf95_popmax=0.001, popmax_af=0.0012, af=0.0005, ac=800, nhomalt=1)
    assert r.triggered and r.strength == S.SUPPORTING
    r = _bs1("ENG", faf95_popmax=0.001, popmax_af=0.0012, af=0.0005, ac=800, nhomalt=2)
    assert r.triggered and r.strength == S.STRONG and "homozygotes" in r.evidence


# --- PM2 clauses --------------------------------------------------------------------

@pytest.mark.parametrize("gene", ["DNM2", "NEB", "SCN1A", "SCN8A"])
def test_pm2_one_allele_allowed(tmp_path, gene):
    ev = _pm2(tmp_path)
    assert ev.evaluate(_V, _ann(gene, gnomad=_gd(ac=1, an=1_000_000, af=1e-6,
                                                 popmax_af=5e-4))).triggered


@pytest.mark.parametrize("gene", ["DNM2", "SCN1A", "SCN2A"])
def test_pm2_two_alleles_not_absent(tmp_path, gene):
    r = _pm2(tmp_path).evaluate(_V, _ann(gene, gnomad=_gd(ac=2, an=1_000_000, af=2e-6,
                                                           popmax_af=2e-6)))
    assert not r.triggered


def test_pm2_minimum_alleles_assessed(tmp_path):
    ev = _pm2(tmp_path)
    r = ev.evaluate(_V, _ann("SCN1A", gnomad=_gd(ac=0, an=5000, af=0.0, popmax_af=0.0)))
    assert not r.triggered and "alleles assessed" in r.evidence
    assert ev.evaluate(_V, _ann("SCN1A", gnomad=_gd(ac=0, an=20000, af=0.0,
                                                    popmax_af=0.0))).triggered
    assert not ev.evaluate(_V, _ann("NEB", gnomad=_gd(ac=0, an=1500, af=0.0,
                                                      popmax_af=0.0))).triggered


def test_mtm1_one_observation_in_females_only(tmp_path):
    ev = _pm2(tmp_path)
    assert ev.evaluate(_V, _ann("MTM1", gnomad=_gd(ac=1, nhemi=0, an=900_000, af=1e-6,
                                                   popmax_af=1e-5))).triggered
    assert not ev.evaluate(_V, _ann("MTM1", gnomad=_gd(ac=1, nhemi=1, an=900_000, af=1e-6,
                                                       popmax_af=1e-5))).triggered


@pytest.mark.parametrize("gene", ["ACVRL1", "ENG"])
def test_hht_pm2_fewer_than_six_alleles(tmp_path, gene):
    ev = _pm2(tmp_path)
    assert ev.evaluate(_V, _ann(gene, gnomad=_gd(ac=5, an=100_000, af=5e-5,
                                                 popmax_af=2e-4))).triggered
    assert not ev.evaluate(_V, _ann(gene, gnomad=_gd(ac=6, an=100_000, af=6e-5,
                                                     popmax_af=2e-4))).triggered


def test_apc_single_allele_cutoff(tmp_path):
    ev = _pm2(tmp_path)
    assert ev.evaluate(_V, _ann("APC", gnomad=_gd(ac=1, an=150_000, af=6.6e-6,
                                                  popmax_af=8e-6))).triggered
    assert not ev.evaluate(_V, _ann("APC", gnomad=_gd(ac=2, an=250_000, af=8e-6,
                                                      popmax_af=8e-6))).triggered


def test_tp53_global_af_and_multi_allele_group(tmp_path):
    ev = _pm2(tmp_path)
    # one-off allele in a small group: judged on the global AF
    assert ev.evaluate(_V, _ann("TP53", gnomad=_gd(ac=1, ac_grpmax=1, an=1_500_000, af=7e-7,
                                                   popmax_af=1e-4))).triggered
    # several alleles in one group above 0.00004 -> withheld
    r = ev.evaluate(_V, _ann("TP53", gnomad=_gd(ac=4, ac_grpmax=4, an=1_500_000, af=2.7e-6,
                                                popmax_af=6e-5)))
    assert not r.triggered and "GrpMax group" in r.evidence
    assert not ev.evaluate(_V, _ann("TP53", gnomad=_gd(ac=60, ac_grpmax=30, an=1_500_000,
                                                       af=4e-5, popmax_af=3e-5))).triggered


def test_cdh1_subpopulation_with_two_alleles(tmp_path):
    ev = _pm2(tmp_path)
    assert ev.evaluate(_V, _ann("CDH1", gnomad=_gd(ac=2, ac_grpmax=2, an=1_500_000, af=1.3e-6,
                                                   popmax_af=1.5e-5))).triggered
    assert not ev.evaluate(_V, _ann("CDH1", gnomad=_gd(ac=3, ac_grpmax=3, an=1_500_000,
                                                       af=2e-6, popmax_af=3e-5))).triggered


def test_dicer1_one_allele_per_subpopulation(tmp_path):
    ev = _pm2(tmp_path)
    assert ev.evaluate(_V, _ann("DICER1", gnomad=_gd(ac=2, ac_grpmax=1, an=1_500_000, af=1.3e-6,
                                                     popmax_af=2e-5))).triggered
    r = ev.evaluate(_V, _ann("DICER1", gnomad=_gd(ac=3, ac_grpmax=2, an=1_500_000, af=2e-6,
                                                  popmax_af=2e-5)))
    assert not r.triggered and "GrpMax group" in r.evidence


# --- BS2 from gnomAD counts ---------------------------------------------------------

def _bs2(gene, **gd):
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.bs2_min_homalt, cfg.bs2_min_hemi, cfg.bs2_min_het = 3, 3, 5
    return BS2Evaluator(cfg).evaluate(_V, _ann(gene, gnomad=_gd(**gd)))


def test_slc6a8_bs2_two_hemizygotes():
    assert _bs2("SLC6A8", nhemi=2, ac=3).triggered
    assert not _bs2("SLC6A8", nhemi=1, ac=3).triggered
    assert _bs2("SLC6A8", nhemi=0, nhomalt=2, ac=4).triggered


def test_pdha1_bs2_hemizygote_tiers():
    assert BS2Applicability(_DP).modes("PDHA1") == frozenset({"XL"})
    assert _bs2("PDHA1", nhemi=16, ac=40).strength == S.STRONG
    assert _bs2("PDHA1", nhemi=4, ac=10).strength == S.SUPPORTING
    assert not _bs2("PDHA1", nhemi=3, ac=10).triggered


def test_slc6a8_bs2_not_with_bs1():
    bs1 = CriteriaResult.met(C.BS1, S.STRONG, "gnomAD FAF95")
    bs2 = CriteriaResult.met(C.BS2, S.STRONG, "gnomAD X-linked: nhemi=5 >= 2")
    _apply_bs1_bs2_exclusion([bs1, bs2], _ann("SLC6A8"))
    assert bs2.suppressed and not bs1.suppressed
    curated = CriteriaResult.met(C.BS2, S.STRONG, "SLC6A8: ClinVar expert-panel applied BS2")
    _apply_bs1_bs2_exclusion([bs1, curated], _ann("SLC6A8"))
    assert not curated.suppressed


# --- PM4 / PM5 ----------------------------------------------------------------------

def test_pm4_pten_neb_rpgr_rows():
    r = PM4Regions(_SHARED / "pm4_regions.tsv")
    assert r.indel_strength("PTEN", 125) == S.MODERATE
    assert r.indel_strength("PTEN", 200) == "not_met"
    assert r.stoploss_strength("PTEN") == S.MODERATE
    assert r.deletion_strength("NEB", 2490) == S.STRONG
    assert r.deletion_strength("NEB", 2600) is None
    assert r.indel_strength("RPGR", 1100) == S.MODERATE and r.indel_strength("RPGR", 800) == "not_met"


def test_neb_exon55_deletion_strong_insertion_default():
    from acmg_classifier.criteria.pathogenic.pm4 import PM4Evaluator
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.pm4_regions_tsv = _SHARED / "pm4_regions.tsv"
    ev = PM4Evaluator(cfg)
    dele = VariantRecord(chrom="chr2", pos=100, ref="ACGT", alt="A", assembly=Assembly.GRCH38)
    ins = VariantRecord(chrom="chr2", pos=100, ref="A", alt="ACGT", assembly=Assembly.GRCH38)
    assert ev.evaluate(dele, _ann("NEB", 2490, cq=CT.INFRAME_DELETION)).strength == S.STRONG
    assert ev.evaluate(ins, _ann("NEB", 2490, cq=CT.INFRAME_INSERTION)).strength == S.MODERATE
    assert ev.evaluate(dele, _ann("NEB", 3000, cq=CT.INFRAME_DELETION)).strength == S.MODERATE


@pytest.mark.parametrize("n_p, n_lp, expected", [
    (2, 0, S.STRONG), (1, 0, S.MODERATE), (0, 2, S.MODERATE), (0, 1, S.SUPPORTING),
    (0, 6, S.MODERATE), (1, 3, S.MODERATE), (0, 0, None),
])
def test_cyp1b1_pm5_tiers(n_p, n_lp, expected):
    assert tier_strength(PM5_RULES["CYP1B1"], n_p, n_lp)[0] == expected


# --- PP3 / BP4 / BP7 ----------------------------------------------------------------

def test_revel_tiers_from_the_audit():
    rv = RevelSpec(_DP)
    assert rv.get("SCN1A").bp4 == {S.SUPPORTING: 0.29, S.MODERATE: 0.183}
    assert rv.get("PAH").bp4 == {S.SUPPORTING: 0.29, S.MODERATE: 0.183, S.STRONG: 0.016,
                                 S.VERY_STRONG: 0.003}
    for g in ("RPGR", "RS1"):
        assert rv.get(g).bp4[S.STRONG] == 0.016
    assert rv.get("AIPL1").pp3 == {S.SUPPORTING: 0.644, S.MODERATE: 0.774, S.STRONG: 0.932}


def test_abca4_splice_pp3_moderate_from_0_8():
    assert _splice_tier("ABCA4", 0.85, S.SUPPORTING) == S.MODERATE
    assert _splice_tier("ABCA4", 0.5, S.SUPPORTING) == S.SUPPORTING
    assert _splice_tier("BRCA1", 0.95, S.SUPPORTING) == S.SUPPORTING
    assert _splice_tier("ABCA4", 0.9, None) is None


def test_bp7_columns():
    rows = _rows()
    for g in ("ACVRL1", "APC", "ATM", "PALB2", "CDH1", "CTLA4", "CYP1B1", "DICER1", "GAMT",
              "IDUA", "LDLR", "MLH1", "PMS2", "OTC", "PAH", "PIK3R1", "PTEN", "SLC6A8", "VWF",
              "ABCA4", "CAPN3", "DYSF"):
        assert rows[g]["bp7_phylop"] == "na", g
    assert rows["BMPR2"]["bp7_phylop"] == "" and rows["KCNQ1"]["bp7_phylop"] == "2.0"
    for g in ("GP1BA", "GP1BB", "GP9", "GCK", "HNF1A", "HNF4A"):
        assert rows[g]["bp7_splice_cutoff"] == "0.2", g
    for g in ("GCK", "HNF1A", "HNF4A"):
        assert rows[g]["bp4_splice_cutoff"] == "0.2", g


def test_abca4_bp7_moderate():
    from acmg_classifier.models.annotation import SpliceScore
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.phylop_bigwig = Path("./missing.bw")
    ann = _ann("ABCA4", cq=CT.SYNONYMOUS)
    ann.splice = SpliceScore(tool="spliceai", max_delta=0.01, is_available=True)
    r = BP7Evaluator(cfg).evaluate(_V, ann)
    assert r.triggered and r.strength == S.MODERATE
    ann2 = _ann("BRCA2", cq=CT.SYNONYMOUS)
    ann2.splice = SpliceScore(tool="spliceai", max_delta=0.01, is_available=True)
    assert BP7Evaluator(cfg).evaluate(_V, ann2).strength == S.SUPPORTING


# =========================== third batch ==========================================

def test_rpgr_rs1_bp4_stand_alone_tier_is_very_strong():
    from acmg_classifier.criteria.benign.bp4 import _revel_bp4
    for g in ("RPGR", "RS1"):
        tiers = RevelSpec(_DP).get(g).bp4
        assert tiers[S.VERY_STRONG] == 0.003 and tiers[S.STRONG] == 0.016
        assert _revel_bp4(0.002, tiers) == S.VERY_STRONG
        assert _revel_bp4(0.01, tiers) == S.STRONG
        assert _revel_bp4(0.1, tiers) == S.MODERATE
        assert _revel_bp4(0.25, tiers) == S.SUPPORTING
        assert _revel_bp4(0.4, tiers) is None


def test_il2rg_pm2_faf_and_no_hemizygotes(tmp_path):
    ev = _pm2(tmp_path)
    ok = ev.evaluate(_V, _ann("IL2RG", gnomad=_gd(ac=3, nhemi=0, faf95_popmax=5e-5,
                                                  popmax_af=9e-5, af=2e-6)))
    assert ok.triggered and ok.strength == S.SUPPORTING and "FAF95" in ok.evidence
    hemi = ev.evaluate(_V, _ann("IL2RG", gnomad=_gd(ac=3, nhemi=1, faf95_popmax=5e-5,
                                                    popmax_af=9e-5, af=2e-6)))
    assert not hemi.triggered and "hemizygote" in hemi.evidence
    assert not ev.evaluate(_V, _ann("IL2RG", gnomad=_gd(ac=90, nhemi=0, faf95_popmax=2e-4,
                                                        popmax_af=3e-4, af=6e-5))).triggered


_RASOPATHY = ["BRAF", "HRAS", "KRAS", "LZTR1", "MAP2K1", "MAP2K2", "MRAS", "NRAS", "PPP1CB",
              "PTPN11", "RAF1", "RIT1", "RRAS2", "SHOC2", "SOS1", "SOS2"]


@pytest.mark.parametrize("gene", _RASOPATHY)
def test_rasopathy_pm5_and_pm1_combine_at_moderate(gene):
    from acmg_classifier.criteria.pm1_hotspots import PM1_CO_REQUIREMENTS
    from acmg_classifier.criteria.pm5_genes import PM5Spec
    from acmg_classifier.criteria.registry import _apply_pm1_co_requirements
    assert PM5Spec(_DP).excludes(gene) == ()
    assert PM1_CO_REQUIREMENTS[gene] == ((), ("PM5:Strong",))
    res = [CriteriaResult.met(C.PM1, S.MODERATE, "domain"),
           CriteriaResult.met(C.PM5, S.MODERATE, "comparator")]
    _apply_pm1_co_requirements(res, _ann(gene))
    assert not res[0].suppressed and not res[1].suppressed
    res = [CriteriaResult.met(C.PM1, S.MODERATE, "domain"),
           CriteriaResult.met(C.PM5, S.STRONG, "curated")]
    _apply_pm1_co_requirements(res, _ann(gene))
    assert res[0].suppressed and not res[1].suppressed


def test_other_vceps_keep_their_pm5_pm1_exclusion():
    from acmg_classifier.criteria.pm5_genes import PM5Spec
    spec = PM5Spec(_DP)
    for g in ("RUNX1", "MYH7", "ACADVL", "DICER1"):
        assert C.PM1 in spec.excludes(g), g


# =========================== PM5 at analogous residues (RASopathy) ================

def _pm5(tmp_path, rows):
    from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
    db = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(db)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, hgvs_p, "
        "amino_acid_change, codon_position, clinical_significance, review_status, "
        "star_rating) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit()
    con.close()
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.clinvar_sqlite = db
    cfg.ps1_paralog_map_tsv = _SHARED / "ps1_paralog_map.tsv"
    cfg.pm5_min_stars = 2
    cfg.pm5_ptc = False
    return PM5Evaluator(cfg)


def _row(vid, gene, p3, aac, pos, sig="Pathogenic"):
    return (vid, "chr9", 5000 + pos, "C", "T", gene, f"NP_x:p.{p3}", aac, pos, sig, "x", 2)


def test_pm5_raf1_from_analogous_braf_residue(tmp_path):
    # BRAF p.Asp594 <-> RAF1 p.Asp486: BRAF D594G (P) supports RAF1 D486N.
    ev = _pm5(tmp_path, [_row("1", "BRAF", "Asp594Gly", "D594G", 594)])
    r = ev.evaluate(_V, _ann("RAF1", 486, "D/N"))
    assert r.triggered and r.strength == S.MODERATE
    assert "analogous residue of BRAF" in r.evidence
    # the same substitution in the paralogue is PS1 material, not PM5
    assert not ev.evaluate(_V, _ann("RAF1", 486, "D/G")).triggered
    # the same residue NUMBER in BRAF is not the analogous residue
    assert not ev.evaluate(_V, _ann("RAF1", 594, "D/N")).triggered


def test_pm5_mras_from_kras_and_lp_counts(tmp_path):
    ev = _pm5(tmp_path, [_row("1", "KRAS", "Gly12Asp", "G12D", 12, "Likely pathogenic")])
    r = ev.evaluate(_V, _ann("MRAS", 22, "G/V"))
    assert r.triggered and r.strength == S.MODERATE


def test_pm5_two_analogous_changes_capped_at_moderate(tmp_path):
    ev = _pm5(tmp_path, [_row("1", "HRAS", "Gly12Asp", "G12D", 12),
                         _row("2", "NRAS", "Gly12Ser", "G12S", 12)])
    r = ev.evaluate(_V, _ann("KRAS", 12, "G/V"))
    assert r.triggered and r.strength == S.MODERATE and "probands" in r.evidence


def test_pm5_wild_type_residue_must_match(tmp_path):
    # MAP2K1 p.Tyr130 <-> MAP2K2 p.Tyr134; a comparator needs the same wild type.
    ev = _pm5(tmp_path, [_row("1", "MAP2K2", "Tyr134Cys", "Y134C", 134)])
    assert ev.evaluate(_V, _ann("MAP2K1", 130, "Y/H")).triggered
    assert not ev.evaluate(_V, _ann("MAP2K1", 131, "Y/H")).triggered


def test_pm5_single_gene_rasopathy_unaffected(tmp_path):
    ev = _pm5(tmp_path, [_row("1", "KRAS", "Gly12Asp", "G12D", 12)])
    assert not ev.evaluate(_V, _ann("PTPN11", 12, "G/V")).triggered
    assert "PTPN11" not in cond.PARALOG_COMPARATORS and "SHOC2" not in cond.PARALOG_COMPARATORS


# =========================== general "VCEP Specifications" text (CSpec page PDFs) ==

from acmg_classifier.criteria.allele_frequency import min_alleles_block  # noqa: E402
from acmg_classifier.criteria.registry import (  # noqa: E402
    _apply_bp7_bp4_exclusion, _apply_combination_caps,
)


def test_pm5_conditions_from_general_text():
    assert cond.SPLICE["TP53"] == (0.20, "candidate", "max")
    assert cond.SPLICE["CYP1B1"] == (0.20, "candidate", "le")
    assert cond.SPLICE["OTC"] == (0.20, "both", "le")
    assert cond.SPLICE["IDUA"] == (0.10, "both", "le")
    assert cond.SPLICE["AIPL1"] == (0.20, "both", "max")
    assert cond.splice_ok("RUNX1", 0.2, "candidate") is True
    assert cond.splice_ok("RUNX1", 0.21, "candidate") is False
    assert {"CYP1B1", "AIPL1"} <= cond.REQUIRES_PP3
    rows = _rows()
    assert rows["TP53"]["pm5_grantham"] == "ge" and rows["CYP1B1"]["pm5_grantham"] == "ge"


@pytest.mark.parametrize("gene, delta, blocked", [
    ("TP53", 0.25, True), ("TP53", 0.1, False), ("MYOC", 0.2, False), ("MYOC", 0.3, True),
    ("IDUA", 0.1, True), ("OTC", 0.2, False), ("ATM", 0.2, True), ("BRCA1", 0.9, False),
])
def test_ps1_candidate_splice_gate(tmp_path, gene, delta, blocked):
    from acmg_classifier.models.annotation import SpliceScore
    ev = _ps1(tmp_path, [])
    ann = _ann(gene, 100, "G/V")
    ann.splice = SpliceScore(tool="spliceai", max_delta=delta, is_available=True)
    r = ev.evaluate(_V, ann)
    assert ("predicted to affect splicing" in r.evidence) is blocked


@pytest.mark.parametrize("gene", ["SCN1A", "SCN2A", "SCN3A", "SCN8A", "ACVRL1", "ENG"])
def test_pm1_not_with_pm5_strong(gene):
    from acmg_classifier.criteria.registry import _apply_pm1_co_requirements
    res = [CriteriaResult.met(C.PM1, S.MODERATE, "x"), CriteriaResult.met(C.PM5, S.STRONG, "y")]
    _apply_pm1_co_requirements(res, _ann(gene))
    assert res[0].suppressed
    res = [CriteriaResult.met(C.PM1, S.MODERATE, "x"), CriteriaResult.met(C.PM5, S.MODERATE, "y")]
    _apply_pm1_co_requirements(res, _ann(gene))
    assert not res[0].suppressed


def test_rpe65_pm1_not_with_bs1():
    from acmg_classifier.criteria.registry import _apply_pm1_co_requirements
    res = [CriteriaResult.met(C.PM1, S.MODERATE, "x"), CriteriaResult.met(C.BS1, S.STRONG, "y")]
    _apply_pm1_co_requirements(res, _ann("RPE65"))
    assert res[0].suppressed


def _combo(gene, **strengths):
    res = [CriteriaResult.met(getattr(C, k), v, k) for k, v in strengths.items()]
    _apply_pp3_pm1_cap(res, _ann(gene))
    _apply_combination_caps(res, _ann(gene))
    return {r.criterion.value: (None if r.suppressed else r.strength) for r in res}


def test_combination_caps():
    assert _combo("OTC", PP3=S.STRONG, PM1=S.MODERATE) == {"PP3": S.STRONG, "PM1": None}
    assert _combo("F8", PM1=S.STRONG, PM5=S.MODERATE) == {"PM1": S.STRONG, "PM5": None}
    assert _combo("F9", PM1=S.MODERATE, PM5=S.MODERATE) == {"PM1": S.MODERATE, "PM5": S.MODERATE}
    assert _combo("TP53", PS1=S.STRONG, PM5=S.MODERATE) == {"PS1": S.STRONG, "PM5": None}
    assert _combo("TP53", PS1=S.MODERATE, PM5=S.MODERATE) == {"PS1": S.MODERATE,
                                                               "PM5": S.MODERATE}
    assert _combo("MYOC", PP3=S.STRONG, PM5=S.MODERATE) == {"PP3": S.STRONG,
                                                            "PM5": S.SUPPORTING}
    assert _combo("MYOC", PP3=S.STRONG, PS1=S.STRONG) == {"PP3": S.MODERATE, "PS1": S.STRONG}
    # CYP1B1: PP3 + PM1 <= 4, then PP3 + PM1 + PM5 <= 5
    assert _combo("CYP1B1", PP3=S.MODERATE, PM1=S.MODERATE, PM5=S.MODERATE) == {
        "PP3": S.MODERATE, "PM1": S.MODERATE, "PM5": S.SUPPORTING}
    assert _combo("CYP1B1", PP3=S.STRONG, PM1=S.MODERATE, PS1=S.STRONG) == {
        "PP3": S.MODERATE, "PM1": None, "PS1": S.STRONG}
    assert _combo("BRCA1", PP3=S.STRONG, PM5=S.STRONG) == {"PP3": S.STRONG, "PM5": S.STRONG}


def test_abca4_bp4_not_with_bp7():
    res = [CriteriaResult.met(C.BP7, S.MODERATE, "x"), CriteriaResult.met(C.BP4, S.SUPPORTING, "y")]
    _apply_bp7_bp4_exclusion(res, _ann("ABCA4", cq=CT.SYNONYMOUS))
    assert res[1].suppressed and not res[0].suppressed
    res = [CriteriaResult.met(C.BP7, S.SUPPORTING, "x"), CriteriaResult.met(C.BP4, S.SUPPORTING, "y")]
    _apply_bp7_bp4_exclusion(res, _ann("PAH", cq=CT.SYNONYMOUS))
    assert not res[1].suppressed


def test_min_alleles_for_point_estimate_ba1_bs1():
    assert "variant allele" in min_alleles_block("CDH1", _gd(ac_grpmax=3, an_grpmax=50000))
    assert "GrpMax population" in min_alleles_block("MYOC", _gd(ac_grpmax=6, an_grpmax=1500))
    assert min_alleles_block("CDH1", _gd(ac_grpmax=6, an_grpmax=50000)) == ""
    assert min_alleles_block("F8", _gd(ac_grpmax=4)) != "" and min_alleles_block("F8", _gd(ac_grpmax=5)) == ""
    assert min_alleles_block("BRCA1", _gd(ac_grpmax=1)) == ""
    r = _bs1("CDH1", faf95_popmax=0.002, popmax_af=0.08, af=0.001, ac=4, ac_grpmax=4,
             an_grpmax=50)
    assert not r.triggered and "CDH1 VCEP" in r.evidence


def test_tsv_values_from_general_text():
    rows = _rows()
    assert rows["CDH1"]["pm2_min_depth"] == "30"
    for g in ("HBB", "HBA2"):
        assert (rows[g]["pm2_min_depth"], rows[g]["bp7_splice_cutoff"], rows[g]["bp7_phylop"]) \
            == ("20", "0.3", "0")
    assert rows["PTEN"]["bp4_splice_cutoff"] == "0.2" and _splice_cutoff("PTEN") == (0.50, False)
    # APC non-cancer subset deliberately not used (companion DB covers chr13/chr17 only)
    assert rows["APC"]["pm2_subset"] == ""
    assert "p.Ala416Ser" in rows["SERPINC1"]["bs1_exclude"]


def test_aipl1_pm4_region_needs_five_residues():
    r = PM4Regions(_SHARED / "pm4_regions.tsv")
    assert r.denied_below("AIPL1", 350, 3) == 5 and r.denied_below("AIPL1", 350, 5) is None
    assert r.denied_below("AIPL1", 100, 1) is None


# =========================== RYR1 PM4_Strong, CDH1 PP3, BRCA1/2 read depth =========

def test_ryr1_pore_region_deletion_strong():
    from acmg_classifier.criteria.pathogenic.pm4 import PM4Evaluator
    r = PM4Regions(_SHARED / "pm4_regions.tsv")
    assert r.deletion_strength("RYR1", 4800) == S.STRONG
    assert r.deletion_strength("RYR1", 4950) == S.STRONG
    assert r.deletion_strength("RYR1", 4951) is None
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.pm4_regions_tsv = _SHARED / "pm4_regions.tsv"
    ev = PM4Evaluator(cfg)
    dele = VariantRecord(chrom="chr19", pos=100, ref="ACGT", alt="A", assembly=Assembly.GRCH38)
    ins = VariantRecord(chrom="chr19", pos=100, ref="A", alt="ACGT", assembly=Assembly.GRCH38)
    assert ev.evaluate(dele, _ann("RYR1", 4850, cq=CT.INFRAME_DELETION)).strength == S.STRONG
    assert ev.evaluate(ins, _ann("RYR1", 4850, cq=CT.INFRAME_INSERTION)).strength == S.MODERATE
    assert ev.evaluate(dele, _ann("RYR1", 2000, cq=CT.INFRAME_DELETION)).strength == S.MODERATE


def test_cdh1_pp3_not_at_last_nucleotide_of_exon3():
    from acmg_classifier.criteria.pathogenic.pp3 import _pp3_excluded
    assert _pp3_excluded("CDH1", "NM_004360.5:c.387G>A")
    assert _pp3_excluded("CDH1", "c.387G>T")
    assert not _pp3_excluded("CDH1", "NM_004360.5:c.387+1G>A")
    assert not _pp3_excluded("CDH1", "NM_004360.5:c.3870G>A")
    assert not _pp3_excluded("CDH1", "NM_004360.5:c.386A>G")
    assert not _pp3_excluded("BRCA1", "c.387G>A") and not _pp3_excluded("CDH1", None)


class _Cov:
    available = True

    def __init__(self, dp):
        self._dp = dp

    def mean_depth(self, chrom, start, end):
        return self._dp


def test_brca_ba1_bs1_need_read_depth_20():
    from acmg_classifier.criteria.allele_frequency import depth_block
    from acmg_classifier.criteria.benign.ba1 import BA1Evaluator
    assert "read depth 12 < 20" in depth_block(_Cov(12), "BRCA1", _V)
    assert depth_block(_Cov(20), "BRCA2", _V) == ""
    assert depth_block(_Cov(None), "BRCA1", _V) == ""       # no value -> gate skipped
    assert depth_block(None, "BRCA1", _V) == ""
    assert depth_block(_Cov(5), "TP53", _V) == ""
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.popmax_af_basis = True
    ann = _ann("BRCA1", gnomad=_gd(faf95_popmax=0.01, popmax_af=0.012, af=0.005, ac=900,
                                  faf95_non_cancer=0.01))
    for cls, crit in ((BA1Evaluator, "BA1"), (BS1Evaluator, "BS1")):
        ev = cls(cfg)
        ev._coverage = _Cov(30)
        assert ev.evaluate(_V, ann).triggered, crit
        ev._coverage = _Cov(10)
        r = ev.evaluate(_V, ann)
        assert not r.triggered and "read depth" in r.evidence, crit


# =========================== PS1 comparator splice (LCA), PM2 population size ======

class _Lookup:
    def __init__(self, scores):
        self._s = scores

    def predict(self, v):
        from acmg_classifier.models.annotation import SpliceScore
        d = self._s.get(v.pos)
        return SpliceScore(tool="spliceai", max_delta=d, is_available=d is not None)


@pytest.mark.parametrize("gene", ["GUCY2D", "AIPL1"])
def test_ps1_comparator_must_not_be_spliceogenic(tmp_path, gene):
    db = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(db)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, hgvs_p, "
        "amino_acid_change, codon_position, clinical_significance, review_status, "
        "star_rating) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        [("1", "chr1", 101, "G", "T", gene, "NP_x:p.Gly100Val", "G100V", 100, "Pathogenic",
          "reviewed by expert panel", 3)])
    con.commit()
    con.close()
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.clinvar_sqlite = db
    cfg.ps1_paralog_map_tsv = _SHARED / "ps1_paralog_map.tsv"
    cfg.assembly = Assembly.GRCH38
    ev = PS1Evaluator(cfg)
    ev._splice_lookup_ready = True
    ev._splice_lookup = _Lookup({101: 0.05})
    assert ev.evaluate(_V, _ann(gene, 100, "G/V")).triggered
    ev._splice_lookup = _Lookup({101: 0.4})
    r = ev.evaluate(_V, _ann(gene, 100, "G/V"))
    assert not r.triggered and "SpliceAI > 0.1" in r.evidence
    ev._splice_lookup = None                      # no SpliceAI files: comparator kept
    r = ev.evaluate(_V, _ann(gene, 100, "G/V"))
    assert r.triggered and "not assessed" in r.evidence


def test_pm2_small_grpmax_population_falls_back_to_global_af(tmp_path):
    ev = _pm2(tmp_path)
    # MYOC cutoff 0.0001: a 2-in-4,000 population (AF 5e-4) is below 10,000 alleles
    r = ev.evaluate(_V, _ann("MYOC", gnomad=_gd(ac=2, ac_grpmax=2, an_grpmax=4000, an=1_500_000,
                                                af=1.3e-6, popmax_af=5e-4)))
    assert r.triggered and "global AF" in r.evidence
    r = ev.evaluate(_V, _ann("MYOC", gnomad=_gd(ac=20, ac_grpmax=20, an_grpmax=40000,
                                                an=1_500_000, af=1.3e-5, popmax_af=5e-4)))
    assert not r.triggered
    # other genes keep the GrpMax frequency whatever the population size
    r = ev.evaluate(_V, _ann("HBB", gnomad=_gd(ac=2, ac_grpmax=2, an_grpmax=4000, an=1_500_000,
                                               af=1.3e-6, popmax_af=5e-4)))
    assert not r.triggered


# =========================== CDH1 (Luo et al. J Med Genet 2023 supplement) =========

def _cdh1_pc(cq, hgvs_c, tx="NM_004360.5"):
    return ConsequenceInfo(transcript_id=tx, gene_id="E", gene_symbol="CDH1", consequence=cq,
                           biotype="protein_coding", is_mane_select=True,
                           hgvs_c=f"{tx}:{hgvs_c}")


@pytest.mark.parametrize("hgvs, cq, expected", [
    ("c.387+1G>A", CT.SPLICE_DONOR, S.MODERATE),        # exon 3 donor (in-frame by RNA)
    ("c.387+2T>C", CT.SPLICE_DONOR, S.MODERATE),
    ("c.1566-1G>C", CT.SPLICE_ACCEPTOR, S.MODERATE),    # exon 11 acceptor
    ("c.1712-2A>C", CT.SPLICE_ACCEPTOR, S.MODERATE),    # exon 12 acceptor
    ("c.1565+1G>A", CT.SPLICE_DONOR, S.STRONG),         # default
    ("c.163+1G>A", CT.SPLICE_DONOR, S.STRONG),
])
def test_cdh1_pvs1_splicing_table(hgvs, cq, expected):
    from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1
    assert evaluate_vcep_pvs1(_cdh1_pc(cq, hgvs))[0] == expected


def test_cdh1_initiation_codon_pvs1():
    from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1
    assert evaluate_vcep_pvs1(_cdh1_pc(CT.START_LOST, "c.1A>G"))[0] == S.VERY_STRONG


@pytest.mark.parametrize("hgvs, cq, expected", [
    ("c.48+1G>A", CT.SPLICE_DONOR, S.SUPPORTING),
    ("c.49-2A>G", CT.SPLICE_ACCEPTOR, S.SUPPORTING),
    ("c.1137+2T>C", CT.SPLICE_DONOR, S.SUPPORTING),
    ("c.2440-2A>G", CT.SPLICE_ACCEPTOR, S.SUPPORTING),
    ("c.163+1G>A", CT.SPLICE_DONOR, S.NOT_MET),         # PVS1_Strong only
    ("c.1009-1G>A", CT.SPLICE_ACCEPTOR, S.NOT_MET),
    ("c.2295+2T>C", CT.SPLICE_DONOR, S.NOT_MET),
    ("c.387+1G>A", CT.SPLICE_DONOR, S.NOT_MET),         # PVS1_Moderate site
])
def test_cdh1_pm5_supporting_only_at_listed_splice_sites(hgvs, cq, expected):
    from acmg_classifier.criteria.pm5_ptc import evaluate_pm5_ptc
    assert evaluate_pm5_ptc(_cdh1_pc(cq, hgvs), AnnotationData())[0] == expected


def test_filtering_af_9999():
    from acmg_classifier.criteria.allele_frequency import filtering_af
    import math
    f = filtering_af(10, 100000)
    lam = f * 100000
    # P(X >= 10 | lam) == 1e-4 at the bound
    p = 1 - sum(math.exp(-lam) * lam ** k / math.factorial(k) for k in range(10))
    assert abs(p - 1e-4) < 1e-6
    assert filtering_af(10, 100000, 0.95) > f                 # FAF95 is higher
    assert filtering_af(0, 1000) == 0.0 and filtering_af(None, 1000) is None


@pytest.mark.parametrize("gene", ["CDH1", "F8", "F9"])
def test_ba1_bs1_on_99_99_ci(gene):
    from acmg_classifier.criteria.benign.ba1 import BA1Evaluator
    assert DiseaseThresholds(_DP).get(gene).af_basis == "faf9999"
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _DP
    cfg.popmax_af_basis = True
    # point AF 0.0025 (>= BA1 0.002 for CDH1) but the 99.99% bound of 10/4000 is lower
    ann = _ann(gene, gnomad=_gd(faf95_popmax=0.0013, popmax_af=0.0025, af=0.0001, ac=10,
                               ac_grpmax=10, an_grpmax=4000))
    r = BA1Evaluator(cfg).evaluate(_V, ann)
    assert "99.99% CI" in r.evidence
    if gene == "CDH1":
        assert not r.triggered
    ann = _ann(gene, gnomad=_gd(faf95_popmax=0.0013, popmax_af=0.0025, af=0.0001, ac=10))
    r = BS1Evaluator(cfg).evaluate(_V, ann)
    assert "not computable" in r.evidence


# =========================== CDH1 v3.1 PVS1 decision tree (c. ranges) ===============

@pytest.mark.parametrize("cq, hgvs, expected", [
    (CT.STOP_GAINED, "c.2386C>T", S.VERY_STRONG),       # nonsense c.4-2388
    (CT.STOP_GAINED, "c.2389G>T", S.STRONG),            # c.2389-2508
    (CT.STOP_GAINED, "c.2506G>T", S.STRONG),
    (CT.STOP_GAINED, "c.2512C>T", S.MODERATE),          # c.2509-2646
    (CT.STOP_GAINED, "c.2547_2548insTAA", S.MODERATE),
    (CT.FRAMESHIFT, "c.2370del", S.VERY_STRONG),        # -1: c.4-2375
    (CT.FRAMESHIFT, "c.2376del", S.STRONG),             # -1: c.2376-2443
    (CT.FRAMESHIFT, "c.2450del", S.MODERATE),           # -1: c.2444-2646
    (CT.FRAMESHIFT, "c.2352dup", S.VERY_STRONG),        # +1: c.4-2352
    (CT.FRAMESHIFT, "c.2490dup", S.STRONG),             # +1: c.2353-2499
    (CT.FRAMESHIFT, "c.2499_2500insA", S.STRONG),
    (CT.FRAMESHIFT, "c.2520dup", S.MODERATE),           # +1: c.2500-2646
    (CT.FRAMESHIFT, "c.2504_2505del", S.MODERATE),      # 2-nt deletion = +1 frame
    (CT.FRAMESHIFT, "c.2440_2441delinsT", S.STRONG),    # net -1
])
def test_cdh1_pvs1_truncation_ranges(cq, hgvs, expected):
    from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1
    assert evaluate_vcep_pvs1(_cdh1_pc(cq, hgvs))[0] == expected


def test_trunc_kind_parsing():
    from acmg_classifier.pvs1.vcep_pvs1 import _trunc_kind_and_start
    k = lambda cq, h: _trunc_kind_and_start(_cdh1_pc(cq, h))
    assert k(CT.FRAMESHIFT, "c.100_101insAT") == ("fs-1", 100)
    assert k(CT.FRAMESHIFT, "c.100_103del") == ("fs-1", 100)
    assert k(CT.FRAMESHIFT, "c.100_101dup") == ("fs-1", 100)
    assert k(CT.FRAMESHIFT, "c.2388+1del") == (None, None)
    assert k(CT.STOP_GAINED, "c.12G>A") == ("nonsense", 12)


# ---------------------------------------------------------------------------
# VCEP exceptions: Likely benign on a single BS (CDH1 v3.1, IDUA, RUNX1, DDX41)
# ---------------------------------------------------------------------------

def _lb_r(crit, strength=None):
    from acmg_classifier.models.criteria import CriteriaResult
    from acmg_classifier.models.enums import ACMGCriterion
    c = getattr(ACMGCriterion, crit)
    return CriteriaResult.met(c, strength) if strength else CriteriaResult.met(c)


@pytest.mark.parametrize("gene,crit,expected", [
    ("CDH1", "BS1", "LIKELY_BENIGN"),
    ("CDH1", "BS2", "LIKELY_BENIGN"),
    ("IDUA", "BS1", "LIKELY_BENIGN"),
    ("IDUA", "BS2", "VUS"),
    ("RUNX1", "BS1", "LIKELY_BENIGN"),
    ("DDX41", "BS1", "LIKELY_BENIGN"),
    ("BRCA1", "BS1", "VUS"),
    (None, "BS1", "VUS"),
])
def test_lb_on_single_bs(gene, crit, expected):
    from acmg_classifier.classification.classifier_2015 import Classifier2015
    from acmg_classifier.models.enums import Pathogenicity
    path, _ = Classifier2015().classify([_lb_r(crit)], gene=gene)
    assert path == getattr(Pathogenicity, expected)


def test_lb_on_single_bs_contradictory_evidence():
    from acmg_classifier.classification.classifier_2015 import Classifier2015
    from acmg_classifier.models.enums import CriterionStrength as S, Pathogenicity
    clf = Classifier2015()
    pp = _lb_r("PP3", S.SUPPORTING)
    # RUNX1/DDX41: "if there is no contradictory evidence supporting pathogenicity"
    assert clf.classify([_lb_r("BS1"), pp], gene="RUNX1")[0] == Pathogenicity.VUS
    # CDH1: no such clause -> Table 5 conflict rule only (no pathogenic rule met)
    assert clf.classify([_lb_r("BS1"), pp], gene="CDH1")[0] == Pathogenicity.LIKELY_BENIGN
    # CDH1 BS1 at Supporting is a BP-level item, not "BS1 alone"
    assert clf.classify([_lb_r("BS1", S.SUPPORTING)], gene="CDH1")[0] == Pathogenicity.VUS
    # pathogenic rule met -> VUS (Table 5 contradiction)
    lp = [_lb_r("PVS1", S.VERY_STRONG), _lb_r("PM2", S.MODERATE), _lb_r("BS2")]
    assert clf.classify(lp, gene="CDH1")[0] == Pathogenicity.VUS
