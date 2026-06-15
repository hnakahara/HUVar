"""Gene-specific VCEP PVS1 trees (RPE65, CYP1B1, VHL, GCK, RAG1, ATM, GP9,
IDUA, ACVRL1): codon-range truncation gates + initiation-codon / splice /
deletion strength overrides."""
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pvs1 import PVS1Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import (
    Assembly, ACMGCriterion, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1

S = CriterionStrength


def _pc(gene, consequence, protein_position=None, exon=None, intron=None):
    return ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=consequence, biotype="protein_coding",
        protein_position=protein_position, exon=exon, intron=intron,
    )


class TestTruncationBands:
    def test_rpe65_main_body_very_strong(self):
        s, _ = evaluate_vcep_pvs1(_pc("RPE65", ConsequenceType.STOP_GAINED, 300))
        assert s == S.VERY_STRONG

    def test_rpe65_last_residues_strong(self):
        s, _ = evaluate_vcep_pvs1(_pc("RPE65", ConsequenceType.STOP_GAINED, 530))
        assert s == S.STRONG

    def test_cyp1b1_haem_domain_through_493_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("CYP1B1", ConsequenceType.FRAMESHIFT, 493))[0] == S.VERY_STRONG

    def test_cyp1b1_after_haem_domain_moderate(self):
        assert evaluate_vcep_pvs1(_pc("CYP1B1", ConsequenceType.STOP_GAINED, 494))[0] == S.MODERATE

    def test_vhl_before_codon54_na(self):
        s, ev = evaluate_vcep_pvs1(_pc("VHL", ConsequenceType.FRAMESHIFT, 30))
        assert s == S.NOT_MET
        assert "outside" in ev

    def test_vhl_critical_domain_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("VHL", ConsequenceType.STOP_GAINED, 100))[0] == S.VERY_STRONG

    def test_vhl_after_2nd_beta_moderate(self):
        assert evaluate_vcep_pvs1(_pc("VHL", ConsequenceType.STOP_GAINED, 210))[0] == S.MODERATE

    def test_gck_exon10_ptc_very_strong(self):
        # Last-exon PTC that escapes NMD is still Very Strong for GCK.
        assert evaluate_vcep_pvs1(_pc("GCK", ConsequenceType.STOP_GAINED, 460))[0] == S.VERY_STRONG

    def test_rag1_core_domain_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("RAG1", ConsequenceType.STOP_GAINED, 1011))[0] == S.VERY_STRONG

    def test_rag1_after_core_na(self):
        assert evaluate_vcep_pvs1(_pc("RAG1", ConsequenceType.FRAMESHIFT, 1012))[0] == S.NOT_MET

    def test_atm_through_r3047_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("ATM", ConsequenceType.STOP_GAINED, 3047))[0] == S.VERY_STRONG

    def test_atm_after_r3047_na(self):
        assert evaluate_vcep_pvs1(_pc("ATM", ConsequenceType.STOP_GAINED, 3050))[0] == S.NOT_MET

    def test_gp9_tm_domain_strong(self):
        assert evaluate_vcep_pvs1(_pc("GP9", ConsequenceType.STOP_GAINED, 150))[0] == S.STRONG

    def test_gp9_after_tm_moderate(self):
        assert evaluate_vcep_pvs1(_pc("GP9", ConsequenceType.STOP_GAINED, 175))[0] == S.MODERATE

    def test_idua_before_c1778_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("IDUA", ConsequenceType.STOP_GAINED, 500))[0] == S.VERY_STRONG

    def test_idua_after_c1778_moderate(self):
        assert evaluate_vcep_pvs1(_pc("IDUA", ConsequenceType.STOP_GAINED, 600))[0] == S.MODERATE

    def test_acvrl1_nmd_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("ACVRL1", ConsequenceType.STOP_GAINED, 400))[0] == S.VERY_STRONG

    def test_acvrl1_critical_region_strong(self):
        assert evaluate_vcep_pvs1(_pc("ACVRL1", ConsequenceType.FRAMESHIFT, 480))[0] == S.STRONG

    def test_acvrl1_after_490_moderate(self):
        assert evaluate_vcep_pvs1(_pc("ACVRL1", ConsequenceType.STOP_GAINED, 500))[0] == S.MODERATE

    def test_unknown_codon_defers(self):
        assert evaluate_vcep_pvs1(_pc("RPE65", ConsequenceType.FRAMESHIFT, None)) is None

    def test_non_vcep_gene_defers(self):
        assert evaluate_vcep_pvs1(_pc("BRCA1", ConsequenceType.STOP_GAINED, 100)) is None

    # --- batch 2 genes ---------------------------------------------------
    def test_pah_upstream_c1285_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("PAH", ConsequenceType.STOP_GAINED, 428))[0] == S.VERY_STRONG

    def test_pah_downstream_c1285_strong(self):
        assert evaluate_vcep_pvs1(_pc("PAH", ConsequenceType.FRAMESHIFT, 429))[0] == S.STRONG

    def test_hnf1a_5prime_c1768_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("HNF1A", ConsequenceType.STOP_GAINED, 589))[0] == S.VERY_STRONG

    def test_hnf1a_nonsense_vs_frameshift_differ(self):
        # codon 610: nonsense > p.601 → Supporting; frameshift ≤ p.618 → Strong.
        assert evaluate_vcep_pvs1(_pc("HNF1A", ConsequenceType.STOP_GAINED, 610))[0] == S.SUPPORTING
        assert evaluate_vcep_pvs1(_pc("HNF1A", ConsequenceType.FRAMESHIFT, 610))[0] == S.STRONG

    def test_hnf1a_frameshift_distal_supporting(self):
        assert evaluate_vcep_pvs1(_pc("HNF1A", ConsequenceType.FRAMESHIFT, 625))[0] == S.SUPPORTING

    def test_gjb2_over_10pct_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("GJB2", ConsequenceType.FRAMESHIFT, 12))[0] == S.VERY_STRONG

    def test_gjb2_last_10pct_moderate(self):
        assert evaluate_vcep_pvs1(_pc("GJB2", ConsequenceType.STOP_GAINED, 220))[0] == S.MODERATE

    def test_foxg1_bands(self):
        assert evaluate_vcep_pvs1(_pc("FOXG1", ConsequenceType.STOP_GAINED, 468))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("FOXG1", ConsequenceType.STOP_GAINED, 475))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("FOXG1", ConsequenceType.STOP_GAINED, 485))[0] == S.MODERATE

    def test_dicer1_nmd_cutoff(self):
        assert evaluate_vcep_pvs1(_pc("DICER1", ConsequenceType.STOP_GAINED, 1850))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("DICER1", ConsequenceType.STOP_GAINED, 1851))[0] == S.MODERATE

    def test_palb2_all_truncations_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("PALB2", ConsequenceType.STOP_GAINED, 1183))[0] == S.VERY_STRONG

    def test_fbn1_nmd_based(self):
        # Mid-gene exon → NMD predicted → Very Strong.
        assert evaluate_vcep_pvs1(_pc("FBN1", ConsequenceType.STOP_GAINED, 1000, exon="30/66"))[0] == S.VERY_STRONG
        # Last exon → NMD escape → Strong (critical C-terminus).
        assert evaluate_vcep_pvs1(_pc("FBN1", ConsequenceType.STOP_GAINED, 2860, exon="66/66"))[0] == S.STRONG

    # --- batch 3 genes ---------------------------------------------------
    def test_gp1ba_tm_strong(self):
        assert evaluate_vcep_pvs1(_pc("GP1BA", ConsequenceType.STOP_GAINED, 540))[0] == S.STRONG

    def test_gp1ba_after_588_moderate(self):
        assert evaluate_vcep_pvs1(_pc("GP1BA", ConsequenceType.FRAMESHIFT, 600))[0] == S.MODERATE

    def test_cdh1_nmd_vs_escape(self):
        assert evaluate_vcep_pvs1(_pc("CDH1", ConsequenceType.STOP_GAINED, 300, exon="3/16"))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("CDH1", ConsequenceType.STOP_GAINED, 836, exon="16/16"))[0] == S.STRONG

    def test_cdh1_splice_default_strong(self):
        # CDH1 caveat: canonical splice default is PVS1_Strong, not Very Strong.
        assert evaluate_vcep_pvs1(_pc("CDH1", ConsequenceType.SPLICE_DONOR, intron="2/15"))[0] == S.STRONG

    def test_aipl1_nonsense_vs_frameshift_bands(self):
        # codon 333: nonsense → Strong (329-346); frameshift → Very Strong (≤337).
        assert evaluate_vcep_pvs1(_pc("AIPL1", ConsequenceType.STOP_GAINED, 333))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("AIPL1", ConsequenceType.FRAMESHIFT, 333))[0] == S.VERY_STRONG

    def test_aipl1_nonsense_distal_moderate(self):
        assert evaluate_vcep_pvs1(_pc("AIPL1", ConsequenceType.STOP_GAINED, 360))[0] == S.MODERATE

    def test_acadvl_nmd_and_escape_10pct(self):
        # Mid-gene → NMD → Very Strong.
        assert evaluate_vcep_pvs1(_pc("ACADVL", ConsequenceType.STOP_GAINED, 200, exon="5/20"))[0] == S.VERY_STRONG
        # Last exon, <10% removed → Moderate.
        assert evaluate_vcep_pvs1(_pc("ACADVL", ConsequenceType.STOP_GAINED, 640, exon="20/20"))[0] == S.MODERATE

    def test_acadvl_intron8_donor_excluded(self):
        s, ev = evaluate_vcep_pvs1(_pc("ACADVL", ConsequenceType.SPLICE_DONOR, intron="8/19"))
        assert s == S.NOT_MET
        assert "intron 8" in ev

    def test_acadvl_other_donor_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("ACADVL", ConsequenceType.SPLICE_DONOR, intron="5/19"))[0] == S.VERY_STRONG

    def test_acadvl_acceptor_not_excluded(self):
        # Exclusion is donor-only; an acceptor (even near intron 8) still fires.
        assert evaluate_vcep_pvs1(_pc("ACADVL", ConsequenceType.SPLICE_ACCEPTOR, intron="8/19"))[0] == S.VERY_STRONG

    def test_tp53_bands(self):
        assert evaluate_vcep_pvs1(_pc("TP53", ConsequenceType.STOP_GAINED, 350))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("TP53", ConsequenceType.STOP_GAINED, 353))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("TP53", ConsequenceType.FRAMESHIFT, 380))[0] == S.MODERATE

    def test_gaa_codon916_cutoff(self):
        assert evaluate_vcep_pvs1(_pc("GAA", ConsequenceType.STOP_GAINED, 915))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("GAA", ConsequenceType.STOP_GAINED, 916))[0] == S.MODERATE

    def test_gamt_nmd_and_escape(self):
        assert evaluate_vcep_pvs1(_pc("GAMT", ConsequenceType.FRAMESHIFT, 100, exon="3/6"))[0] == S.VERY_STRONG
        # Last exon, removes <10% → Moderate.
        assert evaluate_vcep_pvs1(_pc("GAMT", ConsequenceType.STOP_GAINED, 230, exon="6/6"))[0] == S.MODERATE

    # --- batch 4 genes ---------------------------------------------------
    def test_hnf4a_codon419_cutoff(self):
        assert evaluate_vcep_pvs1(_pc("HNF4A", ConsequenceType.STOP_GAINED, 419))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("HNF4A", ConsequenceType.FRAMESHIFT, 420))[0] == S.SUPPORTING

    def test_runx1_nmd_and_escape(self):
        assert evaluate_vcep_pvs1(_pc("RUNX1", ConsequenceType.STOP_GAINED, 100, exon="4/8"))[0] == S.VERY_STRONG
        # C-terminal last exon → NMD escape → Strong.
        assert evaluate_vcep_pvs1(_pc("RUNX1", ConsequenceType.FRAMESHIFT, 440, exon="8/8"))[0] == S.STRONG

    def test_cdkl5_r948_cutoff(self):
        assert evaluate_vcep_pvs1(_pc("CDKL5", ConsequenceType.STOP_GAINED, 948))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("CDKL5", ConsequenceType.STOP_GAINED, 955))[0] == S.MODERATE

    def test_rpgr_orf15_all_very_strong(self):
        # ORF15 isoform (NM_001034853.2): NMD truncations and ORF15
        # glutamylation-disrupting truncations are all PVS1.
        assert evaluate_vcep_pvs1(_pc("RPGR", ConsequenceType.STOP_GAINED, 300, exon="8/15"))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("RPGR", ConsequenceType.FRAMESHIFT, 900, exon="15/15"))[0] == S.VERY_STRONG

    # --- batch 5 genes ---------------------------------------------------
    def test_il2rg_all_truncations_very_strong(self):
        # Last-exon escape still in TM/cytoplasmic domain → Very Strong.
        assert evaluate_vcep_pvs1(_pc("IL2RG", ConsequenceType.STOP_GAINED, 100, exon="3/8"))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("IL2RG", ConsequenceType.FRAMESHIFT, 360, exon="8/8"))[0] == S.VERY_STRONG

    def test_mecp2_e484_cutoff(self):
        # MeCP2_e1 (MANE) numbering: e2 p.E472 -> e1 p.484.
        assert evaluate_vcep_pvs1(_pc("MECP2", ConsequenceType.STOP_GAINED, 484))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("MECP2", ConsequenceType.FRAMESHIFT, 490))[0] == S.MODERATE

    def test_f9_all_truncations_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("F9", ConsequenceType.STOP_GAINED, 50, exon="2/8"))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("F9", ConsequenceType.FRAMESHIFT, 455, exon="8/8"))[0] == S.VERY_STRONG

    def test_abcd1_nmd_and_escape_10pct(self):
        assert evaluate_vcep_pvs1(_pc("ABCD1", ConsequenceType.STOP_GAINED, 300, exon="3/10"))[0] == S.VERY_STRONG
        # Last exon, >10% removed → Strong.
        assert evaluate_vcep_pvs1(_pc("ABCD1", ConsequenceType.STOP_GAINED, 655, exon="10/10"))[0] == S.STRONG
        # Last exon, <10% removed → Moderate.
        assert evaluate_vcep_pvs1(_pc("ABCD1", ConsequenceType.STOP_GAINED, 720, exon="10/10"))[0] == S.MODERATE

    # --- batch 6 genes (cspec re-examination) ----------------------------
    def test_ada_nmd_escape_10pct(self):
        assert evaluate_vcep_pvs1(_pc("ADA", ConsequenceType.STOP_GAINED, 100, exon="4/12"))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("ADA", ConsequenceType.STOP_GAINED, 360, exon="12/12"))[0] == S.MODERATE

    def test_il7r_escape_strong(self):
        # NMD escape (TM domain begins aa240) → Strong.
        assert evaluate_vcep_pvs1(_pc("IL7R", ConsequenceType.STOP_GAINED, 100, exon="3/8"))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("IL7R", ConsequenceType.FRAMESHIFT, 450, exon="8/8"))[0] == S.STRONG

    def test_foxn1_bands(self):
        assert evaluate_vcep_pvs1(_pc("FOXN1", ConsequenceType.STOP_GAINED, 400))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("FOXN1", ConsequenceType.STOP_GAINED, 540))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("FOXN1", ConsequenceType.STOP_GAINED, 600))[0] == S.MODERATE

    def test_rag2_single_exon_critical(self):
        assert evaluate_vcep_pvs1(_pc("RAG2", ConsequenceType.STOP_GAINED, 400))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("RAG2", ConsequenceType.FRAMESHIFT, 500))[0] == S.NOT_MET

    def test_ctla4_exon3_bands(self):
        assert evaluate_vcep_pvs1(_pc("CTLA4", ConsequenceType.STOP_GAINED, 172))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("CTLA4", ConsequenceType.STOP_GAINED, 180))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("CTLA4", ConsequenceType.STOP_GAINED, 210))[0] == S.MODERATE

    def test_kcnq1_bands(self):
        assert evaluate_vcep_pvs1(_pc("KCNQ1", ConsequenceType.STOP_GAINED, 581))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("KCNQ1", ConsequenceType.STOP_GAINED, 600))[0] == S.MODERATE
        assert evaluate_vcep_pvs1(_pc("KCNQ1", ConsequenceType.STOP_GAINED, 650))[0] == S.SUPPORTING

    def test_lynch_cutoffs(self):
        assert evaluate_vcep_pvs1(_pc("MLH1", ConsequenceType.STOP_GAINED, 753))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("MLH1", ConsequenceType.STOP_GAINED, 755))[0] == S.MODERATE
        assert evaluate_vcep_pvs1(_pc("MSH2", ConsequenceType.FRAMESHIFT, 900))[0] == S.MODERATE
        assert evaluate_vcep_pvs1(_pc("MSH6", ConsequenceType.STOP_GAINED, 1341))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("PMS2", ConsequenceType.STOP_GAINED, 800))[0] == S.MODERATE

    def test_otc_cterminal_strong(self):
        assert evaluate_vcep_pvs1(_pc("OTC", ConsequenceType.STOP_GAINED, 300))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("OTC", ConsequenceType.FRAMESHIFT, 350))[0] == S.STRONG

    def test_slc9a6_offset_bands(self):
        assert evaluate_vcep_pvs1(_pc("SLC9A6", ConsequenceType.STOP_GAINED, 573))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("SLC9A6", ConsequenceType.STOP_GAINED, 600))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("SLC9A6", ConsequenceType.STOP_GAINED, 650))[0] == S.MODERATE

    def test_tcf4_cutoff(self):
        assert evaluate_vcep_pvs1(_pc("TCF4", ConsequenceType.STOP_GAINED, 643))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("TCF4", ConsequenceType.FRAMESHIFT, 660))[0] == S.MODERATE

    def test_ube3a_offset_bands(self):
        assert evaluate_vcep_pvs1(_pc("UBE3A", ConsequenceType.STOP_GAINED, 861))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("UBE3A", ConsequenceType.STOP_GAINED, 865))[0] == S.STRONG

    def test_gucy2d_lca_bands(self):
        assert evaluate_vcep_pvs1(_pc("GUCY2D", ConsequenceType.STOP_GAINED, 1000))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("GUCY2D", ConsequenceType.STOP_GAINED, 1080))[0] == S.STRONG

    def test_rs1_bands(self):
        assert evaluate_vcep_pvs1(_pc("RS1", ConsequenceType.STOP_GAINED, 200))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("RS1", ConsequenceType.STOP_GAINED, 224))[0] == S.STRONG

    # --- batch 7 genes (supplied decision-tree files) --------------------
    def test_eng_bands(self):
        assert evaluate_vcep_pvs1(_pc("ENG", ConsequenceType.STOP_GAINED, 601))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("ENG", ConsequenceType.FRAMESHIFT, 620))[0] == S.MODERATE

    def test_gp1bb_bands(self):
        assert evaluate_vcep_pvs1(_pc("GP1BB", ConsequenceType.STOP_GAINED, 160))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("GP1BB", ConsequenceType.STOP_GAINED, 195))[0] == S.MODERATE
        assert evaluate_vcep_pvs1(_pc("GP1BB", ConsequenceType.SPLICE_DONOR, intron="1/1")) is None

    def test_scn1b_bands(self):
        assert evaluate_vcep_pvs1(_pc("SCN1B", ConsequenceType.STOP_GAINED, 204))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("SCN1B", ConsequenceType.STOP_GAINED, 210))[0] == S.MODERATE

    def test_scn2a_bands(self):
        assert evaluate_vcep_pvs1(_pc("SCN2A", ConsequenceType.STOP_GAINED, 1591))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("SCN2A", ConsequenceType.STOP_GAINED, 1700))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("SCN2A", ConsequenceType.STOP_GAINED, 1900))[0] == S.MODERATE

    def test_scn8a_bands(self):
        assert evaluate_vcep_pvs1(_pc("SCN8A", ConsequenceType.STOP_GAINED, 1582))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("SCN8A", ConsequenceType.FRAMESHIFT, 1700))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("SCN8A", ConsequenceType.STOP_GAINED, 1900))[0] == S.MODERATE

    def test_neb_bands(self):
        assert evaluate_vcep_pvs1(_pc("NEB", ConsequenceType.STOP_GAINED, 5000))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("NEB", ConsequenceType.STOP_GAINED, 8500))[0] == S.MODERATE

    def test_f8_bands(self):
        assert evaluate_vcep_pvs1(_pc("F8", ConsequenceType.STOP_GAINED, 2000))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("F8", ConsequenceType.FRAMESHIFT, 2300))[0] == S.STRONG

    def test_pten_bands(self):
        assert evaluate_vcep_pvs1(_pc("PTEN", ConsequenceType.STOP_GAINED, 375))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("PTEN", ConsequenceType.FRAMESHIFT, 390))[0] == S.MODERATE
        assert evaluate_vcep_pvs1(_pc("PTEN", ConsequenceType.START_LOST))[0] == S.VERY_STRONG

    def test_mybpc3_bands(self):
        assert evaluate_vcep_pvs1(_pc("MYBPC3", ConsequenceType.STOP_GAINED, 1253))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("MYBPC3", ConsequenceType.FRAMESHIFT, 1260))[0] == S.MODERATE

    def test_hbb_bands(self):
        # Early PTC escapes NMD but removes >10% → Strong.
        assert evaluate_vcep_pvs1(_pc("HBB", ConsequenceType.STOP_GAINED, 10))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("HBB", ConsequenceType.STOP_GAINED, 50))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("HBB", ConsequenceType.STOP_GAINED, 100))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("HBB", ConsequenceType.STOP_GAINED, 140))[0] == S.MODERATE

    def test_hba2_bands(self):
        assert evaluate_vcep_pvs1(_pc("HBA2", ConsequenceType.STOP_GAINED, 80))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("HBA2", ConsequenceType.STOP_GAINED, 100))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("HBA2", ConsequenceType.STOP_GAINED, 135))[0] == S.MODERATE

    def test_eng_start_strong(self):
        assert evaluate_vcep_pvs1(_pc("ENG", ConsequenceType.START_LOST))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("HBB", ConsequenceType.START_LOST))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("F8", ConsequenceType.START_LOST))[0] == S.MODERATE


class TestStartLost:
    def test_rpe65_start_lost_strong(self):
        assert evaluate_vcep_pvs1(_pc("RPE65", ConsequenceType.START_LOST))[0] == S.STRONG

    def test_gck_start_lost_supporting(self):
        assert evaluate_vcep_pvs1(_pc("GCK", ConsequenceType.START_LOST))[0] == S.SUPPORTING

    def test_vhl_start_lost_na(self):
        # Met1 loss is rescued by the downstream p19 start at Met54.
        assert evaluate_vcep_pvs1(_pc("VHL", ConsequenceType.START_LOST))[0] == S.NOT_MET

    def test_rag1_start_lost_defers(self):
        assert evaluate_vcep_pvs1(_pc("RAG1", ConsequenceType.START_LOST)) is None

    def test_pah_start_lost_strong(self):
        assert evaluate_vcep_pvs1(_pc("PAH", ConsequenceType.START_LOST))[0] == S.STRONG

    def test_hnf1a_start_lost_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("HNF1A", ConsequenceType.START_LOST))[0] == S.VERY_STRONG

    def test_foxg1_start_lost_supporting(self):
        assert evaluate_vcep_pvs1(_pc("FOXG1", ConsequenceType.START_LOST))[0] == S.SUPPORTING

    def test_dicer1_start_lost_na(self):
        assert evaluate_vcep_pvs1(_pc("DICER1", ConsequenceType.START_LOST))[0] == S.NOT_MET

    def test_fbn1_start_lost_moderate(self):
        assert evaluate_vcep_pvs1(_pc("FBN1", ConsequenceType.START_LOST))[0] == S.MODERATE

    def test_aipl1_start_lost_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("AIPL1", ConsequenceType.START_LOST))[0] == S.VERY_STRONG

    def test_acadvl_start_lost_strong(self):
        assert evaluate_vcep_pvs1(_pc("ACADVL", ConsequenceType.START_LOST))[0] == S.STRONG

    def test_tp53_start_lost_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("TP53", ConsequenceType.START_LOST))[0] == S.VERY_STRONG

    def test_gaa_start_lost_strong(self):
        assert evaluate_vcep_pvs1(_pc("GAA", ConsequenceType.START_LOST))[0] == S.STRONG

    def test_gamt_start_lost_moderate(self):
        assert evaluate_vcep_pvs1(_pc("GAMT", ConsequenceType.START_LOST))[0] == S.MODERATE

    def test_gp1ba_start_lost_moderate(self):
        assert evaluate_vcep_pvs1(_pc("GP1BA", ConsequenceType.START_LOST))[0] == S.MODERATE

    def test_hnf4a_start_lost_strong(self):
        assert evaluate_vcep_pvs1(_pc("HNF4A", ConsequenceType.START_LOST))[0] == S.STRONG

    def test_runx1_start_lost_defers(self):
        assert evaluate_vcep_pvs1(_pc("RUNX1", ConsequenceType.START_LOST)) is None

    def test_cdkl5_start_lost_supporting(self):
        assert evaluate_vcep_pvs1(_pc("CDKL5", ConsequenceType.START_LOST))[0] == S.SUPPORTING

    def test_rpgr_start_lost_moderate(self):
        assert evaluate_vcep_pvs1(_pc("RPGR", ConsequenceType.START_LOST))[0] == S.MODERATE

    def test_mecp2_start_lost_na(self):
        # Initiation codon N/A (MECP2_E1 alternate isoform/start).
        assert evaluate_vcep_pvs1(_pc("MECP2", ConsequenceType.START_LOST))[0] == S.NOT_MET

    def test_f9_start_lost_supporting(self):
        assert evaluate_vcep_pvs1(_pc("F9", ConsequenceType.START_LOST))[0] == S.SUPPORTING

    def test_abcd1_start_lost_moderate(self):
        assert evaluate_vcep_pvs1(_pc("ABCD1", ConsequenceType.START_LOST))[0] == S.MODERATE

    def test_il2rg_start_lost_defers(self):
        assert evaluate_vcep_pvs1(_pc("IL2RG", ConsequenceType.START_LOST)) is None

    def test_batch6_start_lost(self):
        assert evaluate_vcep_pvs1(_pc("CTLA4", ConsequenceType.START_LOST))[0] == S.SUPPORTING
        assert evaluate_vcep_pvs1(_pc("MSH6", ConsequenceType.START_LOST))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("PMS2", ConsequenceType.START_LOST))[0] == S.STRONG
        assert evaluate_vcep_pvs1(_pc("TCF4", ConsequenceType.START_LOST))[0] == S.SUPPORTING
        assert evaluate_vcep_pvs1(_pc("SLC9A6", ConsequenceType.START_LOST))[0] == S.SUPPORTING
        assert evaluate_vcep_pvs1(_pc("UBE3A", ConsequenceType.START_LOST))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("GUCY2D", ConsequenceType.START_LOST))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("RS1", ConsequenceType.START_LOST))[0] == S.VERY_STRONG
        assert evaluate_vcep_pvs1(_pc("MLH1", ConsequenceType.START_LOST)) is None


class TestSpliceAndDeletion:
    def test_rpe65_canonical_splice_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("RPE65", ConsequenceType.SPLICE_DONOR))[0] == S.VERY_STRONG

    def test_gp9_splice_defers(self):
        assert evaluate_vcep_pvs1(_pc("GP9", ConsequenceType.SPLICE_ACCEPTOR)) is None

    def test_vhl_exon_deletion_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("VHL", ConsequenceType.TRANSCRIPT_ABLATION))[0] == S.VERY_STRONG

    def test_rag1_full_gene_deletion_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("RAG1", ConsequenceType.TRANSCRIPT_ABLATION))[0] == S.VERY_STRONG


# --- End-to-end through the evaluator (confirms caps are bypassed) -----------

def _ev():
    cfg = MagicMock()
    cfg.disease_prevalence_tsv.exists.return_value = False
    return PVS1Evaluator(cfg)


def _ann(pc):
    # loeuf high + no ClinVar P/LP nulls would normally cap PVS1 to Moderate;
    # the VCEP handler must bypass that.
    return AnnotationData(consequences=[pc], gnomad=GnomADData(loeuf=0.5))


def _snv():
    return VariantRecord(chrom="chr3", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)


class TestEvaluatorRouting:
    def test_rag1_single_exon_truncation_fires_very_strong(self):
        # Single-exon gene: generic tree never predicts NMD and would withhold
        # PVS1 absent a domain. VCEP handler grants Very Strong.
        r = _ev().evaluate(_snv(), _ann(_pc("RAG1", ConsequenceType.STOP_GAINED, 500)))
        assert r.triggered
        assert r.strength == S.VERY_STRONG
        assert r.criterion == ACMGCriterion.PVS1

    def test_vhl_early_truncation_withheld(self):
        r = _ev().evaluate(_snv(), _ann(_pc("VHL", ConsequenceType.FRAMESHIFT, 30)))
        assert not r.triggered

    def test_gck_start_lost_supporting_routes(self):
        r = _ev().evaluate(_snv(), _ann(_pc("GCK", ConsequenceType.START_LOST)))
        assert r.triggered
        assert r.strength == S.SUPPORTING

    def test_vhl_start_lost_withheld(self):
        r = _ev().evaluate(_snv(), _ann(_pc("VHL", ConsequenceType.START_LOST)))
        assert not r.triggered

    def test_non_vcep_missense_not_lof(self):
        r = _ev().evaluate(_snv(), _ann(_pc("RPE65", ConsequenceType.MISSENSE)))
        assert not r.triggered
        assert "Not a LoF" in r.evidence

    def test_mecp2_start_lost_withheld(self):
        r = _ev().evaluate(_snv(), _ann(_pc("MECP2", ConsequenceType.START_LOST)))
        assert not r.triggered

    def test_f9_last_exon_truncation_very_strong(self):
        # Single-exon-escape would normally cap; F9 keeps Very Strong.
        r = _ev().evaluate(_snv(), _ann(_pc("F9", ConsequenceType.FRAMESHIFT, 455, exon="8/8")))
        assert r.triggered
        assert r.strength == S.VERY_STRONG
