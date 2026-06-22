"""User-supplied VCEP data batch: GALT PVS1, VHL PM1/PM4, LDLR PM1.

- GALT PVS1 (GN158 decision tree): NMD codon bands, initiation Strong, in-frame
  exon 6/7/9/10 splice downgrades.
- VHL PM1 (GN078): germline + somatic hotspot residues → Moderate.
- LDLR PM1 (GN013): exon 4 (codons 105-232) + 60 conserved Cys residues → Moderate.
- VHL PM4 (GN078): in-frame indel in the beta/alpha domains (63-204) → Moderate;
  outside → N/A; stop-loss → Moderate.
"""
from pathlib import Path

from acmg_classifier.criteria.pm1_hotspots import PM1Hotspots
from acmg_classifier.criteria.pm4_regions import PM4Regions
from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType, CriterionStrength
from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1, _SPECS
from acmg_classifier.pvs1.vcep_pvs1_exons import SpliceExonOverrides

S = CriterionStrength
_ROOT = Path(__file__).resolve().parents[2]
_PM1 = _ROOT / "resources" / "shared" / "pm1_hotspots.tsv"
_PM4 = _ROOT / "resources" / "shared" / "pm4_regions.tsv"
_SPLICE = _ROOT / "resources" / "shared" / "vcep_pvs1_splice_exons.tsv"


def _pc(gene, consequence, protein_position=None, intron=None):
    return ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=consequence, biotype="protein_coding",
        protein_position=protein_position, intron=intron,
    )


class TestGALTpvs1:
    def test_spec_registered(self):
        assert _SPECS["GALT"].transcript == "NM_000155.4"
        assert _SPECS["GALT"].aa_len == 379

    def test_nmd_region_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("GALT", ConsequenceType.STOP_GAINED, 100))[0] == S.VERY_STRONG

    def test_escape_over_10pct_very_strong(self):
        # codon 340 → removes 380-340=40 aa (>38, >10%) → PVS1.
        assert evaluate_vcep_pvs1(_pc("GALT", ConsequenceType.FRAMESHIFT, 340))[0] == S.VERY_STRONG

    def test_escape_under_10pct_moderate(self):
        # codon 350 → removes 30 aa (<=38) → PVS1_Moderate.
        assert evaluate_vcep_pvs1(_pc("GALT", ConsequenceType.STOP_GAINED, 350))[0] == S.MODERATE

    def test_initiation_strong(self):
        assert evaluate_vcep_pvs1(_pc("GALT", ConsequenceType.START_LOST))[0] == S.STRONG

    def test_deletion_very_strong(self):
        assert evaluate_vcep_pvs1(_pc("GALT", ConsequenceType.TRANSCRIPT_ABLATION))[0] == S.VERY_STRONG

    def test_splice_exon6_strong(self):
        ov = SpliceExonOverrides(_SPLICE)
        # donor of intron 6 skips exon 6 (critical active site) → Strong.
        s, ev = evaluate_vcep_pvs1(_pc("GALT", ConsequenceType.SPLICE_DONOR, intron="6/10"), ov)
        assert s == S.STRONG and "exon 6" in ev

    def test_splice_exon7_moderate(self):
        ov = SpliceExonOverrides(_SPLICE)
        assert evaluate_vcep_pvs1(_pc("GALT", ConsequenceType.SPLICE_DONOR, intron="7/10"), ov)[0] == S.MODERATE

    def test_splice_frame_disrupting_exon_flat_very_strong(self):
        ov = SpliceExonOverrides(_SPLICE)
        # exon 1 skip (frame-disrupting) has no override → flat PVS1.
        assert evaluate_vcep_pvs1(_pc("GALT", ConsequenceType.SPLICE_DONOR, intron="1/10"), ov)[0] == S.VERY_STRONG


class TestVHLpm1:
    def _h(self):
        return PM1Hotspots(_PM1)

    def test_germline_residue_moderate(self):
        assert self._h().lookup("VHL", 167) == CriterionStrength.MODERATE   # R167 germline
        assert self._h().lookup("VHL", 65) == CriterionStrength.MODERATE    # bare germline pos

    def test_somatic_residue_moderate(self):
        assert self._h().lookup("VHL", 89) == CriterionStrength.MODERATE    # L89 somatic
        assert self._h().lookup("VHL", 169) == CriterionStrength.MODERATE   # L169 somatic

    def test_non_hotspot_none(self):
        assert self._h().lookup("VHL", 200) is None


class TestLDLRpm1:
    def _h(self):
        return PM1Hotspots(_PM1)

    def test_exon4_region_moderate(self):
        assert self._h().lookup("LDLR", 150) == CriterionStrength.MODERATE   # inside exon 4
        assert self._h().lookup("LDLR", 105) == CriterionStrength.MODERATE   # exon-4 start
        assert self._h().lookup("LDLR", 232) == CriterionStrength.MODERATE   # exon-4 end

    def test_cys_outside_exon4_moderate(self):
        assert self._h().lookup("LDLR", 27) == CriterionStrength.MODERATE    # Cys27 (before exon 4)
        assert self._h().lookup("LDLR", 711) == CriterionStrength.MODERATE   # Cys711 (after exon 4)

    def test_non_hotspot_none(self):
        assert self._h().lookup("LDLR", 500) is None


class TestFBN1pm1:
    def _h(self):
        return PM1Hotspots(_PM1)

    def test_cbegf_cys_strong(self):
        # 250 is a cysteine in a cbEGF (EGF-like calcium-binding) domain.
        assert self._h().lookup("FBN1", 250) == CriterionStrength.STRONG

    def test_egf_tb_hybrid_cys_moderate(self):
        assert self._h().lookup("FBN1", 85) == CriterionStrength.MODERATE   # EGF-like non-cb
        assert self._h().lookup("FBN1", 186) == CriterionStrength.MODERATE  # TB domain

    def test_non_cys_position_none(self):
        # 300 is not a curated cysteine residue.
        assert self._h().lookup("FBN1", 300) is None


class TestFBN1GlyMotif:
    def _h(self):
        return PM1Hotspots(_PM1)

    def test_gly_between_cys_moderate(self):
        # 259 is a critical Gly between Cys2-Cys3 of a cbEGF domain → Moderate.
        assert self._h().lookup("FBN1", 259) == CriterionStrength.MODERATE
        assert self._h().lookup("FBN1", 301) == CriterionStrength.MODERATE

    def test_cbegf_cys_still_strong(self):
        assert self._h().lookup("FBN1", 250) == CriterionStrength.STRONG


class TestFBN1CysCreating:
    def _ev(self):
        from unittest.mock import MagicMock
        from acmg_classifier.criteria.pathogenic.pm1 import PM1Evaluator
        cfg = MagicMock()
        cfg.pm1_hotspots_tsv = _PM1
        cfg.clinvar_sqlite = Path("does/not/exist.sqlite")
        return PM1Evaluator(cfg)

    def _ann(self, pos, amino_acids):
        from acmg_classifier.models.annotation import AnnotationData, GnomADData
        return AnnotationData(gnomad=GnomADData(), consequences=[ConsequenceInfo(
            transcript_id="NM", gene_id="ENSG", gene_symbol="FBN1",
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
            is_mane_select=True, protein_position=pos, amino_acids=amino_acids,
        )])

    def _snv(self):
        from acmg_classifier.models.variant import VariantRecord
        from acmg_classifier.models.enums import Assembly
        return VariantRecord(chrom="chr15", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38)

    def test_cys_creating_in_domain_moderate(self):
        # Residue 90 is in a disulfide domain (81-112) but not a curated residue;
        # a missense to Cys (alt=C) → PM1_Moderate.
        r = self._ev().evaluate(self._snv(), self._ann(90, "S/C"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE
        assert "Cys-creating" in r.evidence

    def test_non_cys_alt_not_met(self):
        # Same position, alt != Cys, and 90 is not a curated hotspot residue.
        r = self._ev().evaluate(self._snv(), self._ann(90, "S/A"))
        assert not r.triggered

    def test_cys_creating_outside_domain_not_met(self):
        # Position 50 is outside all disulfide domains → Cys-creating does not apply.
        r = self._ev().evaluate(self._snv(), self._ann(50, "S/C"))
        assert not r.triggered


class TestRPGRpm4:
    def _r(self):
        return PM4Regions(_PM4)

    def test_exon1_14_and_orf15_moderate(self):
        assert self._r().indel_strength("RPGR", 300) == CriterionStrength.MODERATE   # exon 1-14
        assert self._r().indel_strength("RPGR", 800) == CriterionStrength.MODERATE   # ORF15 585-1078

    def test_orf15_repeat_not_met(self):
        assert self._r().indel_strength("RPGR", 1100) == "not_met"   # repetitive C-terminus

    def test_stoploss_strong(self):
        assert self._r().stoploss_strength("RPGR") == CriterionStrength.STRONG


class TestVHLpm4:
    def _r(self):
        return PM4Regions(_PM4)

    def test_in_domain_moderate(self):
        assert self._r().indel_strength("VHL", 100) == CriterionStrength.MODERATE
        assert self._r().indel_strength("VHL", 200) == CriterionStrength.MODERATE

    def test_outside_domain_not_met(self):
        assert self._r().indel_strength("VHL", 50) == "not_met"   # before beta domain (63)

    def test_stoploss_moderate(self):
        assert self._r().stoploss_strength("VHL") == CriterionStrength.MODERATE
