"""PM1 corrections after re-checking each rule against the ClinGen CSpec text
(Genome Medicine revision): committed regions, transcript renumbering
(MECP2, UBE3A), exon-defined regions (F9, BRAF, RAF1), and the VCEP caveats
(LDLR/KCNQ1 need PM2; RAG1/RAG2 not with BA1/BS1/BS2; PAH not with PP3_Strong;
RYR1 Moderate->Supporting with PS1/PM5; FBN1 PM1_Strong excludes PS1/PM5;
FBN1 N>S / G>A not given PM1)."""
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.pathogenic.pm1 import PM1Evaluator
from acmg_classifier.criteria.pm1_hotspots import PM1Hotspots
from acmg_classifier.criteria.registry import _apply_pm1_co_requirements
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion, Assembly, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.variant import VariantRecord

S = CriterionStrength
_TSV = Path(__file__).resolve().parents[2] / "resources" / "shared" / "pm1_hotspots.tsv"


@pytest.fixture(scope="module")
def h():
    return PM1Hotspots(_TSV)


@pytest.mark.parametrize("gene,pos,exp", [
    ("MECP2", 101, None), ("MECP2", 102, S.MODERATE), ("MECP2", 174, S.MODERATE),
    ("MECP2", 175, None), ("MECP2", 314, S.MODERATE), ("MECP2", 305, None),
    ("UBE3A", 840, S.MODERATE), ("UBE3A", 820, None),
    ("GP1BB", 147, S.MODERATE),
    ("F9", 85, S.MODERATE), ("F9", 174, S.MODERATE), ("F9", 175, None),
    ("F9", 267, S.STRONG),
    ("HNF4A", 35, None), ("HNF4A", 37, S.SUPPORTING), ("HNF4A", 43, S.MODERATE),
    ("TNNT2", 85, None), ("TNNT2", 189, S.SUPPORTING),
    ("ENG", 360, None), ("ENG", 363, S.MODERATE),
    ("RPE65", 118, S.MODERATE),
    ("GUCY2D", 874, S.SUPPORTING), ("GUCY2D", 883, S.MODERATE), ("GUCY2D", 960, None),
    ("DICER1", 1750, S.SUPPORTING), ("DICER1", 1705, S.MODERATE),
    ("ACADVL", 218, S.MODERATE), ("ACADVL", 490, S.MODERATE), ("ACADVL", 562, S.MODERATE),
    ("ACADVL", 300, None),
    ("PTPN11", 4, S.MODERATE),
    ("F8", 2283, None), ("F8", 2282, S.MODERATE),
    ("BRAF", 245, S.MODERATE), ("BRAF", 450, S.MODERATE), ("BRAF", 581, None),
    ("RAF1", 491, S.MODERATE), ("RAF1", 644, S.MODERATE), ("RAF1", 550, None),
])
def test_regions(h, gene, pos, exp):
    assert h.lookup(gene, pos) == exp


def _ann(gene, pos, aa="R/C"):
    return AnnotationData(gnomad=GnomADData(), consequences=[ConsequenceInfo(
        transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        is_mane_select=True, protein_position=pos, amino_acids=aa,
    )])


def _res(pm1=S.MODERATE, **others):
    out = [CriteriaResult.met(ACMGCriterion.PM1, strength=pm1, evidence="x")]
    for code, st in others.items():
        out.append(CriteriaResult.met(ACMGCriterion(code), strength=st, evidence="y"))
    return out


class TestCoRequirements:
    def test_ldlr_requires_pm2(self):
        r = _res()
        _apply_pm1_co_requirements(r, _ann("LDLR", 150))
        assert r[0].suppressed
        r = _res(PM2=S.MODERATE)
        _apply_pm1_co_requirements(r, _ann("LDLR", 150))
        assert not r[0].suppressed

    def test_kcnq1_requires_pm2(self):
        r = _res(PM2=S.SUPPORTING)
        _apply_pm1_co_requirements(r, _ann("KCNQ1", 310))
        assert not r[0].suppressed

    @pytest.mark.parametrize("gene", ["RAG1", "RAG2"])
    def test_scid_not_with_bs1(self, gene):
        r = _res(BS1=S.STRONG)
        _apply_pm1_co_requirements(r, _ann(gene, 420))
        assert r[0].suppressed and "BS1" in r[0].evidence

    def test_pah_not_with_pp3_strong(self):
        r = _res(PP3=S.STRONG)
        _apply_pm1_co_requirements(r, _ann("PAH", 158))
        assert r[0].suppressed and "PP3_Strong" in r[0].evidence
        r = _res(PP3=S.MODERATE)
        _apply_pm1_co_requirements(r, _ann("PAH", 158))
        assert not r[0].suppressed

    def test_ryr1_downgrade_with_pm5(self):
        r = _res(PM5=S.MODERATE)
        _apply_pm1_co_requirements(r, _ann("RYR1", 330))
        assert r[0].strength == S.SUPPORTING and not r[0].suppressed
        r = _res()
        _apply_pm1_co_requirements(r, _ann("RYR1", 330))
        assert r[0].strength == S.MODERATE
        # pore region (Congenital Myopathies VCEP) is not downgraded
        r = _res(PS1=S.STRONG)
        _apply_pm1_co_requirements(r, _ann("RYR1", 4861))
        assert r[0].strength == S.MODERATE

    def test_fbn1_strong_excludes_ps1_pm5(self):
        r = _res(pm1=S.STRONG, PS1=S.STRONG, PM5=S.MODERATE)
        _apply_pm1_co_requirements(r, _ann("FBN1", 250))
        assert r[1].suppressed and r[2].suppressed and not r[0].suppressed
        r = _res(pm1=S.MODERATE, PM5=S.MODERATE)
        _apply_pm1_co_requirements(r, _ann("FBN1", 85))
        assert not r[1].suppressed


class TestFbn1Tolerated:
    def _eval(self, pos, aa):
        cfg = MagicMock()
        cfg.pm1_hotspots_tsv = _TSV
        return PM1Evaluator(cfg).evaluate(
            VariantRecord(chrom="chr15", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38),
            _ann("FBN1", pos, aa),
        )

    def test_n_to_s_not_applied(self):
        r = self._eval(1489, "N/S")
        assert not r.triggered and "tolerated" in r.evidence

    def test_other_change_at_same_residue_applied(self):
        assert self._eval(1489, "N/K").triggered


@pytest.mark.parametrize("pos,alt,exp", [
    (220, "C", S.MODERATE),      # Y220C: 127 samples
    (220, "S", S.SUPPORTING),    # Y220S: 8
    (179, "R", S.MODERATE),      # H179R: 63
    (344, "P", S.SUPPORTING),    # L344P: 3
    (273, "W", S.MODERATE),      # VCEP codon list: any missense at 273
    (220, "W", None),            # change not seen in cancerhotspots
    (100, "A", None),            # residue not a hotspot
])
def test_tp53_cancerhotspots(h, pos, alt, exp):
    assert h.lookup("TP53", pos, alt) == exp


@pytest.mark.parametrize("pos,alt,exp", [
    (88, "C", S.MODERATE),     # W88 struck from the somatic list -> germline (Chiorean)
    (178, "P", S.MODERATE),    # L178 (Stebbins)
    (80, "R", S.MODERATE),     # S80 germline: any missense, cancerhotspots not needed
    (89, "H", S.SUPPORTING),   # L89H somatic, 9 samples (<10)
    (68, "P", S.SUPPORTING),   # S68P somatic, 1 sample
    (68, "L", None),           # change not seen in cancerhotspots
    (188, "P", None),          # in cancerhotspots but not in the VCEP table
    (64, "P", None),
])
def test_vhl_hotspot_table(h, pos, alt, exp):
    assert h.lookup("VHL", pos, alt) == exp
