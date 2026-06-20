"""PM4 region / strength rules (pm4_regions.tsv + PM4Regions + evaluator).

Covers PM4_Strong residues (RUNX1), allow-list regions with an N/A default
(MYOC), Moderate-in-domain / Supporting-outside / denied-repeat (DICER1),
deny regions combined with the size-Supporting tier (MECP2), stop-loss-only
genes (CDH1), and the stop-loss N/A gene (CYP1B1).
"""
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pm4_regions import PM4Regions
from acmg_classifier.criteria.pathogenic.pm4 import PM4Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_REGIONS = (
    "gene_symbol\tstrength\tregions\tresidues\n"
    "RUNX1\tstrong\t\t107,110,201\n"
    "RUNX1\tsupporting\t89-204\t\n"
    "RUNX1\tregion_default\tnot_met\t\n"
    "MYOC\tmoderate\t246-502\t\n"
    "MYOC\tregion_default\tnot_met\t\n"
    "DICER1\tmoderate\t1682-1846\t\n"
    "DICER1\tdeny\t606-609\t\n"
    "DICER1\tregion_default\tsupporting\t\n"
    "MECP2\tdeny\t381-405\t\n"
    "MECP2\tregion_default\tmoderate\t\n"
    "CDH1\tregion_default\tnot_met\t\n"
    "CDH1\tstoploss\tmoderate\t\n"
    "CYP1B1\tstoploss\tnot_applicable\t\n"
)

# pm4_supporting_max_aa: MECP2 <3aa -> Supporting (size tier).
_DP = (
    "gene_symbol\tpm4\tpm4_supporting_max_aa\n"
    "MECP2\tapplicable\t2\n"
)


def _paths(tmp_path):
    reg = tmp_path / "pm4_regions.tsv"
    reg.write_text(_REGIONS, encoding="utf-8")
    dp = tmp_path / "dp.tsv"
    dp.write_text(_DP, encoding="utf-8")
    return reg, dp


def _cfg(tmp_path):
    reg, dp = _paths(tmp_path)
    cfg = MagicMock()
    cfg.pm4_regions_tsv = reg
    cfg.disease_prevalence_tsv = dp
    return cfg


def _ann(gene, consequence, pos=None):
    return AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=consequence, biotype="protein_coding", protein_position=pos,
    )])


def _var(ref, alt):
    return VariantRecord(chrom="chr1", pos=100, ref=ref, alt=alt, assembly=Assembly.GRCH38)


_DEL_1AA = _var("ACGT", "A")     # 1 aa
_DEL_2AA = _var("ACGTACG", "A")  # 2 aa
_DEL_3AA = _var("ACGTACGTAC", "A")  # 3 aa
_INS = _var("A", "AGCT")          # 1 aa insertion


class TestLoader:
    def test_runx1_tiers(self, tmp_path):
        reg, _ = _paths(tmp_path)
        h = PM4Regions(reg)
        assert h.indel_strength("RUNX1", 107) == CriterionStrength.STRONG  # residue
        assert h.indel_strength("RUNX1", 150) == CriterionStrength.SUPPORTING  # 89-204
        assert h.indel_strength("RUNX1", 300) == "not_met"  # outside RHD

    def test_dicer1(self, tmp_path):
        reg, _ = _paths(tmp_path)
        h = PM4Regions(reg)
        assert h.indel_strength("DICER1", 1700) == CriterionStrength.MODERATE
        assert h.indel_strength("DICER1", 607) == "not_met"   # repeat deny
        assert h.indel_strength("DICER1", 100) == CriterionStrength.SUPPORTING  # default

    def test_stoploss(self, tmp_path):
        reg, _ = _paths(tmp_path)
        h = PM4Regions(reg)
        assert h.stoploss_strength("CDH1") == CriterionStrength.MODERATE
        assert h.stoploss_strength("CYP1B1") == "not_applicable"
        assert h.stoploss_strength("RUNX1") is None


class TestEvaluator:
    def _ev(self, tmp_path):
        return PM4Evaluator(_cfg(tmp_path))

    def test_runx1_strong_residue(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_DEL_1AA, _ann("RUNX1", ConsequenceType.INFRAME_DELETION, 107))
        assert r.triggered and r.strength == CriterionStrength.STRONG

    def test_runx1_supporting_region(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_DEL_1AA, _ann("RUNX1", ConsequenceType.INFRAME_DELETION, 150))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_runx1_outside_not_met(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_DEL_1AA, _ann("RUNX1", ConsequenceType.INFRAME_DELETION, 300))
        assert not r.triggered

    def test_myoc_in_domain_moderate(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_DEL_2AA, _ann("MYOC", ConsequenceType.INFRAME_DELETION, 300))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_myoc_outside_not_met(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_DEL_2AA, _ann("MYOC", ConsequenceType.INFRAME_DELETION, 100))
        assert not r.triggered

    def test_dicer1_outside_supporting(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_DEL_2AA, _ann("DICER1", ConsequenceType.INFRAME_INSERTION, 100))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_dicer1_repeat_denied(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_DEL_2AA, _ann("DICER1", ConsequenceType.INFRAME_DELETION, 607))
        assert not r.triggered

    def test_mecp2_small_indel_supporting_even_in_deny(self, tmp_path):
        # 1-aa indel in the Pro-rich deny region → Supporting (size tier wins).
        r = self._ev(tmp_path).evaluate(_DEL_1AA, _ann("MECP2", ConsequenceType.INFRAME_DELETION, 390))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_mecp2_large_indel_in_deny_not_met(self, tmp_path):
        # 3-aa indel in the deny region → withheld.
        r = self._ev(tmp_path).evaluate(_DEL_3AA, _ann("MECP2", ConsequenceType.INFRAME_DELETION, 390))
        assert not r.triggered

    def test_mecp2_large_indel_outside_moderate(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_DEL_3AA, _ann("MECP2", ConsequenceType.INFRAME_DELETION, 100))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_cdh1_indel_not_met_but_stoploss_moderate(self, tmp_path):
        ev = self._ev(tmp_path)
        indel = ev.evaluate(_DEL_2AA, _ann("CDH1", ConsequenceType.INFRAME_DELETION, 200))
        assert not indel.triggered
        sl = ev.evaluate(_var("A", "T"), _ann("CDH1", ConsequenceType.STOP_LOST))
        assert sl.triggered and sl.strength == CriterionStrength.MODERATE

    def test_cyp1b1_stoploss_na(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_var("A", "T"), _ann("CYP1B1", ConsequenceType.STOP_LOST))
        assert not r.triggered

    def test_cyp1b1_indel_default_moderate(self, tmp_path):
        # CYP1B1 has only a stoploss rule → in-frame indel uses the flat default.
        r = self._ev(tmp_path).evaluate(_DEL_2AA, _ann("CYP1B1", ConsequenceType.INFRAME_DELETION, 100))
        assert r.triggered and r.strength == CriterionStrength.MODERATE


def test_committed_pm4_regions_present():
    import csv
    tsv = Path(__file__).resolve().parents[2] / "resources" / "shared" / "pm4_regions.tsv"
    with tsv.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    genes = {r["gene_symbol"] for r in rows}
    assert {"RUNX1", "MYOC", "DICER1", "MECP2", "CDKL5", "FOXG1", "CDH1", "ATM", "CYP1B1"} <= genes
    runx1_strong = next(r for r in rows if r["gene_symbol"] == "RUNX1" and r["strength"] == "strong")
    assert "107" in runx1_strong["residues"].split(",")
