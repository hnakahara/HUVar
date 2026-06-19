"""BS2 count→strength tiers (``bs2_strength`` column).

Several VCEPs tier BS2 on the observed gnomAD count (e.g. GUCY2D Strong>=6 /
Supporting>=3; BMPR2 Strong>=3 / Moderate>=2 / Supporting>=1). The
``bs2_strength`` column encodes these as ``Strength:count`` pairs; the evaluator
fires at the strongest tier whose mode-specific count threshold is met, instead
of the flat always-Strong default. These tests cover the loader, the evaluator,
and the build-script extraction.
"""
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
from acmg_classifier.criteria.bs2_genes import BS2Applicability
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)

# GUCY2D (AR): Strong>=6 / Supporting>=3 homozygotes.
# IL2RG  (XL): Strong>=3 / Supporting>=2 hemizygotes.
# BMPR2  (AD, hom_only): Strong>=3 / Moderate>=2 / Supporting>=1 homozygotes.
_TSV = (
    "gene_symbol\tbs2\tinheritance\tbs2_count\tbs2_strength\tbs2_female_only\tbs2_hom_only\n"
    "GUCY2D\tapplicable\tAR\t3\tStrong:6,Supporting:3\t\t\n"
    "IL2RG\tapplicable\tXL\t3\tStrong:3,Supporting:2\t\t\n"
    "BMPR2\tapplicable\tAD\t1\tStrong:3,Moderate:2,Supporting:1\t\t1\n"
    "RAG1\tapplicable\tAR\t1\tStrong:3,Supporting:1\t\t\n"
    "GENE0\tapplicable\tAR\t2\t\t\t\n"   # no tiers → flat Strong at bs2_count
)


def _tsv(tmp_path):
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_TSV, encoding="utf-8")
    return p


def _cfg(tmp_path):
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _tsv(tmp_path)
    cfg.bs2_min_homalt = 2
    cfg.bs2_min_hemi = 2
    cfg.bs2_min_het = 3
    return cfg


def _ann(gene, **gd):
    return AnnotationData(
        gnomad=GnomADData(filter_pass=True, **gd),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        )],
    )


def _snv():
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


class TestTierLoader:
    def test_parses_two_tiers_sorted_strongest_first(self, tmp_path):
        s = BS2Applicability(_tsv(tmp_path))
        assert s.tiers("GUCY2D") == (
            (CriterionStrength.STRONG, 6), (CriterionStrength.SUPPORTING, 3),
        )

    def test_parses_three_tiers(self, tmp_path):
        s = BS2Applicability(_tsv(tmp_path))
        assert s.tiers("BMPR2") == (
            (CriterionStrength.STRONG, 3),
            (CriterionStrength.MODERATE, 2),
            (CriterionStrength.SUPPORTING, 1),
        )

    def test_no_tiers_is_empty(self, tmp_path):
        s = BS2Applicability(_tsv(tmp_path))
        assert s.tiers("GENE0") == ()
        assert s.tiers("UNSEEN") == ()
        assert s.tiers(None) == ()


class TestTieredEvaluator:
    def test_gucy2d_high_count_strong(self, tmp_path):
        # 6 homozygotes -> Strong tier.
        r = BS2Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("GUCY2D", ac=6, nhomalt=6))
        assert r.triggered and r.strength == CriterionStrength.STRONG

    def test_gucy2d_mid_count_supporting(self, tmp_path):
        # 3-5 homozygotes -> only the Supporting tier is met.
        r = BS2Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("GUCY2D", ac=4, nhomalt=4))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_gucy2d_below_lowest_tier_not_met(self, tmp_path):
        r = BS2Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("GUCY2D", ac=2, nhomalt=2))
        assert not r.triggered

    def test_il2rg_hemizygote_tiers(self, tmp_path):
        strong = BS2Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("IL2RG", ac=3, nhemi=3))
        assert strong.triggered and strong.strength == CriterionStrength.STRONG
        sup = BS2Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("IL2RG", ac=2, nhemi=2))
        assert sup.triggered and sup.strength == CriterionStrength.SUPPORTING

    def test_bmpr2_homozygote_moderate(self, tmp_path):
        # AD incomplete-penetrance gene counts homozygotes; 2 hom -> Moderate.
        r = BS2Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("BMPR2", ac=50, nhomalt=2))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_non_tiered_gene_fires_strong_by_default(self, tmp_path):
        # GENE0 has no tiers; flat bs2_count=2 -> Strong (unchanged behaviour).
        r = BS2Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("GENE0", ac=2, nhomalt=2))
        assert r.triggered and r.strength == CriterionStrength.STRONG


class TestBs2StrengthExtraction:
    def _bs2(self, *strengths):
        return {"criteriaCodes": [{"label": "BS2", "evidenceStrengths": [
            {"label": lbl, "applicability": "Applicable", "description": desc}
            for lbl, desc in strengths
        ]}]}

    def test_two_tier_homozygotes(self):
        rs = self._bs2(
            ("Strong", "BS2_Strong: applied if observed in at least 6 homozygotes in gnomAD."),
            ("Supporting", "BS2_Supporting: applied if observed in at least 3 homozygotes."),
        )
        assert bdt._bs2_strength(rs) == "Strong:6,Supporting:3"

    def test_three_tier(self):
        rs = self._bs2(
            ("Strong", "observed in >=3 homozygotes"),
            ("Moderate", "observed in >=2 homozygotes"),
            ("Supporting", "observed in at least 1 homozygote"),
        )
        assert bdt._bs2_strength(rs) == "Strong:3,Moderate:2,Supporting:1"

    def test_single_tier_is_blank(self):
        # One strength with a count is captured by bs2_count, not bs2_strength.
        rs = self._bs2(("Strong", "observed in at least 3 homozygotes"))
        assert bdt._bs2_strength(rs) == ""

    def test_gnomad_anchored_count_wins_over_literature(self):
        # GUCY2D Strong cites BOTH a phenotyped-literature ">=3 homozygotes" and a
        # gnomAD ">=6 homozygotes"; the app uses gnomAD, so Strong must be 6 not 3.
        rs = self._bs2(
            ("Strong",
             "Variant is present in >= 3 homozygotes without any features of the "
             "phenotype ... unaffected by age 40. Alternatively, this strength can "
             "be applied if the variant is present in >= 6 homozygotes in gnomAD "
             "v.4.1.0 or later."),
            ("Supporting", "Variant is present in >= 3 homozygotes in gnomAD v.4.1.0 or later."),
        )
        assert bdt._bs2_strength(rs) == "Strong:6,Supporting:3"
