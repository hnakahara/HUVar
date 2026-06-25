"""REVEL PP3/BP4 thresholds: ClinGen/Pejaver defaults + per-gene VCEP overrides.

Covers the RevelSpec loader, the default Pejaver tier mapping, the per-gene
override (with the Supporting-only cap), and the build-script extraction of
REVEL cutoffs from cspec PP3/BP4 descriptions.
"""
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
from acmg_classifier.criteria.revel_genes import RevelSpec
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, RevelData
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, CriterionStrength, InSilicoTool,
)
from acmg_classifier.models.variant import VariantRecord

_HEADER = (
    "gene_symbol\trevel_pp3_supporting\trevel_pp3_moderate\trevel_pp3_strong\t"
    "revel_bp4_supporting\trevel_bp4_moderate\trevel_bp4_strong\n"
)
_ROWS = (
    "MYH7\t0.7\t\t\t0.4\t\t\n"                       # single Supporting (HCM)
    "MYOC\t0.644\t0.773\t0.932\t0.29\t0.183\t0.016\n"  # full tiered
    "VHL\t0.664\t\t\t\t\t\n"                          # PP3-only, no BP4
)


def _spec_tsv(tmp_path: Path) -> Path:
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_HEADER + _ROWS, encoding="utf-8")
    return p


def _cfg(tmp_path: Path) -> MagicMock:
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _spec_tsv(tmp_path)
    cfg.insilico_tool = InSilicoTool.REVEL
    # Opt-in auxiliary predictors off by default (a bare MagicMock attribute would
    # otherwise be truthy and wrongly activate the BayesDel/CADD per-gene rules).
    cfg.use_bayesdel = False
    cfg.use_cadd = False
    return cfg


def _ann(gene: str, score: float | None) -> AnnotationData:
    return AnnotationData(
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        )],
        revel=RevelData(score=score) if score is not None else None,
    )


def _snv() -> VariantRecord:
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


class TestRevelSpec:
    def test_single_supporting_gene(self, tmp_path):
        rule = RevelSpec(_spec_tsv(tmp_path)).get("MYH7")
        assert rule.pp3 == {CriterionStrength.SUPPORTING: 0.7}
        assert rule.bp4 == {CriterionStrength.SUPPORTING: 0.4}

    def test_tiered_gene(self, tmp_path):
        rule = RevelSpec(_spec_tsv(tmp_path)).get("MYOC")
        assert rule.pp3[CriterionStrength.STRONG] == 0.932
        assert rule.bp4[CriterionStrength.STRONG] == 0.016

    def test_pp3_only_gene_has_no_bp4(self, tmp_path):
        rule = RevelSpec(_spec_tsv(tmp_path)).get("VHL")
        assert rule.pp3 == {CriterionStrength.SUPPORTING: 0.664}
        assert rule.bp4 == {}

    def test_unknown_gene_none(self, tmp_path):
        assert RevelSpec(_spec_tsv(tmp_path)).get("TP53") is None


class TestPP3RevelDefaults:
    def test_default_tiers(self, tmp_path):
        ev = PP3Evaluator(_cfg(tmp_path))
        # TP53 has no per-gene rule → genome-wide Bergquist 2024 tiers apply.
        assert ev.evaluate(_snv(), _ann("TP53", 0.95)).strength == CriterionStrength.STRONG
        assert ev.evaluate(_snv(), _ann("TP53", 0.90)).strength == CriterionStrength.THREE_POINT
        assert ev.evaluate(_snv(), _ann("TP53", 0.80)).strength == CriterionStrength.MODERATE
        assert ev.evaluate(_snv(), _ann("TP53", 0.65)).strength == CriterionStrength.SUPPORTING
        assert not ev.evaluate(_snv(), _ann("TP53", 0.50)).triggered

    def test_no_score_not_met(self, tmp_path):
        ev = PP3Evaluator(_cfg(tmp_path))
        assert not ev.evaluate(_snv(), _ann("TP53", None)).triggered


class TestPP3RevelPerGene:
    def test_supporting_only_gene_is_capped(self, tmp_path):
        ev = PP3Evaluator(_cfg(tmp_path))
        # MYH7 caps PP3 at Supporting: a very high REVEL must NOT reach Strong.
        r = ev.evaluate(_snv(), _ann("MYH7", 0.99))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING
        assert "MYH7 VCEP cutoff" in r.evidence

    def test_supporting_only_below_cutoff_not_met(self, tmp_path):
        ev = PP3Evaluator(_cfg(tmp_path))
        # MYH7 cutoff 0.7; 0.65 (above the genome-wide 0.644) must NOT fire.
        assert not ev.evaluate(_snv(), _ann("MYH7", 0.65)).triggered

    def test_tiered_gene_reaches_strong(self, tmp_path):
        ev = PP3Evaluator(_cfg(tmp_path))
        assert ev.evaluate(_snv(), _ann("MYOC", 0.95)).strength == CriterionStrength.STRONG


class TestBP4RevelPerGene:
    def test_default_tiers(self, tmp_path):
        ev = BP4Evaluator(_cfg(tmp_path))
        # Bergquist 2024 REVEL BP4: Strong <=0.016, ThreePoint <=0.052,
        # Moderate <=0.183, Supporting <=0.290 (no Very Strong tier).
        assert ev.evaluate(_snv(), _ann("TP53", 0.010)).strength == CriterionStrength.STRONG
        assert ev.evaluate(_snv(), _ann("TP53", 0.040)).strength == CriterionStrength.THREE_POINT
        assert ev.evaluate(_snv(), _ann("TP53", 0.100)).strength == CriterionStrength.MODERATE
        assert ev.evaluate(_snv(), _ann("TP53", 0.250)).strength == CriterionStrength.SUPPORTING
        assert not ev.evaluate(_snv(), _ann("TP53", 0.400)).triggered

    def test_supporting_only_gene_is_capped(self, tmp_path):
        ev = BP4Evaluator(_cfg(tmp_path))
        # MYH7 BP4 cutoff 0.4 at Supporting only — a near-zero score stays Supporting.
        r = ev.evaluate(_snv(), _ann("MYH7", 0.001))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_pp3_only_gene_has_no_bp4(self, tmp_path):
        ev = BP4Evaluator(_cfg(tmp_path))
        # VHL carries no BP4 cutoff → falls back to the genome-wide default,
        # under which a low REVEL still yields a benign call.
        assert ev.evaluate(_snv(), _ann("VHL", 0.001)).triggered


# ---- build_disease_thresholds.py REVEL extraction ----

_BUILD = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_bspec = importlib.util.spec_from_file_location("build_disease_thresholds", _BUILD)
b = importlib.util.module_from_spec(_bspec)
_bspec.loader.exec_module(b)


def _rs(label: str, *strengths: tuple[str, str]) -> dict:
    return {"criteriaCodes": [{
        "label": label,
        "evidenceStrengths": [
            {"label": s, "applicability": "Applicable", "description": d}
            for s, d in strengths
        ],
    }]}


class TestRevelExtraction:
    def test_operator_value(self):
        rs = _rs("PP3", ("Supporting", "Use REVEL score >= 0.7 for PP3."))
        assert b._revel_tiers(rs, "PP3") == {"supporting": 0.7}

    def test_far_threshold_after_citation(self):
        # The cutoff sits well past the first REVEL mention (HCM Ioannidis text).
        desc = ("Meta-predictors, such as REVEL, are preferred. Use of REVEL "
                "(Ioannidis et al. 2016[<sup>14</sup>](#pmid_27666373)) is "
                "recommended at thresholds of ≥0.70 for PP3.")
        assert b._revel_tiers(_rs("PP3", ("Supporting", desc)), "PP3") == {"supporting": 0.7}

    def test_operator_free_cutoff(self):
        rs = _rs("BP4", ("Supporting", "REVEL. Use 0.326 as a discriminatory cut-off value."))
        assert b._revel_tiers(rs, "BP4") == {"supporting": 0.326}

    def test_tiered_ranges_take_firing_edge(self):
        rs = _rs(
            "PP3",
            ("Supporting", "REVEL score 0.644 - 0.772."),
            ("Moderate", "REVEL score 0.773 - 0.931."),
            ("Strong", "REVEL score >= 0.932."),
        )
        assert b._revel_tiers(rs, "PP3") == {"supporting": 0.644, "moderate": 0.773, "strong": 0.932}

    def test_bp4_range_takes_upper_edge(self):
        rs = _rs("BP4", ("Supporting", "REVEL between 0.184 and 0.290."))
        assert b._revel_tiers(rs, "BP4") == {"supporting": 0.29}

    def test_other_tool_number_not_captured(self):
        # SpliceAI's number must not be read as the REVEL cutoff.
        rs = _rs("PP3", ("Supporting", "REVEL >= 0.7; otherwise SpliceAI >= 0.2 supports PP3."))
        assert b._revel_tiers(rs, "PP3") == {"supporting": 0.7}

    def test_monotonicity_guard_drops_typo(self):
        # GN208-style typo: Moderate cutoff below Supporting is dropped.
        rs = _rs(
            "PP3",
            ("Supporting", "REVEL >= 0.644."),
            ("Moderate", "REVEL score 0.774 - 0.092."),
        )
        assert b._revel_tiers(rs, "PP3") == {"supporting": 0.644}

    def test_no_number_is_empty(self):
        # SCN-style prose with no numeric cutoff → no extraction (default applies).
        rs = _rs("PP3", ("Moderate", "Follow ClinGen recommendations using REVEL."))
        assert b._revel_tiers(rs, "PP3") == {}
