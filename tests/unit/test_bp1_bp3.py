"""BP1 (VCEP-gated, gene-specific) and BP3 (VCEP applicability gate)."""
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.bp_genes import BPApplicability
from acmg_classifier.criteria.benign.bp1 import BP1Evaluator
from acmg_classifier.criteria.benign.bp3 import BP3Evaluator
from acmg_classifier.models.annotation import (
    AnnotationData, GnomADData, ConsequenceInfo, RepeatMaskerRegion, SpliceScore,
)
from acmg_classifier.models.enums import CriterionStrength
from acmg_classifier.models.enums import Assembly, ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_BDT = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


def _bp1_code(applicability, desc):
    return {"label": "BP1", "evidenceStrengths": [
        {"label": "Supporting", "applicability": applicability, "description": desc},
    ]}


class TestBp1Extraction:
    def test_missense_target(self):
        rs = {"criteriaCodes": [_bp1_code("Applicable", "Apply to all missense variants.")]}
        assert bdt._bp1_applicability(rs) == ("applicable", "missense")

    def test_truncating_gof_target(self):
        rs = {"criteriaCodes": [_bp1_code(
            "Applicable",
            "Disease mechanism is gain-of-function; BP1 should be used for any "
            "truncating variant (nonsense, frameshift).",
        )]}
        assert bdt._bp1_applicability(rs) == ("applicable", "truncating")

    def test_declined(self):
        rs = {"criteriaCodes": [{"label": "BP1", "evidenceStrengths": [
            {"label": "Supporting", "applicability": "Not Applicable for this VCEP"},
        ]}]}
        assert bdt._bp1_applicability(rs) == ("not_applicable", "")

    def test_no_code(self):
        assert bdt._bp1_applicability({"criteriaCodes": [{"label": "PM1"}]}) == ("", "")


class TestBp3Extraction:
    def test_applicable(self):
        rs = {"criteriaCodes": [{"label": "BP3", "evidenceStrengths": [
            {"label": "Supporting", "applicability": "Applicable",
             "description": "In-frame indel in a repetitive region."},
        ]}]}
        assert bdt._bp3_applicability(rs) == "applicable"

    def test_declined(self):
        rs = {"criteriaCodes": [{"label": "BP3", "evidenceStrengths": [
            {"label": "Supporting", "applicability": "Not Applicable for this VCEP"},
        ]}]}
        assert bdt._bp3_applicability(rs) == "not_applicable"


_COLS = "gene_symbol\tbp1\tbp1_target\tbp1_exclude\tbp1_strength\tbp1_no_splice\tbp3\tbp3_regions\n"


def _tsv(tmp_path):
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(
        _COLS
        + "PALB2\tapplicable\tmissense\t\t\t\tnot_applicable\t\n"
        + "NRAS\tapplicable\ttruncating\t\t\t\tnot_applicable\t\n"
        + "APC\tapplicable\tmissense\t1021-1035\t\t\tnot_applicable\t\n"
        + "BRCA1\tapplicable\tbroad\t2-101;1391-1424;1650-1857\tStrong\tyes\tnot_applicable\t\n"
        + "MYH7\tnot_applicable\t\t\t\t\tnot_applicable\t\n"
        + "BMPR2\tnot_applicable\t\t\t\t\tapplicable\t\n"
        + "RPGR\tnot_applicable\t\t\t\t\tapplicable\t585-1078\n",
        encoding="utf-8",
    )
    return p


class TestBPApplicabilityLoader:
    def test_loads(self, tmp_path):
        a = BPApplicability(_tsv(tmp_path))
        assert a.bp1("PALB2") == "applicable" and a.bp1_target("PALB2") == "missense"
        assert a.bp1("NRAS") == "applicable" and a.bp1_target("NRAS") == "truncating"
        assert a.bp1("MYH7") == "not_applicable"
        assert a.bp1("UNSEEN") == "" and a.bp1_target("UNSEEN") == ""
        assert a.bp3("BMPR2") == "applicable" and a.bp3("MYH7") == "not_applicable"

    def test_missing_file(self, tmp_path):
        a = BPApplicability(tmp_path / "nope.tsv")
        assert a.bp1("PALB2") == "" and a.bp3("BMPR2") == ""


def _cfg(tmp_path):
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _tsv(tmp_path)
    return cfg


def _ann(gene, ctype, repeat=False, pos=200, splice=None):
    rep = RepeatMaskerRegion(in_repeat=repeat, repeat_class="SINE", repeat_name="Alu")
    kw = {}
    if splice is not None:
        kw["splice"] = SpliceScore(tool="spliceai", max_delta=splice)
    return AnnotationData(
        gnomad=GnomADData(),
        consequences=[ConsequenceInfo(
            transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene, consequence=ctype,
            biotype="protein_coding", is_mane_select=True, protein_position=pos,
        )],
        repeat=rep,
        **kw,
    )


def _snv():
    return VariantRecord(chrom="chr1", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38)


class TestBP1Evaluator:
    def test_missense_gene_fires_on_missense(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("PALB2", ConsequenceType.MISSENSE))
        assert r.triggered and r.criterion == ACMGCriterion.BP1

    def test_missense_gene_not_on_truncating(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("PALB2", ConsequenceType.STOP_GAINED))
        assert not r.triggered

    def test_truncating_gene_fires_on_truncating(self, tmp_path):
        ev = BP1Evaluator(_cfg(tmp_path))
        assert ev.evaluate(_snv(), _ann("NRAS", ConsequenceType.STOP_GAINED)).triggered
        assert ev.evaluate(_snv(), _ann("NRAS", ConsequenceType.FRAMESHIFT)).triggered

    def test_truncating_gene_not_on_missense(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("NRAS", ConsequenceType.MISSENSE))
        assert not r.triggered

    def test_not_applicable_gene(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("MYH7", ConsequenceType.MISSENSE))
        assert not r.triggered

    def test_no_vcep_gene_not_applied(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann("NOVCEP", ConsequenceType.MISSENSE))
        assert not r.triggered


class TestBP1RegionAndStrength:
    def test_apc_excluded_region_withheld(self, tmp_path):
        # APC missense inside the β-catenin repeat (1021-1035) -> no BP1.
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("APC", ConsequenceType.MISSENSE, pos=1028)
        )
        assert not r.triggered and "excluded region" in r.evidence

    def test_apc_outside_region_fires(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("APC", ConsequenceType.MISSENSE, pos=500)
        )
        assert r.triggered

    def test_brca1_broad_strong_no_splice(self, tmp_path):
        # Silent variant outside domains, SpliceAI 0.0 -> BP1_Strong.
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BRCA1", ConsequenceType.SYNONYMOUS, pos=800, splice=0.0)
        )
        assert r.triggered and r.strength == CriterionStrength.STRONG

    def test_brca1_in_domain_withheld(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BRCA1", ConsequenceType.MISSENSE, pos=50, splice=0.0)  # RING 2-101
        )
        assert not r.triggered

    def test_brca1_splice_impact_withheld(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BRCA1", ConsequenceType.MISSENSE, pos=800, splice=0.5)
        )
        assert not r.triggered and "splice impact" in r.evidence

    def test_brca1_no_splice_data_withheld(self, tmp_path):
        r = BP1Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BRCA1", ConsequenceType.MISSENSE, pos=800, splice=None)
        )
        assert not r.triggered and "unavailable" in r.evidence


class TestBP3Gate:
    def test_not_applicable_gene_suppressed(self, tmp_path):
        # MYH7 BP3 not_applicable: even an in-frame indel in a repeat is withheld.
        r = BP3Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("MYH7", ConsequenceType.INFRAME_DELETION, repeat=True)
        )
        assert not r.triggered and "not applicable" in r.evidence

    def test_applicable_gene_uses_heuristic(self, tmp_path):
        r = BP3Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BMPR2", ConsequenceType.INFRAME_DELETION, repeat=True)
        )
        assert r.triggered

    def test_no_vcep_gene_uses_heuristic(self, tmp_path):
        r = BP3Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("NOVCEP", ConsequenceType.INFRAME_DELETION, repeat=True)
        )
        assert r.triggered

    def test_applicable_but_not_in_repeat(self, tmp_path):
        r = BP3Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BMPR2", ConsequenceType.INFRAME_DELETION, repeat=False)
        )
        assert not r.triggered and "repeat" in r.evidence

    def test_region_restricted_in_region_fires_without_dfam(self, tmp_path):
        # RPGR BP3 restricted to ORF15 585-1078: in-region fires even without Dfam.
        r = BP3Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("RPGR", ConsequenceType.INFRAME_DELETION, repeat=False, pos=700)
        )
        assert r.triggered and "BP3 repetitive region" in r.evidence

    def test_region_restricted_outside_region_withheld(self, tmp_path):
        r = BP3Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("RPGR", ConsequenceType.INFRAME_DELETION, repeat=True, pos=200)
        )
        assert not r.triggered and "outside the VCEP BP3 region" in r.evidence
