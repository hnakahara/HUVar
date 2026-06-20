"""Per-gene SpliceAI no-impact cutoff for BP4 / BP7 (bp4/bp7_splice_cutoff).

Most VCEPs treat SpliceAI <= 0.10 as "no impact"; some override it per gene —
the LGMD panels tighten to 0.05, RUNX1/RPGR loosen to 0.20, GAA/TP53 to 0.2,
HBA2/HBB to 0.3. These tests cover the loader, the BP4 splice branch, the BP7
splice gate, the build-script extraction, and the committed TSV.
"""
import csv
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
from acmg_classifier.criteria.benign.bp7 import BP7Evaluator
from acmg_classifier.criteria.bp_genes import BPApplicability
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, SpliceScore,
)
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, CriterionStrength, InSilicoTool,
)
from acmg_classifier.models.variant import VariantRecord

_ROOT = Path(__file__).resolve().parents[2]
_BDT = _ROOT / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)

_TSV = (
    "gene_symbol\tbp7_phylop\tbp4_splice_cutoff\tbp7_splice_cutoff\n"
    "DYSF\t\t0.05\t0.05\n"    # LGMD: tighter
    "RUNX1\t\t0.20\t0.20\n"   # looser
    "GAA\t\t0.2\t\n"          # BP4-only override
    "NOSPEC\t\t\t\n"          # default 0.10
)


def _tsv(tmp_path):
    p = tmp_path / "dp.tsv"
    p.write_text(_TSV, encoding="utf-8")
    return p


class TestLoader:
    def test_reads_cutoffs(self, tmp_path):
        s = BPApplicability(_tsv(tmp_path))
        assert s.bp4_splice_cutoff("DYSF") == 0.05
        assert s.bp7_splice_cutoff("DYSF") == 0.05
        assert s.bp4_splice_cutoff("GAA") == 0.2
        assert s.bp7_splice_cutoff("GAA") is None   # BP4-only
        assert s.bp4_splice_cutoff("NOSPEC") is None
        assert s.bp7_splice_cutoff("UNSEEN") is None


def _cfg(tmp_path):
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _tsv(tmp_path)
    cfg.insilico_tool = InSilicoTool.REVEL
    cfg.bp7_phylop_max = 2.0
    cfg.phylop_bigwig = None
    return cfg


def _snv():
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


def _ann_splice(gene, max_delta, consequence=ConsequenceType.INTRON, dist=20):
    return AnnotationData(
        splice=SpliceScore(tool="spliceai", is_available=True, max_delta=max_delta),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=consequence, biotype="protein_coding",
            intron_distance_from_splice=dist,
        )],
    )


class TestBP4SpliceCutoff:
    def test_tighter_cutoff_blocks_between(self, tmp_path):
        # SpliceAI 0.08: benign under the 0.10 default but NOT under DYSF's 0.05.
        r = BP4Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann_splice("DYSF", 0.08, ConsequenceType.SYNONYMOUS))
        assert not r.triggered

    def test_tighter_cutoff_fires_below(self, tmp_path):
        r = BP4Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann_splice("DYSF", 0.04, ConsequenceType.SYNONYMOUS))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_looser_cutoff_fires_above_default(self, tmp_path):
        # SpliceAI 0.15: not benign under 0.10 default, IS under RUNX1's 0.20.
        r = BP4Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann_splice("RUNX1", 0.15, ConsequenceType.SYNONYMOUS))
        assert r.triggered

    def test_default_gene_uses_010(self, tmp_path):
        r = BP4Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann_splice("NOSPEC", 0.15, ConsequenceType.SYNONYMOUS))
        assert not r.triggered
        r2 = BP4Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann_splice("NOSPEC", 0.08, ConsequenceType.SYNONYMOUS))
        assert r2.triggered


class TestBP7SpliceGate:
    def test_tighter_cutoff_blocks_intronic(self, tmp_path):
        # Deep intronic (dist 20), SpliceAI 0.08: BP7 fires under default but
        # NOT under DYSF's 0.05 splice gate.
        r = BP7Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann_splice("DYSF", 0.08))
        assert not r.triggered

    def test_bp4_only_gene_keeps_bp7_default(self, tmp_path):
        # GAA overrides BP4 (0.2) but NOT BP7 → BP7 still uses 0.10. SpliceAI 0.15
        # must NOT fire BP7.
        r = BP7Evaluator(_cfg(tmp_path)).evaluate(_snv(), _ann_splice("GAA", 0.15))
        assert not r.triggered


class TestExtraction:
    def _code(self, label, desc):
        return {"criteriaCodes": [{"label": label, "evidenceStrengths": [
            {"label": "Supporting", "applicability": "Applicable", "description": desc},
        ]}]}

    def test_lgmd_005(self):
        rs = self._code("BP4", "splice prediction algorithms (SpliceAI <= 0.05) predict no impact")
        assert bdt._splice_cutoff(rs, "BP4") == "0.05"

    def test_runx1_020(self):
        rs = self._code("BP7", "SpliceAI <= 0.20 and not highly conserved")
        assert bdt._splice_cutoff(rs, "BP7") == "0.20"

    def test_default_010_not_emitted(self):
        rs = self._code("BP4", "SpliceAI <= 0.10 (Walker default)")
        assert bdt._splice_cutoff(rs, "BP4") == ""

    def test_pathogenic_side_ignored(self):
        # "predicts impact >= 0.2" is the pathogenic side, not the no-impact cutoff.
        rs = self._code("BP4", "exclude variants where SpliceAI >= 0.2 (predicted impact)")
        assert bdt._splice_cutoff(rs, "BP4") == ""


def test_committed_tsv_values():
    tsv = _ROOT / "resources" / "shared" / "disease_prevalence.tsv"
    with tsv.open(encoding="utf-8") as f:
        rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
    assert rows["DYSF"]["bp4_splice_cutoff"] == "0.05"
    assert rows["DYSF"]["bp7_splice_cutoff"] == "0.05"
    assert rows["RPGR"]["bp7_splice_cutoff"] == "0.2"
    assert rows["RPGR"]["bp4_splice_cutoff"] == ""
    assert rows["TP53"]["bp4_splice_cutoff"] == "0.2"
