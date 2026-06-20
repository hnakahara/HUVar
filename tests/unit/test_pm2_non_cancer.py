"""PM2 non-cancer subset (ENIGMA BRCA1/2): pm2_subset=non_cancer.

The BRCA1/2 VCEP judges PM2 absence on gnomAD's non-cancer subset, so a variant
present only in cancer cohorts still earns PM2. The evaluator uses the non-cancer
AF (gd.af_non_cancer) when available and falls back to the overall AF when the
gnomAD DB predates the column.
"""
import csv
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
from acmg_classifier.criteria.pm2_genes import PM2Spec
from acmg_classifier.local_db.gnomad_db import _merge_rows
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_ROOT = Path(__file__).resolve().parents[2]
_BDT = _ROOT / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)

_HEADER = ("gene_symbol\tpm2_threshold\tpm2_strength\tpm2_basis\tpm2_subpop\t"
           "pm2_zygosity\tpm2_subset\n")
_ROWS = (
    "BRCA1\t0\t\t\t\t\tnon_cancer\n"   # absent from the non-cancer subset
    "KRAS\t0\t\t\t\t\t\n"              # absent (overall), no subset
)


def _cfg(tmp_path):
    p = tmp_path / "dp.tsv"
    p.write_text(_HEADER + _ROWS, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = p
    cfg.gene_inheritance_tsv = tmp_path / "missing.tsv"
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
    return VariantRecord(chrom="chr17", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


class TestLoader:
    def test_subset_parsed(self, tmp_path):
        s = PM2Spec(_cfg(tmp_path).disease_prevalence_tsv)
        assert s.get("BRCA1").subset == "non_cancer"
        assert s.get("KRAS").subset == ""


class TestEvaluator:
    def test_present_in_cancer_only_still_pm2(self, tmp_path):
        # Present overall (cancer cohorts: af=2e-4, ac=5) but ABSENT in the
        # non-cancer subset → PM2 should still fire for BRCA1.
        r = PM2Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BRCA1", af=2e-4, popmax_af=2e-4, ac=5, af_non_cancer=0.0))
        assert r.triggered and "non-cancer" in r.evidence

    def test_present_in_non_cancer_blocks_pm2(self, tmp_path):
        # Present in the non-cancer subset → not absent → PM2 not met.
        r = PM2Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BRCA1", af=2e-4, popmax_af=2e-4, ac=5, af_non_cancer=2e-4))
        assert not r.triggered

    def test_fallback_to_overall_when_subset_absent(self, tmp_path):
        # gnomAD DB predates the column (af_non_cancer=None) → fall back to overall
        # AF; present overall → not absent → not met, with a fallback note.
        r = PM2Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BRCA1", af=2e-4, popmax_af=2e-4, ac=5, af_non_cancer=None))
        assert not r.triggered and "non-cancer subset unavailable" in r.evidence

    def test_truly_absent_fires(self, tmp_path):
        r = PM2Evaluator(_cfg(tmp_path)).evaluate(
            _snv(), _ann("BRCA1", af=0.0, ac=0, af_non_cancer=0.0))
        assert r.triggered


class TestGnomadMerge:
    def test_merge_reads_non_cancer_at_index14(self):
        # Row layout: 0..10 stats, 11 filters, 12/13 grpmax, 14 af_non_cancer.
        row = (2e-4, 1000, 5, 0, 0, 2e-4, "nfe", 1e-4, None, None, None,
               None, None, None, 1e-5)
        gd = _merge_rows([row])
        assert gd.af_non_cancer == 1e-5

    def test_merge_short_tuple_degrades_to_none(self):
        # An old-schema 14-field row (no index 14) → af_non_cancer None.
        row = (2e-4, 1000, 5, 0, 0, 2e-4, "nfe", 1e-4, None, None, None,
               None, None, None)
        gd = _merge_rows([row])
        assert gd.af_non_cancer is None


class TestExtraction:
    def _pm2(self, desc):
        return {"criteriaCodes": [{"label": "PM2", "evidenceStrengths": [
            {"label": "Supporting", "applicability": "Applicable", "description": desc},
        ]}]}

    def test_non_cancer_detected(self):
        rs = self._pm2("Absent from controls in gnomAD v2.1 (non-cancer, exome only subset).")
        assert bdt._pm2_subset(rs) == "non_cancer"

    def test_plain_blank(self):
        rs = self._pm2("Absent or extremely rare in gnomAD.")
        assert bdt._pm2_subset(rs) == ""


def test_committed_brca_non_cancer():
    tsv = _ROOT / "resources" / "shared" / "disease_prevalence.tsv"
    with tsv.open(encoding="utf-8") as f:
        rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
    assert rows["BRCA1"]["pm2_subset"] == "non_cancer"
    assert rows["BRCA2"]["pm2_subset"] == "non_cancer"
