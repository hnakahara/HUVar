"""PM2 read-depth gate (ENIGMA BRCA1/2: "region read depth >= 25").

A variant "absent" from gnomAD in a poorly-covered region does not earn PM2 —
absence there is not callable. The gate uses the gnomAD coverage DuckDB
(CoverageDB); it is skipped when coverage is unavailable / unknown.
"""
import csv
import gzip
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
from acmg_classifier.criteria.pm2_genes import PM2Spec
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_ROOT = Path(__file__).resolve().parents[2]
_BDT = _ROOT / "scripts" / "build_disease_thresholds.py"
_s = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT)
bdt = importlib.util.module_from_spec(_s)
_s.loader.exec_module(bdt)

_HEADER = ("gene_symbol\tpm2_threshold\tpm2_strength\tpm2_basis\tpm2_subpop\t"
           "pm2_zygosity\tpm2_subset\tpm2_min_depth\n")
_ROWS = "BRCA1\t0\t\t\t\t\t\t25\n"


def _cfg(tmp_path):
    p = tmp_path / "dp.tsv"
    p.write_text(_HEADER + _ROWS, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = p
    cfg.gene_inheritance_tsv = tmp_path / "missing.tsv"
    cfg.gnomad_coverage_db = tmp_path / "absent.duckdb"
    return cfg


class _Cov:
    def __init__(self, available, dp):
        self._a, self._dp = available, dp

    @property
    def available(self):
        return self._a

    def mean_depth(self, chrom, start, end=None):
        return self._dp


def _ann(gene, **gd):
    return AnnotationData(
        gnomad=GnomADData(filter_pass=True, **gd),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        )],
    )


def _snv():
    return VariantRecord(chrom="chr17", pos=43044295, ref="A", alt="G", assembly=Assembly.GRCH38)


class TestRuleLoader:
    def test_min_depth_parsed(self, tmp_path):
        s = PM2Spec(_cfg(tmp_path).disease_prevalence_tsv)
        assert s.get("BRCA1").min_depth == 25.0


class TestDepthGate:
    def _ev(self, tmp_path, cov):
        ev = PM2Evaluator(_cfg(tmp_path))
        ev._coverage = cov
        return ev

    def test_absent_well_covered_fires(self, tmp_path):
        ev = self._ev(tmp_path, _Cov(True, 30.0))   # >= 25
        r = ev.evaluate(_snv(), _ann("BRCA1", af=0.0, ac=0))
        assert r.triggered

    def test_absent_low_coverage_blocked(self, tmp_path):
        ev = self._ev(tmp_path, _Cov(True, 10.0))   # < 25 → absence not callable
        r = ev.evaluate(_snv(), _ann("BRCA1", af=0.0, ac=0))
        assert not r.triggered and "depth" in r.evidence.lower()

    def test_coverage_unavailable_skips_gate(self, tmp_path):
        ev = self._ev(tmp_path, _Cov(False, None))
        r = ev.evaluate(_snv(), _ann("BRCA1", af=0.0, ac=0))
        assert r.triggered

    def test_depth_unknown_skips_gate(self, tmp_path):
        ev = self._ev(tmp_path, _Cov(True, None))   # no covered locus in span
        r = ev.evaluate(_snv(), _ann("BRCA1", af=0.0, ac=0))
        assert r.triggered


class TestExtraction:
    def _pm2(self, desc):
        return {"criteriaCodes": [{"label": "PM2", "evidenceStrengths": [
            {"label": "Supporting", "applicability": "Applicable", "description": desc},
        ]}]}

    def test_depth_25_detected(self):
        rs = self._pm2("Absent from controls. Region around the variant must have "
                       "an average read depth >=25.")
        assert bdt._pm2_min_depth(rs) == "25"

    def test_no_depth_blank(self):
        assert bdt._pm2_min_depth(self._pm2("Absent or extremely rare.")) == ""


def test_committed_brca_min_depth():
    tsv = _ROOT / "resources" / "shared" / "disease_prevalence.tsv"
    with tsv.open(encoding="utf-8") as f:
        rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
    assert rows["BRCA1"]["pm2_min_depth"] == "25"
    assert rows["BRCA2"]["pm2_min_depth"] == "25"


# --- builder + CoverageDB round-trip (needs duckdb; skipped where absent) -----

def _bgc():
    spec = importlib.util.spec_from_file_location(
        "build_gnomad_coverage", _ROOT / "scripts" / "build_gnomad_coverage.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class TestCoverageRoundTrip:
    def test_v2_layout(self, tmp_path):
        pytest.importorskip("duckdb")
        from acmg_classifier.local_db.coverage_db import CoverageDB
        bgz = tmp_path / "cov.tsv.bgz"
        with gzip.open(bgz, "wt", encoding="utf-8") as f:
            f.write("#chrom\tpos\tmean\tmedian\n17\t43044295\t30.5\t31\n17\t43044296\t10.0\t9\n")
        db = tmp_path / "cov.duckdb"
        assert _bgc().build(bgz, db) == 2
        cov = CoverageDB(db)
        assert cov.available
        assert cov.mean_depth("chr17", 43044295) == 30.5
        assert cov.mean_depth("17", 43044295, 43044296) == pytest.approx(20.25)
        assert cov.mean_depth("17", 99999999) is None

    def test_v4_locus_layout(self, tmp_path):
        pytest.importorskip("duckdb")
        from acmg_classifier.local_db.coverage_db import CoverageDB
        bgz = tmp_path / "cov4.tsv.bgz"
        with gzip.open(bgz, "wt", encoding="utf-8") as f:
            f.write("locus\tmean\tmedian_approx\nchr17:43044295\t40.0\t41\n")
        db = tmp_path / "cov4.duckdb"
        assert _bgc().build(bgz, db) == 1
        assert CoverageDB(db).mean_depth("17", 43044295) == 40.0
