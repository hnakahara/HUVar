"""Reviewer 2, Comment 6: --cspec GENE=CSPEC_ID selects a disease-specific CSpec
(RYR1 / ACTA1 / VWF) in the CLI, as HUVar_app does."""
import csv
from pathlib import Path
from types import SimpleNamespace

import pytest

from acmg_classifier.config import Config
from acmg_classifier.criteria.cspec_overlay import apply_cspec_selection, parse_cspec_options

_RES = Path(__file__).resolve().parents[2] / "resources" / "shared"


def _row(path, gene):
    with Path(path).open(encoding="utf-8") as fh:
        r = next(r for r in csv.DictReader(fh, delimiter="\t") if r["gene_symbol"] == gene)
    return {k: v or "" for k, v in r.items()}  # ragged rows: None == ""


def _cfg(tmp_path):
    (tmp_path / "shared").mkdir()
    (tmp_path / "shared" / "disease_prevalence.tsv").write_bytes(
        (_RES / "disease_prevalence.tsv").read_bytes())
    return Config(data_dir=tmp_path)


def test_parse_and_errors():
    ms = _RES / "disease_prevalence_multispec.tsv"
    assert parse_cspec_options(["RYR1=congenital_myopathy_ar", "VWF=vwd_ad"], ms) == {
        "RYR1": "congenital_myopathy_ar", "VWF": "vwd_ad"}
    with pytest.raises(ValueError, match="available"):
        parse_cspec_options(["RYR1=nope"], ms)
    with pytest.raises(ValueError, match="single specification"):
        parse_cspec_options(["BRCA1=x"], ms)
    with pytest.raises(ValueError, match="GENE=CSPEC_ID"):
        parse_cspec_options(["RYR1"], ms)


def test_overlay_replaces_selected_rows_only(tmp_path):
    cfg = _cfg(tmp_path)
    assert cfg.disease_prevalence_multispec_tsv.exists()  # bundled fallback
    new = apply_cspec_selection(cfg, {"RYR1": "malignant_hyperthermia", "VWF": "vwd_ar"})
    ms = cfg.disease_prevalence_multispec_tsv
    for gene, cs in (("RYR1", "malignant_hyperthermia"), ("VWF", "vwd_ar")):
        want = next(r for r in csv.DictReader(ms.open(encoding="utf-8"), delimiter="\t")
                    if r["gene_symbol"] == gene and r["cspec_id"] == cs)
        got = _row(new.disease_prevalence_tsv, gene)
        assert got["bs1_threshold"] == want["bs1_threshold"]
    assert _row(new.disease_prevalence_tsv, "BRCA1") == _row(cfg.disease_prevalence_tsv, "BRCA1")
    assert new.cspec_selection == {"RYR1": "malignant_hyperthermia", "VWF": "vwd_ar"}
    assert apply_cspec_selection(cfg, {}) is cfg


def test_cspec_note():
    from acmg_classifier.pipeline.pipeline import _cspec_note
    cfg = SimpleNamespace(cspec_selection={"RYR1": "congenital_myopathy_ar"})
    assert "congenital_myopathy_ar" in _cspec_note(cfg, SimpleNamespace(gene_symbol="RYR1"))[0]
    assert _cspec_note(cfg, SimpleNamespace(gene_symbol="ACTA1")) == []


def test_cli_option_listed():
    from click.testing import CliRunner
    from acmg_classifier.cli import cli
    for cmd in ("classify", "explain"):
        assert "--cspec" in CliRunner().invoke(cli, [cmd, "--help"]).output
