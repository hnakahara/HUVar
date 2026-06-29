"""Per-CSpec multispec overlay: RYR1/ACTA1/VWF disease-specific evaluation.

These cover the data path that lets the app evaluate a variant under a chosen
ClinGen CSpec without touching CLI/batch behaviour:

* ``cspec_overlay.available_cspecs`` / ``overlay_tsv`` (pure overlay logic), and
* the committed ``disease_prevalence_multispec.tsv`` (built by
  ``build_disease_thresholds.py --multispec-out``) actually flips a gene's
  applicability across CSpecs when read by a real threshold loader.
"""
import csv
from pathlib import Path

import pytest

from acmg_classifier.criteria.cspec_overlay import available_cspecs, overlay_tsv
from acmg_classifier.criteria.pvs1_genes import PVS1Applicability
from acmg_classifier.criteria.allele_frequency import DiseaseThresholds

_REPO_ROOT = Path(__file__).resolve().parents[2]
_MULTISPEC = _REPO_ROOT / "resources" / "shared" / "disease_prevalence_multispec.tsv"

_BASE = (
    "gene_symbol\tba1_threshold\tbs1_threshold\tpvs1\n"
    "RYR1\t0.00697\t0.000697\tnot_applicable\n"
    "TP53\t0.001\t0.0001\tapplicable\n"
)
_MS = (
    "cspec_id\tdisease_label\tsource_gn\tgene_symbol\tba1_threshold\tbs1_threshold\tpvs1\n"
    "malignant_hyperthermia\t悪性高熱症\tGN012\tRYR1\t0.0038\t0.0008\tnot_applicable\n"
    "congenital_myopathy_ar\t先天性ミオパチー(AR)\tGN179\tRYR1\t0.00697\t0.000697\tapplicable\n"
)


def _synthetic(tmp_path: Path) -> tuple[Path, Path]:
    base = tmp_path / "base.tsv"
    base.write_text(_BASE, encoding="utf-8")
    ms = tmp_path / "multispec.tsv"
    ms.write_text(_MS, encoding="utf-8")
    return base, ms


class TestAvailableCspecs:
    def test_lists_in_table_order_with_labels(self, tmp_path):
        _, ms = _synthetic(tmp_path)
        cs = available_cspecs("RYR1", ms)
        assert [c["cspec_id"] for c in cs] == [
            "malignant_hyperthermia", "congenital_myopathy_ar"]
        assert cs[0]["label"] == "悪性高熱症"
        assert cs[0]["source_gn"] == "GN012"

    def test_gene_without_cspecs_is_empty(self, tmp_path):
        _, ms = _synthetic(tmp_path)
        assert available_cspecs("TP53", ms) == []

    def test_none_gene_and_missing_file_are_empty(self, tmp_path):
        _, ms = _synthetic(tmp_path)
        assert available_cspecs(None, ms) == []
        assert available_cspecs("RYR1", tmp_path / "absent.tsv") == []


class TestOverlayTsv:
    def test_replaces_only_target_gene(self, tmp_path):
        base, ms = _synthetic(tmp_path)
        out = overlay_tsv(base, ms, "RYR1", "congenital_myopathy_ar")
        try:
            with out.open(encoding="utf-8") as fh:
                rows = {r["gene_symbol"]: r
                        for r in csv.DictReader(fh, delimiter="\t")}
        finally:
            out.unlink()
        assert rows["RYR1"]["pvs1"] == "applicable"
        assert rows["RYR1"]["ba1_threshold"] == "0.00697"
        assert rows["TP53"]["pvs1"] == "applicable"  # other genes untouched

    def test_overlay_columns_restricted_to_base(self, tmp_path):
        base, ms = _synthetic(tmp_path)
        out = overlay_tsv(base, ms, "RYR1", "malignant_hyperthermia")
        try:
            with out.open(encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                assert reader.fieldnames == ["gene_symbol", "ba1_threshold",
                                             "bs1_threshold", "pvs1"]
                row = next(r for r in reader if r["gene_symbol"] == "RYR1")
        finally:
            out.unlink()
        # multispec-only columns must not leak into the overlay
        assert "cspec_id" not in row
        assert row["ba1_threshold"] == "0.0038"

    def test_unknown_cspec_raises(self, tmp_path):
        base, ms = _synthetic(tmp_path)
        with pytest.raises(KeyError):
            overlay_tsv(base, ms, "RYR1", "does_not_exist")

    def test_overlay_flips_loader_applicability(self, tmp_path):
        # The core promise: same gene, different CSpec → different PVS1/BA1 when a
        # real loader reads the overlay table.
        base, ms = _synthetic(tmp_path)
        mh = overlay_tsv(base, ms, "RYR1", "malignant_hyperthermia")
        ar = overlay_tsv(base, ms, "RYR1", "congenital_myopathy_ar")
        try:
            assert PVS1Applicability(mh).is_not_applicable("RYR1") is True
            assert PVS1Applicability(ar).is_not_applicable("RYR1") is False
            assert DiseaseThresholds(mh).get("RYR1").ba1 == 0.0038
            assert DiseaseThresholds(ar).get("RYR1").ba1 == 0.00697
        finally:
            mh.unlink()
            ar.unlink()


@pytest.mark.skipif(not _MULTISPEC.exists(),
                    reason="multispec table not built (run build_disease_thresholds --multispec-out)")
class TestCommittedMultispecTable:
    """The committed multispec table built from the real ClinGen specs must carry
    the curated CSpecs and the disease-distinguishing applicability flips."""

    def test_curated_cspecs_present(self):
        ryr1 = {c["cspec_id"] for c in available_cspecs("RYR1", _MULTISPEC)}
        assert ryr1 == {"malignant_hyperthermia",
                        "congenital_myopathy_ad", "congenital_myopathy_ar"}
        assert {c["cspec_id"] for c in available_cspecs("ACTA1", _MULTISPEC)} == {
            "congenital_myopathy_ad", "congenital_myopathy_ar"}
        assert {c["cspec_id"] for c in available_cspecs("VWF", _MULTISPEC)} == {
            "vwd_ad", "vwd_ar"}

    def test_ryr1_pvs1_differs_by_cspec(self, tmp_path):
        base = tmp_path / "base.tsv"
        base.write_text(_BASE, encoding="utf-8")
        mh = overlay_tsv(base, _MULTISPEC, "RYR1", "malignant_hyperthermia")
        ar = overlay_tsv(base, _MULTISPEC, "RYR1", "congenital_myopathy_ar")
        try:
            # Malignant hyperthermia (GoF) declines PVS1; recessive congenital
            # myopathy (LoF) applies it.
            assert PVS1Applicability(mh).is_not_applicable("RYR1") is True
            assert PVS1Applicability(ar).is_not_applicable("RYR1") is False
        finally:
            mh.unlink()
            ar.unlink()
