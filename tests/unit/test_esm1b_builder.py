"""Unit tests for the Brandes 2023 ESM1b matrix parser and SQLite builder."""
import sqlite3
import zipfile
from pathlib import Path

from acmg_classifier.setup.esm1b_builder import (
    _iter_matrix_rows,
    _parse_header,
    _uniprot_from_name,
    build_esm1b_sqlite,
)


# Minimal Brandes-format CSV: 4 positions (M K E L), 3 alt rows (K R A).
_SAMPLE_CSV = (
    ",M 1,K 2,E 3,L 4\n"
    "K,-11.210,0.000,-7.500,-4.500\n"
    "R,-12.401,-3.200,-6.500,-5.000\n"
    "A,-9.000,-6.000,-3.500,-2.000\n"
)


def test_parse_header_extracts_wt_position_pairs():
    info = _parse_header([",", "M 1", "K 2", "E 3", "L 4"])
    assert info == [("M", 1), ("K", 2), ("E", 3), ("L", 4)]


def test_parse_header_skips_malformed_cells():
    info = _parse_header([",", "M 1", "junk", "E 3"])
    assert info == [("M", 1), None, ("E", 3)]


def test_iter_matrix_rows_skips_wt_diagonal():
    rows = list(_iter_matrix_rows("P12345", _SAMPLE_CSV))
    # K's diagonal at position 2 (WT=K) is 0.000 — skipped.
    assert ("P12345", 2, "K", 0.0) not in rows
    # All other 11 cells should be present (3 rows × 4 positions − 1 WT cell).
    assert len(rows) == 11


def test_iter_matrix_rows_yields_correct_llr_values():
    rows = list(_iter_matrix_rows("P12345", _SAMPLE_CSV))
    by_key = {(r[1], r[2]): r[3] for r in rows}
    assert by_key[(1, "K")] == -11.210
    assert by_key[(1, "R")] == -12.401
    assert by_key[(3, "A")] == -3.500
    assert by_key[(4, "K")] == -4.500


def test_uniprot_from_name_handles_canonical_and_isoform():
    assert _uniprot_from_name("P38398_LLR.csv") == "P38398"
    assert _uniprot_from_name("P38398-2_LLR.csv") == "P38398-2"
    assert _uniprot_from_name("subdir/Q14524_LLR.csv") == "Q14524"


def test_uniprot_from_name_rejects_unexpected_names():
    assert _uniprot_from_name("README.csv") is None
    assert _uniprot_from_name("notes.txt") is None


def test_build_esm1b_sqlite_end_to_end(tmp_path: Path):
    zip_path = tmp_path / "ALL_hum_isoforms_ESM1b_LLR.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("P12345_LLR.csv", _SAMPLE_CSV)
        zf.writestr("P12345-2_LLR.csv", _SAMPLE_CSV)

    dest = tmp_path / "esm1b.sqlite"
    build_esm1b_sqlite(zip_path, dest, batch_size=5)

    conn = sqlite3.connect(str(dest))
    try:
        # Two isoforms × 11 rows each = 22 rows.
        total = conn.execute("SELECT COUNT(*) FROM scores").fetchone()[0]
        assert total == 22

        # Spot-check a specific lookup.
        row = conn.execute(
            "SELECT llr FROM scores "
            "WHERE uniprot_id = ? AND aa_pos = ? AND alt_aa = ?",
            ("P12345", 1, "R"),
        ).fetchone()
        assert row is not None
        assert row[0] == -12.401

        # Isoform suffix must be preserved.
        row = conn.execute(
            "SELECT llr FROM scores WHERE uniprot_id = ? AND aa_pos = ?",
            ("P12345-2", 4),
        ).fetchone()
        assert row is not None
    finally:
        conn.close()
