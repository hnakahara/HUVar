"""Unit tests for the Brandes 2023 ESM1b matrix parser and SQLite builder."""
import sqlite3
import zipfile
from pathlib import Path

from acmg_classifier.local_db.esm1b_db import query_esm1b
from acmg_classifier.setup.esm1b_builder import (
    _iter_aliases,
    _iter_matrix_rows,
    _parse_header,
    _uniprot_from_name,
    build_esm1b_sqlite,
)


# Brandes isoform_list.csv: `id,txt` where txt is "GENE (MNEMONIC) | ACCESSION".
_SAMPLE_ISOFORM_LIST = (
    "id,txt\n"
    "P12345,PIK3CD (PK3CD) | P12345\n"
    "A0A024RBG1,NUDT4B (NUD4B) | A0A024RBG1\n"
    "Q5SR53,  (CA200) | Q5SR53\n"
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


def test_iter_aliases_parses_mnemonic_to_entry_name(tmp_path: Path):
    (tmp_path / "isoform_list.csv").write_text(_SAMPLE_ISOFORM_LIST, encoding="utf-8")
    zip_path = tmp_path / "ALL_hum_isoforms_ESM1b_LLR.zip"
    zip_path.touch()  # _iter_aliases prefers the sidecar CSV next to the zip

    aliases = dict(_iter_aliases(zip_path))
    assert aliases["PK3CD_HUMAN"] == "P12345"
    assert aliases["NUD4B_HUMAN"] == "A0A024RBG1"
    # Empty gene symbol before the parens must not break mnemonic extraction.
    assert aliases["CA200_HUMAN"] == "Q5SR53"


def test_iter_aliases_missing_file_yields_nothing(tmp_path: Path):
    zip_path = tmp_path / "ALL_hum_isoforms_ESM1b_LLR.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("P12345_LLR.csv", _SAMPLE_CSV)  # no isoform_list.csv anywhere
    assert list(_iter_aliases(zip_path)) == []


def _build_db_with_aliases(tmp_path: Path) -> Path:
    (tmp_path / "isoform_list.csv").write_text(_SAMPLE_ISOFORM_LIST, encoding="utf-8")
    zip_path = tmp_path / "ALL_hum_isoforms_ESM1b_LLR.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("P12345_LLR.csv", _SAMPLE_CSV)
    dest = tmp_path / "esm1b.sqlite"
    build_esm1b_sqlite(zip_path, dest, batch_size=5)
    return dest


def test_build_populates_aliases_table(tmp_path: Path):
    dest = _build_db_with_aliases(tmp_path)
    conn = sqlite3.connect(str(dest))
    try:
        row = conn.execute(
            "SELECT uniprot_id FROM aliases WHERE entry_name = ?", ("PK3CD_HUMAN",)
        ).fetchone()
        assert row is not None and row[0] == "P12345"
    finally:
        conn.close()


def test_query_esm1b_resolves_swissprot_entry_name(tmp_path: Path):
    """GRCh37 cache emits 'PK3CD_HUMAN'; it must resolve to accession P12345."""
    dest = _build_db_with_aliases(tmp_path)
    res = query_esm1b(dest, "PK3CD_HUMAN", 1, "R")
    assert res is not None
    assert res.llr == -12.401


def test_query_esm1b_accepts_plain_accession(tmp_path: Path):
    """GRCh38 cache emits the accession directly — must pass through unchanged."""
    dest = _build_db_with_aliases(tmp_path)
    res = query_esm1b(dest, "P12345", 1, "R")
    assert res is not None
    assert res.llr == -12.401


def test_query_esm1b_strips_human_for_trembl_without_alias(tmp_path: Path):
    """TrEMBL entry names are '<accession>_HUMAN' — strip when no alias row."""
    dest = _build_db_with_aliases(tmp_path)
    # 'P12345_HUMAN' has no alias row, so it strips to accession 'P12345'.
    res = query_esm1b(dest, "P12345_HUMAN", 1, "R")
    assert res is not None
    assert res.llr == -12.401
