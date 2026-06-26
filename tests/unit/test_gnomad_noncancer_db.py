"""GnomADDB non-cancer companion fallback (GRCh38 v3.1.2 overlay).

gnomAD v4.1 (GRCh38) dropped the non-cancer subset, so the main build's
af_non_cancer is always NULL. GnomADDB consults a small companion DB (built from
v3.1.2 genomes) to backfill af_non_cancer for a variant present in the main DB.
"""
from pathlib import Path

import pytest

duckdb = pytest.importorskip("duckdb")

from acmg_classifier.local_db.gnomad_db import GnomADDB  # noqa: E402

# Minimal `variants` schema = exactly the columns query() degrades against when
# af_xy / ac_xx / grpmax / af_non_cancer are absent (older-style build). With
# af_non_cancer omitted here, the main record's af_non_cancer is NULL → the
# companion fallback must supply it.
_MAIN_COLS = (
    "chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, af DOUBLE, an INTEGER, "
    "ac INTEGER, nhomalt INTEGER, nhemi INTEGER, popmax_af DOUBLE, "
    "popmax_pop TEXT, faf95_popmax DOUBLE, filters TEXT"
)


def _main_db(tmp_path: Path, rows: list[tuple]) -> Path:
    p = tmp_path / "main.duckdb"
    con = duckdb.connect(str(p))
    con.execute(f"CREATE TABLE variants ({_MAIN_COLS})")
    con.executemany(
        "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", rows
    )
    con.close()
    return p


def _noncancer_db(tmp_path: Path, rows: list[tuple]) -> Path:
    """Legacy companion schema (no faf95_non_cancer column) — exercises the
    schema-probe backward-compat path."""
    p = tmp_path / "non_cancer.duckdb"
    con = duckdb.connect(str(p))
    con.execute(
        "CREATE TABLE non_cancer (chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, "
        "af_non_cancer DOUBLE)"
    )
    if rows:  # DuckDB's executemany rejects an empty parameter list
        con.executemany("INSERT INTO non_cancer VALUES (?,?,?,?,?)", rows)
    con.close()
    return p


def _noncancer_db_v2(tmp_path: Path, rows: list[tuple]) -> Path:
    """Current companion schema, carrying the recomputed faf95_non_cancer."""
    p = tmp_path / "non_cancer.duckdb"
    con = duckdb.connect(str(p))
    con.execute(
        "CREATE TABLE non_cancer (chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, "
        "af_non_cancer DOUBLE, faf95_non_cancer DOUBLE)"
    )
    if rows:
        con.executemany("INSERT INTO non_cancer VALUES (?,?,?,?,?,?)", rows)
    con.close()
    return p


# A present BRCA1-region variant: overall af=2e-4, popmax_af=2e-4, filters PASS.
_PRESENT = ("17", 43000000, "A", "G", 2e-4, 100000, 20, 0, 0, 2e-4, "nfe", 1e-4, "PASS")
_NO_CONSTRAINT = Path("does-not-exist.tsv")


def test_backfills_af_non_cancer_from_companion(tmp_path):
    main = _main_db(tmp_path, [_PRESENT])
    nc = _noncancer_db(tmp_path, [("17", 43000000, "A", "G", 1e-5)])
    db = GnomADDB(main, _NO_CONSTRAINT, nc)
    gd = db.query("17", 43000000, "A", "G")
    assert gd is not None
    assert gd.af_non_cancer == pytest.approx(1e-5)
    # The overall AF is untouched — the companion only supplies the subset value.
    assert gd.af == pytest.approx(2e-4)


def test_backfills_faf95_non_cancer_from_companion(tmp_path):
    # Current companion schema: both af and the recomputed popmax FAF95 backfill.
    main = _main_db(tmp_path, [_PRESENT])
    nc = _noncancer_db_v2(tmp_path, [("17", 43000000, "A", "G", 1e-5, 8e-6)])
    db = GnomADDB(main, _NO_CONSTRAINT, nc)
    gd = db.query("17", 43000000, "A", "G")
    assert gd is not None
    assert gd.af_non_cancer == pytest.approx(1e-5)
    assert gd.faf95_non_cancer == pytest.approx(8e-6)


def test_legacy_companion_without_faf95_degrades_to_none(tmp_path):
    # A companion DB built before the faf95_non_cancer column: the schema probe
    # selects NULL for it, so faf95_non_cancer stays None (BA1/BS1 then fall back
    # to the overall FAF95) while af_non_cancer still backfills.
    main = _main_db(tmp_path, [_PRESENT])
    nc = _noncancer_db(tmp_path, [("17", 43000000, "A", "G", 1e-5)])
    db = GnomADDB(main, _NO_CONSTRAINT, nc)
    gd = db.query("17", 43000000, "A", "G")
    assert gd is not None
    assert gd.af_non_cancer == pytest.approx(1e-5)
    assert gd.faf95_non_cancer is None


def test_absent_in_companion_stays_none(tmp_path):
    main = _main_db(tmp_path, [_PRESENT])
    nc = _noncancer_db(tmp_path, [])  # variant not in the companion
    db = GnomADDB(main, _NO_CONSTRAINT, nc)
    gd = db.query("17", 43000000, "A", "G")
    assert gd is not None and gd.af_non_cancer is None


def test_no_companion_configured_is_noop(tmp_path):
    main = _main_db(tmp_path, [_PRESENT])
    db = GnomADDB(main, _NO_CONSTRAINT)  # GRCh37-style: no companion
    gd = db.query("17", 43000000, "A", "G")
    assert gd is not None and gd.af_non_cancer is None


def test_variant_absent_from_main_skips_companion(tmp_path):
    # Absent from the main DB (no row) → synthetic absent record; the companion
    # is NOT consulted (an absent variant is absent in every subset), so PM2's
    # own "absent" path handles it. Documents the rows-present restriction.
    main = _main_db(tmp_path, [_PRESENT])
    nc = _noncancer_db(tmp_path, [("17", 43999999, "C", "T", 5e-5)])
    db = GnomADDB(main, _NO_CONSTRAINT, nc)
    gd = db.query("17", 43999999, "C", "T")
    assert gd is not None and gd.ac == 0 and gd.af_non_cancer is None
