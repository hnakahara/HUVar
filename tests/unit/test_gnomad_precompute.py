"""GnomADDB.precompute() must be byte-for-byte equivalent to per-variant query().

precompute() replaces the connection-per-variant access pattern with a single
JOIN, so its only contract is that cached() returns exactly what query() would
have. These tests build real tmp DuckDBs and assert equality across the tricky
cases: present / absent / filtered / multi-row merge / chrom-spelling / the
non-cancer companion backfill.
"""
from pathlib import Path

import pytest

duckdb = pytest.importorskip("duckdb")

from acmg_classifier.local_db.gnomad_db import GnomADDB  # noqa: E402
from acmg_classifier.models.enums import Assembly  # noqa: E402
from acmg_classifier.models.variant import VariantRecord  # noqa: E402

_NO_CONSTRAINT = Path("does-not-exist.tsv")

# Full modern schema (af_xy / ac_xx / grpmax / af_non_cancer all present) so the
# merge exercises every field.
_MAIN_COLS = (
    "chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, af DOUBLE, an INTEGER, "
    "ac INTEGER, nhomalt INTEGER, nhemi INTEGER, popmax_af DOUBLE, "
    "popmax_pop TEXT, faf95_popmax DOUBLE, af_xy DOUBLE, ac_xx INTEGER, "
    "nhomalt_xx INTEGER, filters TEXT, ac_grpmax INTEGER, an_grpmax INTEGER, "
    "af_non_cancer DOUBLE"
)


def _row(chrom, pos, ref, alt, af, popmax_af, faf95, *, nhomalt=0, pop="afr",
         filters="PASS", af_xy=None, ac_xx=None, nhomalt_xx=None,
         ac_grpmax=None, an_grpmax=None, af_non_cancer=None):
    return (chrom, pos, ref, alt, af, 100000, int(af * 100000), nhomalt, 0,
            popmax_af, pop, faf95, af_xy, ac_xx, nhomalt_xx, filters,
            ac_grpmax, an_grpmax, af_non_cancer)


def _main_db(tmp_path: Path, rows: list[tuple]) -> Path:
    p = tmp_path / "main.duckdb"
    con = duckdb.connect(str(p))
    con.execute(f"CREATE TABLE variants ({_MAIN_COLS})")
    if rows:
        placeholders = ",".join("?" * 19)
        con.executemany(f"INSERT INTO variants VALUES ({placeholders})", rows)
    con.close()
    return p


def _noncancer_db(tmp_path: Path, rows: list[tuple]) -> Path:
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


def _variant(chrom, pos, ref, alt) -> VariantRecord:
    return VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt,
                         assembly=Assembly.GRCH38)


def _assert_equivalent(db: GnomADDB, variants: list[VariantRecord]) -> None:
    """precompute() then cached() must equal a fresh query() for every variant."""
    expected = {
        v.key: db.query(v.chrom, v.pos, v.ref, v.alt) for v in variants
    }
    db.precompute(variants)
    for v in variants:
        got = db.cached(v)
        exp = expected[v.key]
        assert got == exp, f"mismatch for {v.key}: {got!r} != {exp!r}"


def test_present_absent_and_filtered(tmp_path):
    # "1"-style chrom in the DB, "chr"-normalised variants — exercises the
    # chrom-candidate join. Mix of present / absent / all-filtered rows.
    rows = [
        _row("1", 100, "A", "G", 0.001, 0.002, 0.0015, af_xy=0.0009,
             ac_xx=10, nhomalt_xx=1, ac_grpmax=200, an_grpmax=100000),
        _row("1", 200, "C", "T", 0.02, 0.03, 0.025, filters="AC0"),  # filtered
    ]
    db = GnomADDB(_main_db(tmp_path, rows), _NO_CONSTRAINT)
    variants = [
        _variant("chr1", 100, "A", "G"),   # present, PASS
        _variant("chr1", 200, "C", "T"),   # present, filter-failed
        _variant("chr1", 300, "G", "A"),   # absent → synthetic record
    ]
    _assert_equivalent(db, variants)


def test_multi_row_merge_equivalence(tmp_path):
    # Same variant twice (exomes+genomes style) → merge by per-field MAX.
    rows = [
        _row("chr17", 43000, "A", "G", 0.0003, 0.00032, 0.00022, nhomalt=2, pop="afr"),
        _row("chr17", 43000, "A", "G", 0.0001, 0.00015, 0.00012, nhomalt=9, pop="nfe"),
    ]
    db = GnomADDB(_main_db(tmp_path, rows), _NO_CONSTRAINT)
    _assert_equivalent(db, [_variant("chr17", 43000, "A", "G")])


def test_noncancer_backfill_equivalence(tmp_path):
    main = _main_db(tmp_path, [
        _row("17", 43000000, "A", "G", 2e-4, 2e-4, 1e-4),  # af_non_cancer NULL
    ])
    nc = _noncancer_db(tmp_path, [("17", 43000000, "A", "G", 1e-5, 8e-6)])
    db = GnomADDB(main, _NO_CONSTRAINT, nc)
    variants = [
        _variant("chr17", 43000000, "A", "G"),  # present → backfills subset
        _variant("chr17", 43999999, "C", "T"),  # absent → companion skipped
    ]
    _assert_equivalent(db, variants)


def test_missing_db_falls_back_to_query(tmp_path):
    # precompute() on a missing DB leaves the cache empty; cached() must fall
    # back to query() (which returns None) exactly as before.
    db = GnomADDB(tmp_path / "nope.duckdb", _NO_CONSTRAINT)
    v = _variant("chr1", 100, "A", "G")
    db.precompute([v])
    assert db.cached(v) is None
