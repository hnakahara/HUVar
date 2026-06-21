#!/usr/bin/env python3
"""Build a per-locus coverage DuckDB from the gnomAD exomes coverage summary.

Input: the gnomAD coverage summary TSV (bgz), downloaded by
``scripts/setup_data.py --with-gnomad-coverage``. The column layout differs by
release:
  * v2.1 exomes: ``#chrom  pos  mean  median  over_1 ...`` (chrom without "chr").
  * v4.0 exomes: ``locus  mean  median_approx  ...`` (locus = "chr1:12345").

Output: a DuckDB with one table ``coverage(chrom TEXT, pos INTEGER,
mean_dp DOUBLE)`` indexed on (chrom, pos) — consumed by
:class:`acmg_classifier.local_db.coverage_db.CoverageDB` for the PM2 read-depth
gate.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def _header_cols(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        first = fh.readline().rstrip("\n")
    return first.lstrip("#").split("\t")


def _select_sql(src: str, cols: list[str]) -> str:
    """A SELECT that normalises either layout to (chrom, pos, mean_dp)."""
    lower = [c.lower() for c in cols]
    read = (
        f"read_csv('{src}', delim='\t', header=true, compression='gzip', "
        "ignore_errors=true, auto_detect=true)"
    )
    if "locus" in lower:
        # v4.x: locus = "chr1:12345"
        return (
            "SELECT split_part(locus, ':', 1) AS chrom, "
            "CAST(split_part(locus, ':', 2) AS INTEGER) AS pos, "
            "CAST(mean AS DOUBLE) AS mean_dp "
            f"FROM {read} WHERE mean IS NOT NULL"
        )
    # v2.x: separate chrom + pos columns (the chrom header may be "#chrom",
    # which read_csv exposes as "chrom" after the '#').
    chrom_col = "chrom" if "chrom" in lower else cols[0]
    return (
        f'SELECT CAST("{chrom_col}" AS VARCHAR) AS chrom, '
        'CAST(pos AS INTEGER) AS pos, CAST(mean AS DOUBLE) AS mean_dp '
        f"FROM {read} WHERE mean IS NOT NULL"
    )


def build(coverage_tsv: Path, out_db: Path) -> int:
    import duckdb

    cols = _header_cols(coverage_tsv)
    out_db.parent.mkdir(parents=True, exist_ok=True)
    if out_db.exists():
        out_db.unlink()
    con = duckdb.connect(str(out_db))
    src = str(coverage_tsv).replace("\\", "/")
    con.execute(f"CREATE TABLE coverage AS {_select_sql(src, cols)}")
    con.execute("CREATE INDEX idx_cov ON coverage(chrom, pos)")
    n = con.execute("SELECT count(*) FROM coverage").fetchone()[0]
    con.close()
    return int(n)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coverage-tsv", required=True,
                    help="gnomAD coverage summary .tsv.bgz")
    ap.add_argument("--out", required=True, help="output coverage DuckDB")
    args = ap.parse_args()
    n = build(Path(args.coverage_tsv), Path(args.out))
    print(f"coverage rows: {n} | written → {args.out}")


if __name__ == "__main__":
    main()
