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
    """The header column names **exactly as read_csv exposes them** — the leading
    "#" of a "#chrom" header is kept, because DuckDB does not strip it."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        first = fh.readline().rstrip("\n")
    return first.split("\t")


def _select_sql(src: str, cols: list[str]) -> str:
    """A SELECT that normalises either layout to (chrom, pos, mean_dp).

    Each source column is referenced by its **actual** name (qualified with the
    table alias ``t``): the output aliases chrom/pos/mean_dp would otherwise
    collide with the source column names, and newer DuckDB then treats the bare
    column as a forward lateral reference to the not-yet-defined output alias and
    errors. The chrom header may be "#chrom" (the "#" is retained by read_csv)."""
    # normalised (lower-case, "#"-stripped) name -> actual column name
    norm = {c.lower().lstrip("#"): c for c in cols}
    read = (
        f"read_csv('{src}', delim='\t', header=true, compression='gzip', "
        "ignore_errors=true, auto_detect=true) AS t"
    )
    mean_col = norm.get("mean", "mean")
    if "locus" in norm:
        # v4.x: locus = "chr1:12345"
        locus_col = norm["locus"]
        return (
            f"SELECT split_part(t.\"{locus_col}\", ':', 1) AS chrom, "
            f"CAST(split_part(t.\"{locus_col}\", ':', 2) AS INTEGER) AS pos, "
            f'CAST(t."{mean_col}" AS DOUBLE) AS mean_dp '
            f'FROM {read} WHERE t."{mean_col}" IS NOT NULL'
        )
    # v2.x: separate chrom + pos columns (chrom header may be "#chrom").
    chrom_col = norm.get("chrom", cols[0])
    pos_col = norm.get("pos", "pos")
    return (
        f'SELECT CAST(t."{chrom_col}" AS VARCHAR) AS chrom, '
        f'CAST(t."{pos_col}" AS INTEGER) AS pos, '
        f'CAST(t."{mean_col}" AS DOUBLE) AS mean_dp '
        f'FROM {read} WHERE t."{mean_col}" IS NOT NULL'
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
