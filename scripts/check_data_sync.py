#!/usr/bin/env python3
"""Pre-run check of a data directory against this checkout.

    python scripts/check_data_sync.py /path/to/data [--assembly GRCh38]

Reports:
  * the curated tables in ``resources/shared/*.tsv`` whose copy in
    ``<data>/shared/`` is missing or different. The evaluators read the
    ``data/shared`` copies, so a stale copy hides rule changes in the repository
    (refresh with ``cp resources/shared/*.tsv <data>/shared/``);
  * whether the gnomAD DuckDB carries the GrpMax AC / AN columns (used by the
    subpopulation PM2 rules, the BA1 / BS1 minimum allele counts and the PM2
    population-size rule; without them those conditions are skipped);
  * whether the gnomAD coverage DB exists (read-depth conditions of PM2 / BA1 /
    BS1; skipped without it).
"""
from __future__ import annotations

import argparse
import filecmp
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("--assembly", default="GRCh38", choices=("GRCh38", "GRCh37"))
    args = ap.parse_args()
    from acmg_classifier.config import Config
    from acmg_classifier.models.enums import Assembly
    cfg = Config(data_dir=args.data_dir, assembly=Assembly(args.assembly))

    shared = args.data_dir / "shared"
    stale = [f.name for f in sorted((ROOT / "resources" / "shared").glob("*.tsv"))
             if not (shared / f.name).exists()
             or not filecmp.cmp(f, shared / f.name, shallow=False)]
    print(f"data/shared tables missing or different from the repo ({len(stale)}):")
    for name in stale:
        print(f"  {name}")
    if stale:
        print(f"  -> cp {ROOT / 'resources' / 'shared'}/*.tsv {shared}/")

    g = cfg.gnomad_duckdb
    if g.exists():
        import duckdb
        con = duckdb.connect(str(g), read_only=True)
        cols = {r[0] for r in con.execute("DESCRIBE variants").fetchall()}
        con.close()
        ok = {"ac_grpmax", "an_grpmax"} <= cols
        print(f"gnomAD GrpMax AC/AN columns: {'yes' if ok else 'NO (rebuild the gnomAD step)'}")
    else:
        print(f"gnomAD DB not found: {g}")
    cov = cfg.gnomad_coverage_db
    print(f"gnomAD coverage DB: {'yes' if cov.exists() else f'NO ({cov}); --only gnomad-coverage'}")
    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())
