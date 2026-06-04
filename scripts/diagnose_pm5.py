#!/usr/bin/env python3
"""Diagnose why PM5 finds no same-codon comparator for a variant.

PM5 matches a ClinVar P/LP missense at the SAME protein residue (the SQLite
``codon_position`` column, which the ClinVar builder derives from the stored
protein HGVS). A false-negative PM5 means that lookup returned nothing. This
script dumps every ClinVar row this build holds for a gene at (and near) a
residue, so you can see whether the expected comparator is absent, has a NULL /
wrong ``codon_position`` (an HGVS-parse failure), too few stars, or a
non-P/LP significance.

Usage:
    python scripts/diagnose_pm5.py --db data/GRCh38/clinvar/clinvar_ps1_pm5_GRCh38.sqlite \
        --gene PIK3CD --codon 524
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, type=Path, help="clinvar_ps1_pm5_*.sqlite path")
    ap.add_argument("--gene", required=True)
    ap.add_argument("--codon", required=True, type=int)
    ap.add_argument("--window", type=int, default=2,
                    help="also show +/- this many residues (default 2)")
    args = ap.parse_args()

    con = sqlite3.connect(f"file:{args.db.as_posix()}?mode=ro", uri=True)

    # 1) Exact-codon rows (what PM5 actually matches on).
    print(f"== {args.gene} codon_position = {args.codon} (PM5 match key) ==")
    rows = con.execute(
        "SELECT variation_id, hgvs_p, amino_acid_change, codon_position, "
        "clinical_significance, review_status, star_rating "
        "FROM variants WHERE gene_symbol=? AND codon_position=? "
        "ORDER BY star_rating DESC",
        (args.gene, args.codon),
    ).fetchall()
    _dump(rows)
    plp = [r for r in rows if (r[4] or "") in
           ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
           and (r[6] or 0) >= 2]
    print(f"   -> {len(plp)} row(s) qualify as a >=2 star P/LP comparator\n")

    # 2) Same gene, any row whose HGVS text mentions the residue number but whose
    #    parsed codon_position is NULL or different — reveals HGVS-parse misses.
    print(f"== {args.gene} rows mentioning '{args.codon}' in HGVS but codon_position != {args.codon} ==")
    near = con.execute(
        "SELECT variation_id, hgvs_p, amino_acid_change, codon_position, "
        "clinical_significance, review_status, star_rating "
        "FROM variants WHERE gene_symbol=? AND (hgvs_p LIKE ? OR hgvs_c LIKE ?) "
        "AND (codon_position IS NULL OR codon_position<>?) "
        "ORDER BY star_rating DESC LIMIT 40",
        (args.gene, f"%{args.codon}%", f"%{args.codon}%", args.codon),
    ).fetchall()
    _dump(near)

    # 3) Coverage sanity: how many of this gene's rows have a NULL codon_position?
    total, nullc = con.execute(
        "SELECT COUNT(*), SUM(CASE WHEN codon_position IS NULL THEN 1 ELSE 0 END) "
        "FROM variants WHERE gene_symbol=?",
        (args.gene,),
    ).fetchone()
    print(f"\n== {args.gene}: {nullc}/{total} rows have NULL codon_position "
          f"({(nullc or 0) * 100 // (total or 1)}%) ==")
    con.close()


def _dump(rows: list[tuple]) -> None:
    if not rows:
        print("   (none)")
        return
    for vid, hp, aa, cp, sig, rs, st in rows:
        print(f"   id={vid} star={st} cp={cp!s:>5} sig={sig!r} aa={aa} hgvs_p={hp!r}")


if __name__ == "__main__":
    main()
