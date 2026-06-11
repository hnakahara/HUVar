"""Pinpoint why a PM5 same-codon query returns no ClinVar comparators.

PM5 looks for OTHER pathogenic missense at the same codon via
``WHERE gene_symbol=? AND codon_position=? AND clinical_significance IN (P/LP)
AND star_rating>=?``. When that returns nothing despite a known comparator
existing in ClinVar (e.g. TP53 p.Arg248Gln for a p.Arg248Trp candidate), one of
those four stored fields is wrong. This script dumps exactly what the built DB
holds so we can see which one.

Usage (run against the SAME data dir the classifier uses):
    python scripts/diagnose_clinvar_pm5.py --gene TP53 --codon 248
    python scripts/diagnose_clinvar_pm5.py --gene TP53 --codon 248 \
        --db data/GRCh38/clinvar/clinvar_ps1_pm5_GRCh38.sqlite
    # optionally pin a known ClinVar VariationID to look up directly:
    python scripts/diagnose_clinvar_pm5.py --gene TP53 --codon 248 --variation-id 12356
"""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

_P_LP = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")


def _connect(db: Path) -> sqlite3.Connection:
    uri = "file:" + db.as_posix() + "?mode=ro&immutable=1"
    return sqlite3.connect(uri, uri=True)


def _cols(con: sqlite3.Connection) -> set[str]:
    return {r[1] for r in con.execute("PRAGMA table_info(variants)").fetchall()}


def _print_rows(rows: list[tuple], header: list[str]) -> None:
    if not rows:
        print("    (none)")
        return
    for r in rows:
        print("    " + " | ".join(f"{h}={v!r}" for h, v in zip(header, r)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="data/GRCh38/clinvar/clinvar_ps1_pm5_GRCh38.sqlite")
    ap.add_argument("--gene", required=True)
    ap.add_argument("--codon", type=int, required=True)
    ap.add_argument("--min-stars", type=int, default=1)
    ap.add_argument("--variation-id", default=None)
    args = ap.parse_args()

    db = Path(args.db)
    if not db.exists():
        raise SystemExit(f"DB not found: {db.resolve()}")
    con = _connect(db)
    cols = _cols(con)
    print(f"DB: {db.resolve()}")
    print(f"Columns: {sorted(cols)}\n")

    # 1) The EXACT PM5 query (minus the self-exclusion LIKE). If this is empty,
    #    PM5 cannot fire — the rest of the script shows which clause killed it.
    sig_list = ",".join("?" for _ in _P_LP)
    print(f"[1] Exact PM5 match  gene={args.gene} codon={args.codon} "
          f"P/LP star>={args.min_stars}:")
    q1 = con.execute(
        f"""SELECT variation_id, gene_symbol, hgvs_p, codon_position,
                   clinical_significance, star_rating
            FROM variants
            WHERE gene_symbol = ? AND codon_position = ? AND star_rating >= ?
              AND clinical_significance IN ({sig_list})""",
        (args.gene, args.codon, args.min_stars, *_P_LP),
    ).fetchall()
    _print_rows(q1, ["id", "gene", "hgvs_p", "codon", "sig", "stars"])

    # 2) Same codon, ANY gene / ANY significance / ANY stars — does the residue
    #    exist at all, and under what gene_symbol / significance / star?
    print(f"\n[2] ANY row at codon_position={args.codon} (any gene/sig/star):")
    q2 = con.execute(
        """SELECT variation_id, gene_symbol, hgvs_p, codon_position,
                  clinical_significance, star_rating
           FROM variants WHERE codon_position = ?
           ORDER BY gene_symbol LIMIT 40""",
        (args.codon,),
    ).fetchall()
    _print_rows(q2, ["id", "gene", "hgvs_p", "codon", "sig", "stars"])

    # 3) All rows for the gene with codon_position NULL but an hgvs_p that names
    #    this codon — proves a codon_position PARSE failure (the prime suspect).
    print(f"\n[3] {args.gene} rows whose hgvs_p mentions {args.codon} but "
          f"codon_position IS NULL (parse failures):")
    q3 = con.execute(
        """SELECT variation_id, gene_symbol, hgvs_p, codon_position,
                  clinical_significance, star_rating
           FROM variants
           WHERE gene_symbol = ? AND codon_position IS NULL
             AND hgvs_p LIKE ? LIMIT 40""",
        (args.gene, f"%{args.codon}%"),
    ).fetchall()
    _print_rows(q3, ["id", "gene", "hgvs_p", "codon", "sig", "stars"])

    # 4) Distinct gene_symbols actually stored for this codon's hgvs_p — reveals
    #    an overlapping-locus mis-assignment (e.g. TP53 stored under WRAP53).
    print(f"\n[4] Distinct gene_symbols on rows whose hgvs_p mentions "
          f"'{args.codon}':")
    q4 = con.execute(
        """SELECT gene_symbol, COUNT(*) FROM variants
           WHERE hgvs_p LIKE ? GROUP BY gene_symbol ORDER BY 2 DESC LIMIT 20""",
        (f"%{args.codon}%",),
    ).fetchall()
    _print_rows(q4, ["gene", "n"])

    # 5) Direct lookup of a known comparator VariationID, if provided.
    if args.variation_id:
        print(f"\n[5] VariationID {args.variation_id} as stored:")
        q5 = con.execute(
            """SELECT variation_id, gene_symbol, hgvs_p, codon_position,
                      clinical_significance, review_status, star_rating
               FROM variants WHERE variation_id = ?""",
            (str(args.variation_id),),
        ).fetchall()
        _print_rows(q5, ["id", "gene", "hgvs_p", "codon", "sig", "review", "stars"])

    con.close()
    print("\nInterpretation:")
    print("  [1] non-empty  -> PM5 SHOULD fire; the bug is in the evaluator/annotation,")
    print("                    not the DB (re-check pc.gene_symbol / pc.protein_position).")
    print("  [1] empty but [2] shows the variant under a DIFFERENT gene_symbol")
    print("                 -> overlapping-locus mis-assignment in _gene_symbol.")
    print("  [3] non-empty  -> codon_position parse failure (hgvs_p format the")
    print("                    builder's _parse_aa_change does not handle).")
    print("  [2] shows sig NOT in P/LP, or stars < min -> significance/star issue.")


if __name__ == "__main__":
    main()
