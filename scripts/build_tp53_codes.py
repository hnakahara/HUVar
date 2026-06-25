#!/usr/bin/env python3
"""Build resources/shared/tp53_pp3_bp4_codes.tsv from the ClinGen TP53 VCEP's
precomputed PP3/BP4 code spreadsheet.

The TP53 VCEP distributes "Supplementary Table S2 — Bioinformatic predictions and
corresponding PP3 and BP4 codes for every possible p53 missense variant using the
TP53 specifications v2". Each row carries the variant's Align-GVGD class, BayesDel
score and the resulting PP3/BP4 code — so the aGVGD class this pipeline does NOT
compute is already baked into the published code. We extract the per-variant
"Preliminary bioinformatic code (missense only)" column (the protein-level call;
splice is handled separately by the pipeline's SpliceAI/OpenSpliceAI PP3 branch).

Output columns (tab-separated):
    hgvs_c   transcript change on NM_000546.6, e.g. "c.4G>C" (UNIQUE key)
    hgvs_p   protein change on NP_000537.3, e.g. "p.Glu2Gln" (fallback key)
    agvgd    Align-GVGD class, e.g. "Class C65" (evidence: why the code applies)
    bayesdel BayesDel score, e.g. "0.1490" (evidence: why the code applies)
    code     one of PP3 / PP3_moderate / BP4 / BP4_moderate / No evidence

The agvgd/bayesdel columns are NOT used for the call (the VCEP-assigned ``code``
already encodes them); they are carried through so a MET PP3/BP4 can explain *why*
in its evidence trail.

Usage:
    python scripts/build_tp53_codes.py --xlsx /path/to/PP3-BP4-codes.xlsx \
        --out resources/shared/tp53_pp3_bp4_codes.tsv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

# Column indices in the supplementary table (0-based); data starts at row 3.
_COL_C = 0       # Transcript change (NM_000546.6)
_COL_P = 1       # Protein change (NP_000537.3)
_COL_AGVGD = 2   # Align-GVGD class
_COL_BAYESDEL = 3  # BayesDel score
_COL_CODE = 4    # Preliminary bioinformatic code (missense only)
_DATA_START = 3

_VALID_CODES = {"PP3", "PP3_moderate", "BP4", "BP4_moderate", "No evidence"}


def build(xlsx_path: Path, out_path: Path) -> int:
    import openpyxl

    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    ws = wb.worksheets[0]

    rows: list[tuple[str, str, str, str, str]] = []
    seen_c: set[str] = set()
    for i, row in enumerate(ws.iter_rows(values_only=True)):
        if i < _DATA_START:
            continue
        c_change = row[_COL_C]
        p_change = row[_COL_P]
        code = row[_COL_CODE]
        if c_change is None or code is None:
            continue
        c = str(c_change).strip()
        p = str(p_change).strip() if p_change is not None else ""
        agvgd = str(row[_COL_AGVGD]).strip() if row[_COL_AGVGD] is not None else ""
        bd = row[_COL_BAYESDEL]
        bayesdel = f"{float(bd):.4f}" if isinstance(bd, (int, float)) else (str(bd).strip() if bd is not None else "")
        code = str(code).strip()
        if code not in _VALID_CODES:
            raise ValueError(f"Unexpected TP53 code {code!r} at row {i} ({c})")
        if c in seen_c:
            raise ValueError(f"Duplicate transcript change {c!r} — key must be unique")
        seen_c.add(c)
        rows.append((c, p, agvgd, bayesdel, code))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["hgvs_c", "hgvs_p", "agvgd", "bayesdel", "code"])
        w.writerows(rows)
    return len(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", type=Path, required=True,
                    help="ClinGen TP53 VCEP PP3/BP4 supplementary xlsx")
    ap.add_argument("--out", type=Path,
                    default=Path("resources/shared/tp53_pp3_bp4_codes.tsv"))
    args = ap.parse_args()
    n = build(args.xlsx, args.out)
    print(f"wrote {n} TP53 missense codes -> {args.out}")


if __name__ == "__main__":
    main()
