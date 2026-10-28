#!/usr/bin/env python3
"""Convert the ClinGen ENIGMA BRCA1/BRCA2 VCEP Specifications Table 9 (xlsx) into
``resources/shared/enigma_brca12_table9.tsv``.

Table 9 lists, per variant, the PS3 or BS3 code (and strength) the VCEP assigns
from calibrated functional assays, after its own handling of splicing (assays
measuring protein only are not used for variants predicted/observed to alter
splicing). Variants coded ``None`` get neither code. HGVS c. is on the MANE
Select transcripts (BRCA1 NM_007294.4, BRCA2 NM_000059.4).

Usage (needs openpyxl):
  python scripts/build_enigma_table9.py --xlsx Specifications_Table9_V1.2_2024-11-18.xlsx
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

_COLS = ["Gene", "HGVS Nucleotide", "HGVS Protein", "Assigned Code", "Code Weight",
         "Standardised Text", "Predicted or Observed Splicing", "# Pubs",
         "Result1", "Result2", "Result3", "Result4"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--xlsx", type=Path, required=True)
    ap.add_argument("--version", default="1.2 (2024-11-18)")
    ap.add_argument("--out", type=Path, default=Path("resources/shared/enigma_brca12_table9.tsv"))
    a = ap.parse_args()
    import openpyxl

    ws = openpyxl.load_workbook(a.xlsx, read_only=True, data_only=True).worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    hdr_i = next(i for i, r in enumerate(rows) if r and r[0] == "Gene")
    hdr = [str(c) if c is not None else "" for c in rows[hdr_i]]
    idx = {c: hdr.index(c) for c in _COLS}
    n = 0
    with a.out.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# ClinGen ENIGMA BRCA1/BRCA2 VCEP Specifications Table 9, version {a.version} "
                 "(GN092 / GN097 v1.2.0). HGVS c. on NM_007294.4 (BRCA1) / NM_000059.4 (BRCA2).\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "hgvs_c", "hgvs_p", "code", "strength", "splicing", "n_pubs",
                    "results", "text"])
        for r in rows[hdr_i + 2:]:   # skip the column-description row
            if not r or r[0] not in ("BRCA1", "BRCA2"):
                continue
            def g(c: str, r=r) -> str:
                return "" if r[idx[c]] is None else re.sub(r"\s+", " ", str(r[idx[c]])).strip()
            results = " | ".join(x for x in (g("Result1"), g("Result2"), g("Result3"), g("Result4")) if x)
            w.writerow([g("Gene"), g("HGVS Nucleotide"), g("HGVS Protein"), g("Assigned Code"),
                        g("Code Weight"), g("Predicted or Observed Splicing"), g("# Pubs"),
                        results, g("Standardised Text")])
            n += 1
    print(f"{n} variants -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
