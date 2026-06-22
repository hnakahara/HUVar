#!/usr/bin/env python3
"""Build ps1_paralog_map.tsv from the SCN paralogue amino-acid alignment.

Input: the user-supplied "Paralogous-Genes_AminoAcid-Mapping.xlsx" — one row per
alignment column with the residue number of each paralogue (SCN1A/SCN2A/SCN3A/
SCN8A) in that column, or "NA" for a gap. The SCN epilepsy VCEPs (GN067-070)
grant PS1 from the same amino-acid change at the *analogous* residue of a
paralogue gene ("See Paralogous Gene Table").

Output: a TSV (one row per alignment column; columns = the paralogue genes;
empty = gap) consumed by
:class:`acmg_classifier.criteria.ps1_paralog.PS1ParalogMap`.
"""
from __future__ import annotations

import argparse
import csv
import re
import zipfile
from pathlib import Path


def _read_xlsx(path: Path) -> list[list[str]]:
    z = zipfile.ZipFile(path)
    ss = z.read("xl/sharedStrings.xml").decode("utf-8", "ignore")
    strings = [re.sub(r"<[^>]+>", "", m) for m in re.findall(r"<si>(.*?)</si>", ss, re.S)]
    sheet = z.read("xl/worksheets/sheet1.xml").decode("utf-8", "ignore")

    def colnum(ref: str) -> int:
        s = 0
        for ch in re.match(r"([A-Z]+)", ref).group(1):
            s = s * 26 + (ord(ch) - 64)
        return s

    rows: list[list[str]] = []
    for r in re.findall(r"<row[^>]*>(.*?)</row>", sheet, re.S):
        d: dict[int, str] = {}
        for c in re.findall(r"<c[^>]*/>|<c[^>]*>.*?</c>", r, re.S):
            ref = re.search(r'r="([A-Z]+\d+)"', c)
            if not ref:
                continue
            cn = colnum(ref.group(1))
            t = re.search(r't="([^"]+)"', c)
            v = re.search(r"<v>(.*?)</v>", c, re.S)
            if v is None:
                it = re.search(r"<t[^>]*>(.*?)</t>", c, re.S)
                d[cn] = it.group(1) if it else ""
            elif t and t.group(1) == "s":
                d[cn] = strings[int(v.group(1))]
            else:
                d[cn] = v.group(1)
        width = max(d) if d else 0
        rows.append([d.get(i, "") for i in range(1, width + 1)])
    return rows


def build(xlsx: Path, out: Path) -> int:
    rows = _read_xlsx(xlsx)
    header = [h.strip() for h in rows[0]]
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows[1:]:
            cells = [("" if (c.strip().upper() in ("", "NA")) else c.strip())
                     for c in (r + [""] * len(header))[:len(header)]]
            if any(cells):
                w.writerow(cells)
                n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True, help="Paralogous-Genes_AminoAcid-Mapping.xlsx")
    ap.add_argument("--out", default="resources/shared/ps1_paralog_map.tsv")
    args = ap.parse_args()
    n = build(Path(args.xlsx), Path(args.out))
    print(f"paralog alignment rows: {n} | written → {args.out}")


if __name__ == "__main__":
    main()
