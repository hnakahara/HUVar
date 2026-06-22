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


# KCNQ1 PS1 paralogue (GN112) is approved for KCNQ2 only; the paralogue table the
# VCEP references (cardiodb) gives no KCNQ2 residue numbers, so the analogous
# residue is taken from a Needleman-Wunsch global alignment of the two UniProt
# sequences (validated against the spec example KCNQ1 p.Thr144Ala ↔ KCNQ2
# p.Thr114Ala). The alignment is appended to the SCN table as KCNQ1/KCNQ2 rows.
_KCNQ1_ACC, _KCNQ2_ACC = "P51787", "O43526"


def _uniprot_seq(acc: str) -> str:
    import json
    import urllib.request
    j = json.load(urllib.request.urlopen(
        f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60))
    return j["sequence"]["value"]


def _nw_align(a: str, b: str, gap: int = -8) -> dict[int, int]:
    """Needleman-Wunsch global alignment (BLOSUM62). Returns {a_pos: b_pos} for
    aligned (non-gap) columns, 1-based."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from acmg_classifier.criteria.blosum62 import blosum62_score
    n, m = len(a), len(b)
    prev = [gap * j for j in range(m + 1)]
    tb = bytearray((n + 1) * (m + 1))
    for j in range(1, m + 1):
        tb[j] = 2
    cur = [0] * (m + 1)
    for i in range(1, n + 1):
        cur[0] = gap * i
        tb[i * (m + 1)] = 1
        ai = a[i - 1]
        for j in range(1, m + 1):
            s = blosum62_score(ai, b[j - 1])
            s = -4 if s is None else s
            diag = prev[j - 1] + s
            up = prev[j] + gap
            left = cur[j - 1] + gap
            best, d = diag, 0
            if up > best:
                best, d = up, 1
            if left > best:
                best, d = left, 2
            cur[j] = best
            tb[i * (m + 1) + j] = d
        prev, cur = cur, prev
    out: dict[int, int] = {}
    i, j = n, m
    while i > 0 and j > 0:
        d = tb[i * (m + 1) + j]
        if d == 0:
            out[i] = j
            i -= 1
            j -= 1
        elif d == 1:
            i -= 1
        else:
            j -= 1
    return out


def _kcnq_rows(ncols: int) -> list[list[str]]:
    """KCNQ1↔KCNQ2 alignment rows (padded to the SCN column width + KCNQ1/KCNQ2)."""
    try:
        q1, q2 = _uniprot_seq(_KCNQ1_ACC), _uniprot_seq(_KCNQ2_ACC)
    except Exception as exc:  # offline / UniProt down → SCN-only table
        print(f"  [warn] KCNQ1/KCNQ2 alignment skipped ({exc})")
        return []
    amap = _nw_align(q1, q2)
    rows = []
    for p1 in sorted(amap):
        rows.append([""] * ncols + [str(p1), str(amap[p1])])
    return rows


def build(xlsx: Path, out: Path) -> int:
    rows = _read_xlsx(xlsx)
    scn_header = [h.strip() for h in rows[0]]
    header = scn_header + ["KCNQ1", "KCNQ2"]
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows[1:]:
            cells = [("" if (c.strip().upper() in ("", "NA")) else c.strip())
                     for c in (r + [""] * len(scn_header))[:len(scn_header)]]
            if any(cells):
                w.writerow(cells + ["", ""])  # SCN row, KCNQ cols empty
                n += 1
        for kr in _kcnq_rows(len(scn_header)):
            w.writerow(kr)
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
