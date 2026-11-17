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


# Paralogue groups without a machine-readable VCEP residue table: the analogous
# residue is taken from Needleman-Wunsch global alignments of the UniProt
# sequences (= MANE Select proteins) to the group's first ("hub") gene, and
# appended to the SCN table as one row per hub residue.
#   * KCNQ1 / KCNQ2 (GN112): PS1 paralogue approved for KCNQ2 only; validated
#     against the spec example KCNQ1 p.Thr144Ala <-> KCNQ2 p.Thr114Ala.
#   * RASopathy VCEP "Applicable for observed analogous residue positions in ...":
#     MAP2K1 / MAP2K2 (GN045; p.Tyr130 <-> p.Tyr134), SOS1 / SOS2 (GN041; p.Met269
#     <-> p.Met267), BRAF / RAF1 (GN049 / GN040; CR2 p.Ser365 <-> p.Ser259, kinase
#     p.Asp594 <-> p.Asp486) and HRAS / KRAS / NRAS / MRAS / RIT1 / RRAS2 (GN046 ...;
#     HRAS p.Gly12 / p.Gln61 <-> MRAS p.Gly22 / p.Gln71, RIT1 p.Gly30 / p.Gln79,
#     RRAS2 p.Gly23 / p.Gln72). KRAS is the MANE isoform 4B (P01116-2).
_GROUPS = (
    ("KCNQ1", "P51787", (("KCNQ2", "O43526"),)),
    ("MAP2K1", "Q02750", (("MAP2K2", "P36507"),)),
    ("SOS1", "Q07889", (("SOS2", "Q07890"),)),
    ("BRAF", "P15056", (("RAF1", "P04049"),)),
    ("HRAS", "P01112", (("KRAS", "P01116-2"), ("NRAS", "P01111"), ("MRAS", "O14807"),
                        ("RIT1", "Q92963"), ("RRAS2", "P62070"))),
)
_PAIR_GENES = [g for hub, _a, spokes in _GROUPS for g in (hub, *[sg for sg, _ in spokes])]


def _read_fasta(path: Path) -> dict[str, str]:
    """Sequences keyed by accession (">sp|ACC|..." or ">ACC")."""
    seqs: dict[str, str] = {}
    key = None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line.startswith(">"):
            head = line[1:].split()[0]
            key = head.split("|")[1] if "|" in head else head
            seqs[key] = ""
        elif key and line:
            seqs[key] += line
    return seqs


def _uniprot_seq(acc: str) -> str:
    import urllib.request
    text = urllib.request.urlopen(
        f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=60).read().decode()
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


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


def _pair_rows(scn: list[str], fasta: dict[str, str]) -> dict[str, list[list[str]]]:
    """Alignment rows per group (keyed by its hub gene), padded to the header:
    one row per hub residue with the aligned residue of every other member."""
    width = len(scn) + len(_PAIR_GENES)
    out: dict[str, list[list[str]]] = {}
    for hub, hub_acc, spokes in _GROUPS:
        try:
            hub_seq = fasta.get(hub_acc) or _uniprot_seq(hub_acc)
            maps = {g: _nw_align(hub_seq, fasta.get(acc) or _uniprot_seq(acc))
                    for g, acc in spokes}
        except Exception as exc:  # offline / UniProt down
            print(f"  [warn] {hub} group alignment skipped ({exc})")
            continue
        rows = []
        for p in range(1, len(hub_seq) + 1):
            hits = {g: m[p] for g, m in maps.items() if p in m}
            if not hits:
                continue
            row = [""] * width
            row[len(scn) + _PAIR_GENES.index(hub)] = str(p)
            for g, q in hits.items():
                row[len(scn) + _PAIR_GENES.index(g)] = str(q)
            rows.append(row)
        out[hub] = rows
    return out


def _existing(out: Path) -> tuple[list[str], list[list[str]]]:
    with out.open(encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]


def build(xlsx: Path | None, out: Path, fasta: Path | None = None) -> int:
    """Write the table. Without *xlsx* the SCN rows (and any pair that cannot be
    re-aligned) are kept from the existing *out* file."""
    old_header, old_rows = _existing(out) if out.exists() else ([], [])
    if xlsx is not None:
        rows = _read_xlsx(xlsx)
        scn = [h.strip() for h in rows[0]]
        scn_rows = []
        for r in rows[1:]:
            cells = [("" if (c.strip().upper() in ("", "NA")) else c.strip())
                     for c in (r + [""] * len(scn))[:len(scn)]]
            if any(cells):
                scn_rows.append(cells)
    else:
        scn = [h for h in old_header if h not in _PAIR_GENES]
        scn_rows = [r[:len(scn)] for r in old_rows if any(r[:len(scn)])]
    pairs = _pair_rows(scn, _read_fasta(fasta) if fasta else {})
    header = scn + _PAIR_GENES
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for cells in scn_rows:
            w.writerow(cells + [""] * len(_PAIR_GENES))
            n += 1
        for hub, _acc, spokes in _GROUPS:
            rows = pairs.get(hub)
            if rows is None and hub in old_header:   # keep the existing alignment
                members = [g for g in (hub, *[sg for sg, _ in spokes]) if g in old_header]
                idx = {g: old_header.index(g) for g in members}
                rows = []
                for r in old_rows:
                    r = r + [""] * (len(old_header) - len(r))
                    if r[idx[hub]] and any(r[idx[g]] for g in members if g != hub):
                        row = [""] * len(header)
                        for g in members:
                            row[header.index(g)] = r[idx[g]]
                        rows.append(row)
            for row in rows or []:
                w.writerow(row)
                n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", help="Paralogous-Genes_AminoAcid-Mapping.xlsx (SCN); "
                    "omit to keep the SCN rows of the existing table")
    ap.add_argument("--fasta", help="UniProt FASTA with the pair sequences (offline use); "
                    "missing accessions are fetched from UniProt")
    ap.add_argument("--out", default="resources/shared/ps1_paralog_map.tsv")
    args = ap.parse_args()
    n = build(Path(args.xlsx) if args.xlsx else None, Path(args.out),
              Path(args.fasta) if args.fasta else None)
    print(f"paralog alignment rows: {n} | written → {args.out}")


if __name__ == "__main__":
    main()
