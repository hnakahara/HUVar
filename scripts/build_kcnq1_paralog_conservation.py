"""Build resources/shared/kcnq1_paralog_conservation.tsv: for every KCNQ1 residue
(NM_000218.3 / NP_000209.2, ENST00000155840) the residue each human KCNQ paralogue
(KCNQ2-5) has at the aligned position of the cardiodb paralogue-annotation
alignment (T-COFFEE / M-Coffee, 32 Kv channels).

The LQTS VCEP (KCNQ1, GN112) applies PM5 only at a residue "highly conserved across
all 5 human KCNQ paralogues" and checks it on
https://www.cardiodb.org/paralogue_annotation/gene.php?name=KCNQ1; "poor
conservation (ineligibility for PM5) is defined whenever 1 or more KCNQ's show a
different amino acid at the position". A gap counts as a different residue.

usage: python scripts/build_kcnq1_paralog_conservation.py ALIGNMENT OUT.tsv
ALIGNMENT: https://www.cardiodb.org/paralogue_annotation/images/KCNQ1_alignment.htm
           saved as HTML or as page text.
"""
from __future__ import annotations

import html
import re
import sys

GENES = ("KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5")
_LINE = re.compile(r"^(KCNQ[1-5])\s+\d+\s+(\S+)\s+\d+\s*$")


def read_alignment(path: str) -> dict[str, str]:
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = html.unescape(re.sub(r"<[^>]+>", "", fh.read()))
    seqs: dict[str, str] = {}
    for line in text.splitlines():
        m = _LINE.match(line.strip())
        if m:
            seqs[m.group(1)] = seqs.get(m.group(1), "") + m.group(2)
    missing = [g for g in GENES if g not in seqs]
    if missing:
        raise SystemExit(f"{path}: no alignment rows for {', '.join(missing)}")
    if len({len(seqs[g]) for g in GENES}) != 1:
        raise SystemExit(f"{path}: alignment rows of unequal length")
    return seqs


def main() -> None:
    seqs = read_alignment(sys.argv[1])
    q1 = seqs["KCNQ1"]
    if len(q1.replace("-", "")) != 676:
        raise SystemExit("KCNQ1 row is not the 676-residue NP_000209.2")
    with open(sys.argv[2], "w", encoding="utf-8", newline="\n") as fh:
        fh.write("# KCNQ1 (NP_000209.2) residues and the aligned KCNQ2-5 residues of the cardiodb\n"
                 "# paralogue-annotation alignment (KCNQ1_alignment.htm, T-COFFEE). conserved=1\n"
                 "# when all four paralogues have the KCNQ1 residue (LQTS VCEP PM5 eligibility).\n"
                 "# Built by scripts/build_kcnq1_paralog_conservation.py\n")
        fh.write("position\tkcnq1\tkcnq2\tkcnq3\tkcnq4\tkcnq5\tconserved\n")
        pos = 0
        for i, aa in enumerate(q1):
            if aa == "-":
                continue
            pos += 1
            other = [seqs[g][i] for g in GENES[1:]]
            cons = int(all(o == aa for o in other))
            fh.write(f"{pos}\t{aa}\t" + "\t".join(other) + f"\t{cons}\n")


if __name__ == "__main__":
    main()
