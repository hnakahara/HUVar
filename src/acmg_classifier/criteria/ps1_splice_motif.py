"""BMPR2 VCEP (GN125 v2.0.0) PS1 for splice variants — splice-motif rule.

  "PS1_Moderate: For splice variants outside splice donor/acceptor ±1,2
  dinucleotide positions but residing in the same splice donor motif (last 3
  bases of the exon and intronic positions +3 to +6) or acceptor motif (first
  base of the exon and intronic positions -20 to -3) as a previously established
  pathogenic variant." — PS1_Supporting with a likely pathogenic variant.

The variant and the comparator must lie in the SAME motif (the same exon end for
a donor motif, the same exon start for an acceptor motif); the canonical ±1,2
positions are not part of the motif. Exonic positions are identified from the
MANE Select exon structure. Exonic variants qualify only when they do not change
the amino acid (synonymous), and exonic comparators only when they do not change
the protein, so that the comparator acts through splicing.
"""
from __future__ import annotations

import re
from typing import Optional

MOTIF_GENES: frozenset[str] = frozenset({"BMPR2"})

_INTRONIC = re.compile(r"c\.(\d+)([+-])(\d+)(?![0-9])")
_EXONIC = re.compile(r"c\.(\d+)(?:_\d+)?[ACGT]")


def _bounds(transcript_id: str) -> Optional[tuple[set[int], set[int]]]:
    from acmg_classifier.pvs1.last_nt import _exon_bounds
    b = _exon_bounds(transcript_id)
    if not b:
        return None
    ends = {e for _, e in b}
    starts = {e + 1 for _, e in b}
    return ends, starts


def motif_of(hgvs_c: Optional[str], transcript_id: str) -> Optional[tuple[str, int]]:
    """``("D", exon_end)`` / ``("A", exon_start)`` when the c. position lies in the
    BMPR2 donor / acceptor motif (excluding ±1,2), else ``None``."""
    c = (hgvs_c or "").split(":")[-1]
    m = _INTRONIC.match(c)
    if m:
        base, sign, off = int(m.group(1)), m.group(2), int(m.group(3))
        if sign == "+" and 3 <= off <= 6:
            return "D", base
        if sign == "-" and 3 <= off <= 20:
            return "A", base
        return None
    m = _EXONIC.match(c)
    if not m:
        return None
    x = int(m.group(1))
    bd = _bounds(transcript_id)
    if bd is None:
        return None
    ends, starts = bd
    for n in ends:
        if n - 2 <= x <= n:
            return "D", n
    if x in starts:
        return "A", x
    return None


def same_transcript(hgvs_c: Optional[str], transcript_id: str) -> bool:
    """A ClinVar HGVS c. name is comparable when it has no accession or the same
    RefSeq accession (version ignored) as the MANE transcript."""
    acc = (hgvs_c or "").split(":")[0] if ":" in (hgvs_c or "") else ""
    if not acc:
        return True
    return acc.split(".")[0] == (transcript_id or "").split(".")[0]
