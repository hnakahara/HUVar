"""BLOSUM62 (Henikoff & Henikoff 1992) amino-acid substitution matrix.

A few ClinGen VCEPs gate PM5 on a BLOSUM62 comparison rather than Grantham
distance (e.g. PTEN GN003: the variant under evaluation "must have a BLOSUM62
score equal to or less than the known variant"). BLOSUM62 is a *similarity*
score — higher means more similar/conservative — so it runs in the opposite
direction to Grantham *distance*: a candidate that is chemically as-severe-or-
more than the comparator has a BLOSUM62 score **<=** the comparator's.

The published integer matrix is embedded verbatim (NCBI ``BLOSUM62``); only the
20 standard amino acids are encoded (the ambiguity codes B/Z/X and the stop ``*``
are not used by the PM5 gate, which compares concrete missense substitutions).

Reference: Henikoff S, Henikoff JG. "Amino acid substitution matrices from
protein blocks." PNAS. 1992;89(22):10915-10919.
"""
from __future__ import annotations

from acmg_classifier.criteria.grantham import AA3_TO_AA1

# Row/column order of the embedded matrix.
_ORDER = "ARNDCQEGHILKMFPSTWYV"
# Full symmetric BLOSUM62 (20x20), rows in _ORDER. Diagonal included.
_ROWS: tuple[tuple[int, ...], ...] = (
    ( 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0),  # A
    (-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3),  # R
    (-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3),  # N
    (-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3),  # D
    ( 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1),  # C
    (-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2),  # Q
    (-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2),  # E
    ( 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3),  # G
    (-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3),  # H
    (-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3),  # I
    (-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1),  # L
    (-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2),  # K
    (-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1),  # M
    (-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1),  # F
    (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2),  # P
    ( 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2),  # S
    ( 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0),  # T
    (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3),  # W
    (-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1),  # Y
    ( 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4),  # V
)


def _build_matrix() -> dict[tuple[str, str], int]:
    m: dict[tuple[str, str], int] = {}
    for i, row in enumerate(_ROWS):
        for j, score in enumerate(row):
            m[(_ORDER[i], _ORDER[j])] = score
    return m


_MATRIX = _build_matrix()


def _to_aa1(code: str | None) -> str | None:
    """Normalise a 1- or 3-letter amino-acid code to a single upper-case letter."""
    if not code:
        return None
    code = code.strip()
    if len(code) == 1:
        c = code.upper()
        return c if c in _ORDER else None
    if len(code) == 3:
        return AA3_TO_AA1.get(code.capitalize())
    return None


def blosum62_score(aa_a: str | None, aa_b: str | None) -> int | None:
    """BLOSUM62 substitution score between two residues, or ``None`` if either
    code is unknown. The matrix is symmetric; identical residues return the
    diagonal self-score."""
    a = _to_aa1(aa_a)
    b = _to_aa1(aa_b)
    if a is None or b is None:
        return None
    return _MATRIX.get((a, b))
