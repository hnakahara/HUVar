"""Grantham (1974) amino-acid distance matrix.

Several ClinGen VCEPs gate PM5 on a Grantham-distance comparison: the variant
under evaluation must be chemically *as different or more different* from the
wild-type residue than a known pathogenic/likely-pathogenic change at the same
codon (e.g. PIK3CD, PIK3R1, VHL, HNF1A — see ``pm5_genes`` / the ``pm5_grantham``
column in ``disease_prevalence.tsv``). The VCEP texts cite "Grantham, 1974,
Table 2", so the published integer matrix is embedded verbatim here rather than
recomputed from the (c, p, v) formula, whose rounding differs from the table by
±1 at a few pairs (e.g. Cys-Trp: formula 214 vs table 215).

Reference: Grantham R. "Amino acid difference formula to help explain protein
evolution." Science. 1974;185(4154):862-864. Table 2.
"""
from __future__ import annotations

# Three-letter (HGVS) → one-letter amino-acid codes. Used to parse the
# comparator's protein change (ClinVar normalises to 3-letter: p.Arg156His).
AA3_TO_AA1: dict[str, str] = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# Residue order of Grantham 1974 Table 2 and the lower-triangle distances.
# _TRIANGLE[i] lists distance(_ORDER[i], _ORDER[j]) for every j > i.
_ORDER = "SRLPTAVGIFYCHQNKDEMW"
_TRIANGLE: tuple[tuple[int, ...], ...] = (
    (110, 145, 74, 58, 99, 124, 56, 142, 155, 144, 112, 89, 68, 46, 121, 65, 80, 135, 177),
    (102, 103, 71, 112, 96, 125, 97, 97, 77, 180, 29, 43, 86, 26, 96, 54, 91, 101),
    (98, 92, 96, 32, 138, 5, 22, 36, 198, 99, 113, 153, 107, 172, 138, 15, 61),
    (38, 27, 68, 42, 95, 114, 110, 169, 77, 76, 91, 103, 108, 93, 87, 147),
    (58, 69, 59, 89, 103, 92, 149, 47, 42, 65, 78, 85, 65, 81, 128),
    (64, 60, 94, 113, 112, 195, 86, 91, 111, 106, 126, 107, 84, 148),
    (109, 29, 50, 55, 192, 84, 96, 133, 97, 152, 121, 21, 88),
    (135, 153, 147, 159, 98, 87, 80, 127, 94, 98, 127, 184),
    (21, 33, 198, 94, 109, 149, 102, 168, 134, 10, 61),
    (22, 205, 100, 116, 158, 102, 177, 140, 28, 40),
    (194, 83, 99, 143, 85, 160, 122, 36, 37),
    (174, 154, 139, 202, 154, 170, 196, 215),
    (24, 68, 32, 81, 40, 87, 115),
    (46, 53, 61, 29, 101, 130),
    (94, 23, 42, 142, 174),
    (101, 56, 95, 110),
    (45, 160, 181),
    (126, 152),
    (126,),
)


def _build_matrix() -> dict[frozenset[str], int]:
    m: dict[frozenset[str], int] = {}
    for i, row in enumerate(_TRIANGLE):
        for off, dist in enumerate(row):
            j = i + 1 + off
            m[frozenset((_ORDER[i], _ORDER[j]))] = dist
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


def grantham_distance(aa_a: str | None, aa_b: str | None) -> int | None:
    """Grantham (1974) distance between two residues, or ``None`` if either code
    is unknown. Identical residues return 0; the matrix is symmetric."""
    a = _to_aa1(aa_a)
    b = _to_aa1(aa_b)
    if a is None or b is None:
        return None
    if a == b:
        return 0
    return _MATRIX.get(frozenset((a, b)))
