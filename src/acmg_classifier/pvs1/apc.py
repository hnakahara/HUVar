"""APC-specific PVS1 decision tree (InSiGHT / Tayoun et al. 2018, 2023 update).

APC's PVS1 deviates from the generic ClinGen tree in two modelled ways:

* **Nonsense / frameshift codon-range gate** — PVS1 applies only when the
  premature truncation / frameshift lies within codons **49-2645** of
  NM_000038.6; upstream of 49 or downstream of 2645 it is N/A (an N-terminal
  alternative start or a non-critical C-terminus).
* **Canonical GT-AG +/-1,2 splice variants and "G to non-G last nucleotide"
  exonic changes** — an explicit, ALLELE-SPECIFIC strength table (Lists A-E),
  reflecting that e.g. ``c.835-1G>A`` is PVS1 but ``c.835-1G>C/T`` is only
  Moderate, and the "G to non-G last NT" exonic changes are downgraded.

Out of scope (fall back to the generic tree / not modelled): the exon-level
deletion & duplication rules (Fig 1B), and the generic exon-skip splice logic
for splice variants not in the lists.
"""
from __future__ import annotations

from acmg_classifier.models.enums import ConsequenceType, CriterionStrength

# Codon range (NM_000038.6) within which a truncating variant is PVS1.
_APC_TRUNC_MIN = 49
_APC_TRUNC_MAX = 2645

# --- Lists A-E (Fig 1B). Each raw entry "c.<pos><REF>><ALT1>,<ALT2>,..."
# expands to one cDNA key per alt allele (the assignments are allele-specific).
_LIST_A = [  # PVS1
    "c.136-1G>A,C,T", "c.136-2A>C,G,T", "c.220+1G>A,C,T", "c.220+2T>A,C,G",
    "c.221-1G>A,C,T", "c.221-2A>C,G,T", "c.422+1G>A,C,T", "c.422+2T>A,C,G",
    "c.423-1G>A,C,T", "c.423-2A>C,G,T", "c.531+1G>A,C,T", "c.531+2T>A,C,G",
    "c.532-1G>A,C,T", "c.532-2A>C,G,T",
    "c.646-1G>A,C,T", "c.646-2A>C,G,T", "c.730-1G>A,C,T", "c.834+1G>A,C,T",
    "c.834+2T>A,C,G", "c.835-1G>A", "c.933+1G>A,C,T", "c.933+2T>A,C,G",
    "c.1312+1G>A,C,T", "c.1312+2T>A,C,G", "c.1409-1G>A,C,T", "c.1409-2A>C,G,T",
    "c.1548+1G>A,C,T", "c.1548+2T>A,G",
    "c.1549-1G>A,C,T", "c.1549-2A>C,G,T", "c.1626+1G>A,C,T", "c.1626+2T>A,C,G",
    "c.1627-1G>A,C,T", "c.1627-2A>C,G,T", "c.1743+1G>A,C,T", "c.1743+2T>A,C,G",
    "c.1744-1G>A,C,T", "c.1744-2A>C,G,T", "c.1958+1G>A,C,T", "c.1958+2T>A,C,G",
    "c.1959-1G>A",
]
_LIST_B = [  # PVS1_Strong
    "c.220G>A,C,T", "c.422G>A,C,T", "c.834G>A,C,T", "c.1548G>A,C,T",
    "c.1548+2T>C", "c.1626G>A,C,T", "c.1743G>A,C,T", "c.1958G>A,C,T",
]
_LIST_C = [  # PVS1_Moderate
    "c.645+1G>A,C,T", "c.645+2T>A,G", "c.729+1G>A,C,T", "c.729+2T>A,G",
    "c.730-2A>C,G,T", "c.835-1G>C,T", "c.835-2A>C,G,T", "c.1408+1G>A,C,T",
    "c.1408+2T>A,C,G",
]
_LIST_D = [  # PVS1_Supporting
    "c.729+2T>C", "c.933G>A,C,T",
]
_LIST_E = [  # N/A
    "c.-18-1G>A,C,T", "c.-18-2A>C,G,T", "c.135G>A,C,T", "c.135+1G>A,C,T",
    "c.135+2T>A,C,G", "c.645G>A,T,C", "c.645+2T>C", "c.729G>A,T,C",
    "c.934-1G>A,C,T", "c.934-2A>C,G,T", "c.1313-1G>A,C,T", "c.1313-2A>C,G,T",
    "c.1408G>A,C,T", "c.1959-1G>C,T", "c.1959-2A>C,G,T",
]


def _expand(raw: str) -> list[str]:
    """"c.835-1G>A,C,T" -> ["c.835-1G>A", "c.835-1G>C", "c.835-1G>T"]."""
    prefix, _, alts = raw.partition(">")
    return [f"{prefix}>{a}" for a in alts.split(",")]


def _build_table() -> dict[str, CriterionStrength]:
    table: dict[str, CriterionStrength] = {}
    for entries, strength in (
        (_LIST_A, CriterionStrength.VERY_STRONG),
        (_LIST_B, CriterionStrength.STRONG),
        (_LIST_C, CriterionStrength.MODERATE),
        (_LIST_D, CriterionStrength.SUPPORTING),
        (_LIST_E, CriterionStrength.NOT_MET),
    ):
        for raw in entries:
            for key in _expand(raw):
                table[key] = strength
    return table


_APC_LIST = _build_table()
_STRENGTH_LABEL = {
    CriterionStrength.VERY_STRONG: "PVS1 (List A)",
    CriterionStrength.STRONG: "PVS1_Strong (List B)",
    CriterionStrength.MODERATE: "PVS1_Moderate (List C)",
    CriterionStrength.SUPPORTING: "PVS1_Supporting (List D)",
    CriterionStrength.NOT_MET: "N/A (List E)",
}


def _cdna_key(hgvs_c: str | None) -> str:
    """Bare cDNA change ('c.422G>C') from a (possibly transcript-prefixed) HGVS
    coding string, for matching against the APC lists."""
    if not hgvs_c:
        return ""
    return hgvs_c.split(":")[-1].strip()


def evaluate_apc_pvs1(pc) -> tuple[CriterionStrength, str] | None:
    """APC-specific PVS1 strength, or None when APC's special rules do not
    resolve the variant (the caller then falls back to the generic tree).

    Returns ``CriterionStrength.NOT_MET`` for an explicit N/A (List E or a
    truncation outside codons 49-2645)."""
    # 1. Explicit allele-specific list (canonical splice / last-NT exonic). This
    #    runs first because it can fire on changes VEP calls missense/synonymous.
    key = _cdna_key(pc.hgvs_c)
    if key and key in _APC_LIST:
        strength = _APC_LIST[key]
        return strength, f"APC VCEP {_STRENGTH_LABEL[strength]}: {key}"

    # 2. Nonsense / frameshift codon-range gate (NM_000038.6 codons 49-2645).
    if pc.consequence in (ConsequenceType.STOP_GAINED, ConsequenceType.FRAMESHIFT):
        codon = pc.protein_position
        if codon is None:
            return None  # cannot place the truncation → defer to generic tree
        if _APC_TRUNC_MIN <= codon <= _APC_TRUNC_MAX:
            return (
                CriterionStrength.VERY_STRONG,
                f"APC: truncation at codon {codon} within 49-2645 -> PVS1",
            )
        return (
            CriterionStrength.NOT_MET,
            f"APC: truncation at codon {codon} outside 49-2645 -> PVS1 N/A",
        )

    return None
