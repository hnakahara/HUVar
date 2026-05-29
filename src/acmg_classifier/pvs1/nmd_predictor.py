"""NMD (nonsense-mediated decay) prediction rules per ClinGen PVS1 2019."""
from __future__ import annotations
from acmg_classifier.models.annotation import ConsequenceInfo


def predicts_nmd(consequence: ConsequenceInfo) -> bool:
    """
    Return True if the variant is predicted to trigger NMD.

    NMD-escape rules (variant does NOT trigger NMD):
    1. Variant is in the last exon
    2. Variant is in the last 50 bp of the penultimate exon
       (i.e. within 50 bp upstream of the penultimate exon-exon junction)
    3. The exon affected is the only exon (single-exon gene/transcript)

    These rules are approximations from the ClinGen PVS1 2019 guidance.
    The function deliberately defaults to True ("NMD predicted") whenever
    exon information is missing or malformed — that is the conservative
    choice for a *pathogenicity* criterion because a false-positive PVS1
    is more recoverable than a false negative (the latter silently
    under-classifies a real LoF variant).
    """
    exon_str = consequence.exon  # VEP-style "5/12" (exon 5 of 12)
    if exon_str is None:
        return True

    parts = exon_str.split("/")
    if len(parts) != 2:
        return True

    try:
        exon_num = int(parts[0])
        total_exons = int(parts[1])
    except ValueError:
        return True

    # Single-exon transcripts have no downstream junction, so NMD cannot fire.
    if total_exons == 1:
        return False

    # Last exon: by definition no downstream junction, so the canonical NMD
    # surveillance mechanism does not engage.
    if exon_num == total_exons:
        return False

    # Penultimate exon: NMD escape applies only to the LAST 50 bp, which
    # cannot be determined from the exon index alone. We conservatively keep
    # the variant in the NMD-positive bucket — transcript_evaluator.py is
    # the coordinate-precise path when full VEP data is available.
    return True


def is_last_exon(consequence: ConsequenceInfo) -> bool:
    """Strict "exon == last exon" test parsed from VEP's "n/N" field.
    Used by the decision tree to decide whether to apply the last-exon
    Strong/Moderate fallback after NMD has been ruled out."""
    exon_str = consequence.exon
    if exon_str is None:
        return False
    parts = exon_str.split("/")
    if len(parts) != 2:
        return False
    try:
        return int(parts[0]) == int(parts[1])
    except ValueError:
        return False


def is_penultimate_exon(consequence: ConsequenceInfo) -> bool:
    """Penultimate-exon test (exon == last-1). Same conservative parsing as
    is_last_exon — invalid VEP data is treated as "not penultimate" so the
    decision tree falls through to the Supporting fallback rather than
    inflating strength on bad input."""
    exon_str = consequence.exon
    if exon_str is None:
        return False
    parts = exon_str.split("/")
    if len(parts) != 2:
        return False
    try:
        return int(parts[0]) == int(parts[1]) - 1
    except ValueError:
        return False
