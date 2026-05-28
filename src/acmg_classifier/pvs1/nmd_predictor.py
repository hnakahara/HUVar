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
    """
    exon_str = consequence.exon  # e.g. "5/12"
    if exon_str is None:
        return True  # assume NMD if no exon info

    parts = exon_str.split("/")
    if len(parts) != 2:
        return True

    try:
        exon_num = int(parts[0])
        total_exons = int(parts[1])
    except ValueError:
        return True

    # Single-exon transcript
    if total_exons == 1:
        return False

    # Last exon
    if exon_num == total_exons:
        return False

    # Penultimate exon, last 50 bp — approximated by checking if exon_num == total_exons - 1
    # Full implementation would check the genomic coordinates against the exon boundary.
    # Without exact CDS offset we conservatively allow NMD for penultimate exon.
    # (transcript_evaluator.py provides coordinate-precise check when VEP data is complete)

    return True


def is_last_exon(consequence: ConsequenceInfo) -> bool:
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
