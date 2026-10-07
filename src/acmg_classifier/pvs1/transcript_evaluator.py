"""Evaluate transcript- and gene-level properties needed for the PVS1 decision tree."""
from __future__ import annotations
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo


def has_alternative_transcript_rescue(annotation: AnnotationData) -> bool:
    """
    Return True if an alternative MANE/canonical transcript does NOT carry the LoF consequence.

    When a variant causes LoF on one transcript but is tolerated on another clinically
    relevant transcript, PVS1 strength is reduced (e.g., VeryStrong -> Strong).
    """
    primary = annotation.primary_consequence
    if primary is None:
        return False

    from acmg_classifier.models.enums import ConsequenceType
    lof_consequences = {
        ConsequenceType.FRAMESHIFT,
        ConsequenceType.STOP_GAINED,
        ConsequenceType.SPLICE_ACCEPTOR,
        ConsequenceType.SPLICE_DONOR,
        ConsequenceType.START_LOST,
        ConsequenceType.TRANSCRIPT_ABLATION,
    }

    # A "rescue" exists only when BOTH conditions hold: at least one transcript
    # carries the LoF, AND at least one clinically-relevant (MANE/canonical)
    # transcript carries a non-LoF consequence. Minor isoforms are ignored.
    lof_transcripts = 0
    non_lof_mane_canonical = 0

    for c in annotation.consequences:
        if c.consequence in lof_consequences:
            lof_transcripts += 1
        elif c.is_mane_select or c.is_canonical:
            non_lof_mane_canonical += 1

    return non_lof_mane_canonical > 0 and lof_transcripts > 0


# Minimum number of P/LP null variants (nonsense, frameshift, canonical ±1/2
# splice; ClinVar review status >=1 star) that establishes LoF as a disease
# mechanism for a gene without a VCEP PVS1 decision or a ClinGen HI score of 3.
# The use of reported pathogenic null variants follows the approach of the
# Franklin platform; the value 3 is our choice of the minimum indicating
# recurrent pathogenic null variants (not independently calibrated).
_MIN_PLP_NULL = 3


def gene_has_lof_mechanism(
    consequence: ConsequenceInfo | None,
    gnomad_loeuf: float | None = None,
    clinvar_plp_null: int = 0,
    hi_score: str | None = None,
) -> bool:
    """Is loss-of-function an established disease mechanism for the gene?

    Established when the gene has a ClinGen Dosage Sensitivity haploinsufficiency
    score of 3, or when ClinVar reports >= 3 P/LP null variants for the gene.
    (A VCEP that explicitly applies PVS1 is handled by the caller.)

    gnomAD LOEUF is NOT used: population constraint alone does not establish
    LoF as a disease mechanism (ACGS 2023; AutoPVS1 uses pLI for reference only).
    The argument is retained for backward compatibility and reporting.
    """
    return lof_mechanism_reason(gnomad_loeuf, clinvar_plp_null, hi_score)[0]


def lof_mechanism_reason(
    gnomad_loeuf: float | None,
    clinvar_plp_null: int,
    hi_score: str | None,
) -> tuple[bool, str]:
    """``(established, explanation)`` — see :func:`gene_has_lof_mechanism`."""
    loeuf_note = f"; LOEUF={gnomad_loeuf:.2f} (reference only)" if gnomad_loeuf is not None else ""
    if (hi_score or "").strip() == "3":
        return True, f"ClinGen haploinsufficiency score 3{loeuf_note}"
    if clinvar_plp_null >= _MIN_PLP_NULL:
        return True, f"{clinvar_plp_null} ClinVar P/LP null variants (>= {_MIN_PLP_NULL}){loeuf_note}"
    hi_note = f"ClinGen HI score {hi_score}" if hi_score else "no ClinGen HI score"
    return False, (
        f"LoF mechanism not established: {hi_note}, {clinvar_plp_null} ClinVar P/LP "
        f"null variants (< {_MIN_PLP_NULL}){loeuf_note}"
    )
