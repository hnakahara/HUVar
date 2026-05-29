from __future__ import annotations
from pathlib import Path
import pysam  # type: ignore

from acmg_classifier.models.variant import VariantRecord


def _normalize_chrom(fasta: pysam.FastaFile, chrom: str) -> str:
    """Match the input chromosome name to the form actually present in the
    FASTA index. FASTA files in the wild use either "chr1" (UCSC) or "1"
    (Ensembl). We try the input form first, then the toggled form, then
    give up — pysam will surface a clean error if neither exists."""
    if chrom in fasta.references:
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    return alt if alt in fasta.references else chrom


def _fetch_ref_base(fasta: pysam.FastaFile, chrom: str, pos: int) -> str:
    """Single base lookup at a 1-based position (pysam uses 0-based half-open).
    Used by left_align_and_trim while shifting indels upstream."""
    return fasta.fetch(chrom, pos - 1, pos).upper()


def left_align_and_trim(
    variant: VariantRecord,
    fasta_path: Path,
) -> VariantRecord:
    """
    Left-align and minimally represent an indel variant against the reference FASTA.

    SNVs pass through unchanged. Already left-normalised variants also pass through
    unchanged (idempotent). This is the standard VCF normalisation step performed
    before any database lookup so that coordinates match stored entries.

    Implements the vt-style normalisation algorithm (Tan, Abecasis &
    Kang 2015) directly because pysam/cyvcf2 don't expose a normalise
    primitive. Order is right-trim → left-shift → left-trim; all three
    steps are required for minimal-representation idempotency.
    """
    # SNVs are already in minimal representation by definition.
    if variant.is_snv:
        return variant

    fasta = pysam.FastaFile(str(fasta_path))
    try:
        chrom = _normalize_chrom(fasta, variant.chrom)
        pos = variant.pos
        ref = variant.ref
        alt = variant.alt

        # Right-trim shared suffix: e.g. ATC/AGC → AT/AG (the trailing C is
        # redundant). Stop when one side reaches length 1 — we must always
        # keep at least one base on each side for a valid VCF representation.
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]

        # Left-shift: if the last base of either allele equals the genomic
        # base immediately upstream, we can shift everything one position
        # left while keeping the same biological variant. This is what
        # makes the representation unique across equivalent representations
        # (e.g. CAA→CA at pos N is the same as AA→A at pos N-1).
        while len(ref) > 1 or len(alt) > 1:
            prev = _fetch_ref_base(fasta, chrom, pos - 1)
            if ref[-1] != prev and alt[-1] != prev:
                break
            ref = prev + ref[:-1] if len(ref) > 1 else prev
            alt = prev + alt[:-1] if len(alt) > 1 else prev
            pos -= 1

        # Left-trim shared prefix, leaving exactly one anchor base. Equivalent
        # to right-trim on the leading side; needed when the left-shift loop
        # adds extra anchor bases that aren't strictly necessary.
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            pos += 1

    finally:
        fasta.close()

    # Returning a model_copy (not mutating in place) keeps VariantRecord
    # effectively immutable for callers. Setting variant_id=None forces
    # model_post_init to regenerate the canonical key from the new coords.
    return variant.model_copy(
        update={
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_id": None,
        }
    )
