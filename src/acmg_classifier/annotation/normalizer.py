from __future__ import annotations
from pathlib import Path
import pysam  # type: ignore

from acmg_classifier.models.variant import VariantRecord


def _normalize_chrom(fasta: pysam.FastaFile, chrom: str) -> str:
    """FASTA 内の実際の染色体名に合わせる (chr 有無どちらにも対応)。"""
    if chrom in fasta.references:
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    return alt if alt in fasta.references else chrom


def _fetch_ref_base(fasta: pysam.FastaFile, chrom: str, pos: int) -> str:
    """Return single reference base at 1-based position (0-based fetch)."""
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
    """
    if variant.is_snv:
        return variant

    fasta = pysam.FastaFile(str(fasta_path))
    try:
        chrom = _normalize_chrom(fasta, variant.chrom)
        pos = variant.pos
        ref = variant.ref
        alt = variant.alt

        # Right-trim shared suffix
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]

        # Left-shift while the last base of ref/alt equals the preceding genomic base
        while len(ref) > 1 or len(alt) > 1:
            prev = _fetch_ref_base(fasta, chrom, pos - 1)
            if ref[-1] != prev and alt[-1] != prev:
                break
            ref = prev + ref[:-1] if len(ref) > 1 else prev
            alt = prev + alt[:-1] if len(alt) > 1 else prev
            pos -= 1

        # Left-trim shared prefix (anchor base remains)
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            pos += 1

    finally:
        fasta.close()

    return variant.model_copy(
        update={
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_id": None,  # forces recalculation in model_post_init
        }
    )
