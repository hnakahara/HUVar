"""Variant normalisation: minimal representation + left alignment.

VCFs from GATK and similar callers often carry indels with extra padding bases,
e.g. ``GAAA>GAA`` for a single-A deletion (common when a multi-allelic record such
as ``GAAA>GAA,G`` is split). ClinVar, gnomAD and the SpliceAI score files store the
normalised form (``GA>G``, left-aligned), and every lookup here is an exact
CHROM:POS:REF:ALT match, so an un-normalised indel silently misses all of them
(absent from gnomAD -> PM2 met, no ClinVar match, no splice score).

``normalize_alleles`` implements the vt algorithm (Tan, Abecasis & Kang 2015):
repeatedly right-trim the shared last base and, when an allele becomes empty,
extend both alleles one base to the left from the reference; then trim the
shared leading bases down to a single anchor. Without a reference only the
trimming is done (enough for padding, but repeats are not left-shifted).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from acmg_classifier.models.variant import VariantRecord


def _contig(fasta, chrom: str) -> Optional[str]:
    """The FASTA's name for *chrom* ("chr1" / "1" / "chrM" / "MT"), or None."""
    names = fasta.references
    bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
    candidates = [chrom, bare, "chr" + bare]
    if bare.upper() in ("M", "MT"):
        candidates += ["chrM", "MT", "M"]
    for c in candidates:
        if c in names:
            return c
    return None


def normalize_alleles(chrom: str, pos: int, ref: str, alt: str,
                      fasta=None) -> tuple[int, str, str]:
    """Normalised ``(pos, ref, alt)``. ``fasta`` is an open ``pysam.FastaFile`` or None.

    SNVs and identical alleles are returned unchanged. Left alignment needs the
    reference and is only attempted when the REF allele matches it at ``pos``
    (a mismatching REF -- e.g. a VCF of another assembly -- is only trimmed)."""
    ref, alt = ref.upper(), alt.upper()
    if (len(ref) == 1 and len(alt) == 1) or ref == alt:
        return pos, ref, alt
    original = (pos, ref, alt)

    contig = None
    if fasta is not None:
        contig = _contig(fasta, chrom)
        if contig is not None:
            try:
                genome = fasta.fetch(contig, pos - 1, pos - 1 + len(ref)).upper()
            except (ValueError, KeyError, IndexError):
                genome = ""
            if genome != ref:
                contig = None

    if contig is None:
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
    else:
        while True:
            changed = False
            if ref and alt and ref[-1] == alt[-1]:
                ref, alt = ref[:-1], alt[:-1]
                changed = True
            if (not ref or not alt) and pos > 1:
                base = fasta.fetch(contig, pos - 2, pos - 1).upper()
                ref, alt = base + ref, base + alt
                pos -= 1
                changed = True
            if not changed:
                break
        if not ref or not alt:          # at the start of the contig: leave it as it was
            return original

    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def open_fasta(fasta_path: Optional[Path]):
    """An open ``pysam.FastaFile`` for *fasta_path*, or None when unavailable."""
    if fasta_path is None or not Path(fasta_path).exists():
        return None
    try:
        import pysam  # type: ignore
        return pysam.FastaFile(str(fasta_path))
    except Exception:  # noqa: BLE001 -- normalisation falls back to trimming
        return None


def normalize_record(variant: VariantRecord, fasta=None) -> VariantRecord:
    """*variant* in normalised form; ``input_variant`` keeps the original key when
    the representation changed."""
    pos, ref, alt = normalize_alleles(variant.chrom, variant.pos, variant.ref,
                                      variant.alt, fasta)
    if (pos, ref, alt) == (variant.pos, variant.ref, variant.alt):
        return variant
    # model_copy skips model_post_init, so set the derived key explicitly
    return variant.model_copy(update={
        "pos": pos, "ref": ref, "alt": alt,
        "variant_id": f"{variant.chrom}:{pos}:{ref}:{alt}",
        "input_variant": variant.input_variant or variant.key,
    })


def left_align_and_trim(variant: VariantRecord, fasta_path: Path) -> VariantRecord:
    """Normalise *variant* against the FASTA at *fasta_path* (kept for API
    compatibility; batch callers should open the FASTA once and use
    ``normalize_record``)."""
    if variant.is_snv:
        return variant
    fasta = open_fasta(fasta_path)
    try:
        return normalize_record(variant, fasta)
    finally:
        if fasta is not None:
            fasta.close()
