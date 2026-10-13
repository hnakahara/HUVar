from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, field_validator

from acmg_classifier.models.enums import Assembly, VariantType

import re as _re

# Indels longer than this (bp, |len(REF) - len(ALT)| or the longer allele) are at
# the size of structural variants: they are still classified, but flagged, because
# gnomAD/ClinVar allele matching, in-silico predictors and the PVS1/PM4 rules are
# designed for small variants. (Reviewer 1, Genome Medicine revision.)
LARGE_INDEL_BP = 50
_BASES = _re.compile(r"^[ACGTN]+$")


def unsupported_allele_reason(ref: str, alt: str) -> Optional[str]:
    """Reason an allele cannot be classified, or None if it is a supported
    sequence variant (SNV, MNV or indel written with explicit bases)."""
    a = (alt or "").upper()
    if not a or a == ".":
        return "no ALT allele (ALT='.'): not a variant, excluded from annotation/classification"
    if a.startswith("<") and a.endswith(">"):
        return (f"symbolic ALT allele {alt} (CNV/structural variant): not supported, "
                "excluded from classification; use a dedicated CNV/SV classifier")
    if "[" in a or "]" in a:
        return "breakend ALT allele (structural variant): not supported, excluded from classification"
    if a == "*":
        return "spanning-deletion ALT allele '*': classified with the overlapping deletion record, excluded here"
    if not _BASES.match(a) or not _BASES.match((ref or "").upper()):
        return f"non-nucleotide allele (REF={ref}, ALT={alt}): not supported, excluded from classification"
    return None


def large_indel_warning(ref: str, alt: str) -> Optional[str]:
    size = max(len(ref or ""), len(alt or ""))
    if len(ref or "") != len(alt or "") and size > LARGE_INDEL_BP:
        return (f"LARGE_INDEL: {size} bp allele (> {LARGE_INDEL_BP} bp, structural-variant size); "
                "population frequency, in-silico and splice evidence may be unavailable and PVS1/PM4 "
                "assume small variants — review with a CNV/SV-aware tool")
    return None


class VariantRecord(BaseModel):
    """Normalised representation of a single variant from a VCF record."""

    chrom: str
    pos: int          # 1-based, left-normalised
    ref: str
    alt: str
    assembly: Assembly

    # Optional VCF fields
    variant_id: Optional[str] = None   # CHROM:POS:REF:ALT key (set after normalisation)
    vcf_id: Optional[str] = None       # ID column from VCF
    qual: Optional[float] = None
    filter: Optional[str] = None
    sample_id: Optional[str] = None

    @field_validator("chrom")
    @classmethod
    def normalise_chrom(cls, v: str) -> str:
        # GRCh38 source data is inconsistent: gnomAD uses "chr1", legacy VCFs
        # use "1". Normalise to UCSC-style "chr" prefix so every downstream
        # lookup (tabix, SQLite key, RepeatMasker BED) sees the same form.
        return v if v.startswith("chr") else f"chr{v}"

    @field_validator("ref", "alt")
    @classmethod
    def upper_bases(cls, v: str) -> str:
        # VCF allows lowercase bases; uppercase is required for byte-equal
        # comparison against reference FASTA and ClinVar/gnomAD records.
        return v.upper()

    @property
    def key(self) -> str:
        # Canonical CHROM:POS:REF:ALT identifier used as the join key across
        # every per-variant annotation source. Keeping this as a single string
        # lets us index dicts and SQLite without composite keys.
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def variant_type(self) -> VariantType:
        # PVS1/PM4/BP3 dispatch on variant type, so we classify by length:
        # 1bp↔1bp is SNV, unequal-length is INDEL, equal-length >1bp is MNV.
        if len(self.ref) == 1 and len(self.alt) == 1:
            return VariantType.SNV
        if len(self.ref) != len(self.alt):
            return VariantType.INDEL
        return VariantType.MNV

    @property
    def is_snv(self) -> bool:
        return self.variant_type == VariantType.SNV

    @property
    def is_indel(self) -> bool:
        return self.variant_type == VariantType.INDEL

    def model_post_init(self, __context: object) -> None:
        # Default variant_id to the canonical key so every variant has a stable
        # identifier even when the VCF ID column is "." (the common case for
        # research/clinical VCFs without dbSNP rsIDs).
        if self.variant_id is None:
            self.variant_id = self.key
