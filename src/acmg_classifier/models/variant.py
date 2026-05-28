from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, field_validator

from acmg_classifier.models.enums import Assembly, VariantType


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
        return v if v.startswith("chr") else f"chr{v}"

    @field_validator("ref", "alt")
    @classmethod
    def upper_bases(cls, v: str) -> str:
        return v.upper()

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def variant_type(self) -> VariantType:
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
        if self.variant_id is None:
            self.variant_id = self.key
