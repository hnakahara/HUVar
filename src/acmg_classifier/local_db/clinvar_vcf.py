"""tabix-indexed ClinVar VCF query.

Originally created to feed PP5 (which is now retired — see pp5.py). The
loader is still used by the orchestrator to populate
AnnotationData.clinvar_vcf so the per-variant TSV/JSON output can include
the public ClinVar classification for human review, even though PP5 itself
is never automatically applied.
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterable

import structlog

from acmg_classifier.models.annotation import ClinVarRecord
from acmg_classifier.utils.tabix import TabixReader, open_tabix, fetch_region

log = structlog.get_logger()

# ClinVar review-status to star-rating mapping. Star ratings are the
# numeric proxy for assertion confidence used throughout the codebase
# (ACMG/ClinGen require ≥1 star for evidence-based criteria like PS1/PM5).
# Keys must match ClinVar CLNREVSTAT *lowercase* substrings so that
# punctuation/casing changes in upstream files do not break the lookup.
_STAR_MAP = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
}


def _parse_info(info_str: str) -> dict[str, str]:
    """Parse a VCF INFO column into a plain dict.

    Flag-style fields (no `=`) become `"1"` so callers can use a uniform
    dict[str, str] interface instead of dict[str, str | None]."""
    d: dict[str, str] = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, _, v = token.partition("=")
            d[k] = v
        else:
            d[token] = "1"
    return d


def _star_rating(review_status: str) -> int:
    """Map ClinVar review-status text to a 0-4 star rating.

    Uses substring match (not exact equality) because ClinVar occasionally
    appends conflict descriptors to the canonical strings. Returns 0 when
    no pattern matches so unknown statuses cannot inadvertently trigger
    star-gated criteria."""
    for pattern, stars in _STAR_MAP.items():
        if pattern in review_status.lower():
            return stars
    return 0


def _match_clinvar(lines: Iterable[str], pos: int, ref: str, alt: str) -> list[ClinVarRecord]:
    """Collect all ClinVar records from tabix region lines exactly matching
    (pos, ref, alt). Shared by query_clinvar_vcf and its reader variant.

    fetch returns ALL variants overlapping [pos, pos]; ALT alleles must be
    filtered manually because tabix has no concept of allele identity."""
    records: list[ClinVarRecord] = []
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        vcf_chrom, vcf_pos_str, vcf_id, vcf_ref, vcf_alt = fields[:5]
        if int(vcf_pos_str) != pos or vcf_ref != ref or vcf_alt != alt:
            continue
        info = _parse_info(fields[7])
        # ClinVar VCF uses underscores in INFO values (e.g. "Likely_pathogenic")
        # — restore spaces to match the canonical sig strings used elsewhere.
        clnsig = info.get("CLNSIG", "").replace("_", " ")
        clnrevstat = info.get("CLNREVSTAT", "").replace("_", " ")
        stars = _star_rating(clnrevstat)
        records.append(ClinVarRecord(
            variation_id=vcf_id if vcf_id != "." else None,
            clinical_significance=clnsig,
            review_status=clnrevstat,
            star_rating=stars,
        ))
    return records


def query_clinvar_vcf(
    vcf_gz_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> list[ClinVarRecord]:
    """Fetch all ClinVar VCF records exactly matching the (chrom, pos, ref, alt)
    quadruple via tabix region lookup.

    A single position can carry multiple ClinVar records (different
    submitters, different alleles), so we filter the tabix region by exact
    REF/ALT match. Returns an empty list on any error — annotation should
    not abort the pipeline because ClinVar is unavailable."""
    if not vcf_gz_path.exists():
        log.warning("clinvar_vcf_missing", path=str(vcf_gz_path))
        return []

    try:
        with open_tabix(vcf_gz_path) as tf:
            return _match_clinvar(fetch_region(tf, chrom, pos, pos), pos, ref, alt)
    except Exception as exc:
        log.error("clinvar_vcf_error", error=str(exc))
    return []


def query_clinvar_vcf_reader(
    reader: TabixReader,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> list[ClinVarRecord]:
    """Batch variant of query_clinvar_vcf using a persistent thread-local handle."""
    if not reader.exists():
        log.warning("clinvar_vcf_missing", path=str(reader.path))
        return []
    try:
        return _match_clinvar(reader.fetch(chrom, pos, pos), pos, ref, alt)
    except Exception as exc:
        log.error("clinvar_vcf_error", error=str(exc))
    return []
