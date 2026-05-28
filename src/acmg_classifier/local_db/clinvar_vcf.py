"""tabix-indexed ClinVar VCF query for PP5 (P/LP classification lookup)."""
from __future__ import annotations
from pathlib import Path

import structlog

from acmg_classifier.models.annotation import ClinVarRecord
from acmg_classifier.utils.tabix import open_tabix, fetch_region

log = structlog.get_logger()

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
    d: dict[str, str] = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, _, v = token.partition("=")
            d[k] = v
        else:
            d[token] = "1"
    return d


def _star_rating(review_status: str) -> int:
    for pattern, stars in _STAR_MAP.items():
        if pattern in review_status.lower():
            return stars
    return 0


def query_clinvar_vcf(
    vcf_gz_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> list[ClinVarRecord]:
    if not vcf_gz_path.exists():
        log.warning("clinvar_vcf_missing", path=str(vcf_gz_path))
        return []

    records: list[ClinVarRecord] = []
    try:
        with open_tabix(vcf_gz_path) as tf:
            for line in fetch_region(tf, chrom, pos, pos):
                fields = line.split("\t")
                if len(fields) < 8:
                    continue
                vcf_chrom, vcf_pos_str, vcf_id, vcf_ref, vcf_alt = fields[:5]
                if int(vcf_pos_str) != pos or vcf_ref != ref or vcf_alt != alt:
                    continue
                info = _parse_info(fields[7])
                clnsig = info.get("CLNSIG", "").replace("_", " ")
                clnrevstat = info.get("CLNREVSTAT", "").replace("_", " ")
                stars = _star_rating(clnrevstat)
                records.append(ClinVarRecord(
                    variation_id=vcf_id if vcf_id != "." else None,
                    clinical_significance=clnsig,
                    review_status=clnrevstat,
                    star_rating=stars,
                ))
    except Exception as exc:
        log.error("clinvar_vcf_error", error=str(exc))
    return records
