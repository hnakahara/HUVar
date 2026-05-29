"""tabix query for AlphaMissense precomputed scores (PP3/BP4)."""
from __future__ import annotations
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.models.annotation import AlphaMissenseData
from acmg_classifier.utils.tabix import open_tabix, fetch_region

log = structlog.get_logger()


def query_alphamissense(
    tsv_gz_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> Optional[AlphaMissenseData]:
    """
    Query AlphaMissense TSV (tabix-indexed).

    Column order (AlphaMissense hg38/hg19 TSV):
    #CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant
    am_pathogenicity  am_class
    """
    if not tsv_gz_path.exists():
        log.warning("alphamissense_missing", path=str(tsv_gz_path))
        return None

    try:
        with open_tabix(tsv_gz_path) as tf:
            # AlphaMissense files contain one row per (REF, ALT) combination
            # at each position, so we must filter on all three to pick the
            # correct allele. The header row begins with '#' and is skipped.
            for line in fetch_region(tf, chrom, pos, pos):
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 9:
                    continue
                _, f_pos, f_ref, f_alt = fields[0], fields[1], fields[2], fields[3]
                if int(f_pos) != pos or f_ref != ref or f_alt != alt:
                    continue
                try:
                    score = float(fields[8])
                except (ValueError, IndexError):
                    # Treat unparseable rows as "no score" rather than aborting
                    # — a corrupt single row should not break annotation.
                    continue
                am_class = fields[9].strip() if len(fields) > 9 else None
                return AlphaMissenseData(score=score, classification=am_class)
    except Exception as exc:
        log.error("alphamissense_error", error=str(exc))
    return None
