"""tabix query for REVEL precomputed scores (PP3/BP4)."""
from __future__ import annotations
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.models.annotation import RevelData
from acmg_classifier.utils.tabix import open_tabix, fetch_region

log = structlog.get_logger()


def query_revel(
    tsv_gz_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> Optional[RevelData]:
    """Query the REVEL TSV (tabix-indexed) for one SNV.

    scripts/setup_data.py normalises REVEL into a 5-column, tab-separated,
    bgzipped file indexed on the assembly's position column:

        chrom  pos  ref  alt  REVEL

    REVEL covers only single-nucleotide missense substitutions, so a row is
    matched on (pos, ref, alt). Returns None when the file is absent, the
    variant is uncovered, or the lookup raises — PP3/BP4 then simply skip the
    REVEL branch rather than aborting annotation.
    """
    if not tsv_gz_path.exists():
        log.warning("revel_missing", path=str(tsv_gz_path))
        return None

    try:
        with open_tabix(tsv_gz_path) as tf:
            for line in fetch_region(tf, chrom, pos, pos):
                fields = line.split("\t")
                if len(fields) < 5:
                    continue
                f_pos, f_ref, f_alt = fields[1], fields[2], fields[3]
                if f_ref != ref or f_alt != alt:
                    continue
                try:
                    if int(f_pos) != pos:
                        continue
                    score = float(fields[4])
                except ValueError:
                    # A single malformed row must not break annotation.
                    continue
                return RevelData(score=score)
    except Exception as exc:
        log.error("revel_error", error=str(exc))
    return None
