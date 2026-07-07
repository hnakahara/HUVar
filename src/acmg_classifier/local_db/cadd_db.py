"""tabix query for CADD precomputed scores (auxiliary PP3/BP4)."""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional

import structlog

from acmg_classifier.models.annotation import CADDData
from acmg_classifier.utils.tabix import TabixReader, open_tabix, fetch_region

log = structlog.get_logger()


def _match_cadd(lines: Iterable[str], pos: int, ref: str, alt: str) -> Optional[CADDData]:
    """Pick the CADD PHRED from tabix region lines matching (pos, ref, alt)."""
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        f_pos, f_ref, f_alt = fields[1], fields[2], fields[3]
        if f_ref != ref or f_alt != alt:
            continue
        try:
            if int(f_pos) != pos:
                continue
            phred = float(fields[4])
        except ValueError:
            # A single malformed row must not break annotation.
            continue
        return CADDData(phred=phred)
    return None


def query_cadd(
    tsv_gz_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> Optional[CADDData]:
    """Query the CADD TSV (tabix-indexed) for one variant.

    scripts/setup_data.py normalises CADD into a 5-column, tab-separated,
    bgzipped file indexed on the assembly's position column:

        chrom  pos  ref  alt  CADD_PHRED

    A row is matched on (pos, ref, alt); the 5th column is the PHRED-scaled
    score. Returns None when the file is absent, the variant is uncovered, or the
    lookup raises — the criteria then simply skip the CADD branch rather than
    aborting annotation.
    """
    if not tsv_gz_path.exists():
        log.warning("cadd_missing", path=str(tsv_gz_path))
        return None

    try:
        with open_tabix(tsv_gz_path) as tf:
            return _match_cadd(fetch_region(tf, chrom, pos, pos), pos, ref, alt)
    except Exception as exc:
        log.error("cadd_error", error=str(exc))
    return None


def query_cadd_reader(
    reader: TabixReader,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> Optional[CADDData]:
    """Batch variant of query_cadd using a persistent thread-local handle."""
    if not reader.exists():
        log.warning("cadd_missing", path=str(reader.path))
        return None
    try:
        return _match_cadd(reader.fetch(chrom, pos, pos), pos, ref, alt)
    except Exception as exc:
        log.error("cadd_error", error=str(exc))
    return None
