"""tabix query for Dfam-based RepeatMasker BED (BP3 repeat region check)."""
from __future__ import annotations
from pathlib import Path
from typing import Iterable

import structlog

from acmg_classifier.models.annotation import RepeatMaskerRegion
from acmg_classifier.utils.tabix import TabixReader, open_tabix, fetch_region

log = structlog.get_logger()


def _match_repeat(lines: Iterable[str]) -> RepeatMaskerRegion:
    """First overlapping repeat element from tabix region lines, or "not in
    repeat". Shared by query_repeat and its reader variant.

    First overlapping record wins — repeat elements rarely overlap each other,
    and even when they do the answer to "is this position in a repeat?" doesn't
    change after the first hit."""
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        rep_name = fields[3] if len(fields) > 3 else None
        rep_class = fields[6] if len(fields) > 6 else None
        return RepeatMaskerRegion(
            in_repeat=True,
            repeat_class=rep_class,
            repeat_name=rep_name,
        )
    return RepeatMaskerRegion(in_repeat=False)


def query_repeat(
    bed_gz_path: Path,
    chrom: str,
    pos: int,
) -> RepeatMaskerRegion:
    """
    Check if position falls within a Dfam repeat element.

    BED format: chrom  start  end  name  score  strand  repClass  repFamily
    (UCSC RepeatMasker BED-like, but with Dfam annotations).
    """
    # Missing repeat BED is non-fatal — we return "not in repeat" so PM4
    # cannot be silently down-graded to BP3 by an absent data file. The
    # warning is emitted once via structlog rather than raising.
    if not bed_gz_path.exists():
        log.warning("repeatmasker_missing", path=str(bed_gz_path))
        return RepeatMaskerRegion(in_repeat=False)

    try:
        with open_tabix(bed_gz_path) as tf:
            return _match_repeat(fetch_region(tf, chrom, pos, pos))
    except Exception as exc:
        log.error("repeatmasker_error", error=str(exc))
    return RepeatMaskerRegion(in_repeat=False)


def query_repeat_reader(
    reader: TabixReader,
    chrom: str,
    pos: int,
) -> RepeatMaskerRegion:
    """Batch variant of query_repeat using a persistent thread-local handle."""
    if not reader.exists():
        log.warning("repeatmasker_missing", path=str(reader.path))
        return RepeatMaskerRegion(in_repeat=False)
    try:
        return _match_repeat(reader.fetch(chrom, pos, pos))
    except Exception as exc:
        log.error("repeatmasker_error", error=str(exc))
    return RepeatMaskerRegion(in_repeat=False)
