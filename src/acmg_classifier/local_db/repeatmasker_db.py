"""tabix query for Dfam-based RepeatMasker BED (BP3 repeat region check)."""
from __future__ import annotations
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.models.annotation import RepeatMaskerRegion
from acmg_classifier.utils.tabix import open_tabix, fetch_region

log = structlog.get_logger()


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
    if not bed_gz_path.exists():
        log.warning("repeatmasker_missing", path=str(bed_gz_path))
        return RepeatMaskerRegion(in_repeat=False)

    try:
        with open_tabix(bed_gz_path) as tf:
            for line in fetch_region(tf, chrom, pos, pos):
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
    except Exception as exc:
        log.error("repeatmasker_error", error=str(exc))
    return RepeatMaskerRegion(in_repeat=False)
