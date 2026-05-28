from __future__ import annotations
from contextlib import contextmanager
from pathlib import Path
from typing import Generator, Iterator
import pysam  # type: ignore


@contextmanager
def open_tabix(path: Path) -> Generator[pysam.TabixFile, None, None]:
    """Context manager for a pysam tabix file."""
    tf = pysam.TabixFile(str(path))
    try:
        yield tf
    finally:
        tf.close()


def fetch_region(
    tf: pysam.TabixFile,
    chrom: str,
    start: int,
    end: int,
) -> Iterator[str]:
    """
    Yield raw tab-separated lines overlapping [start, end] (1-based, inclusive).
    Silently skips if contig not in index.
    """
    try:
        yield from tf.fetch(chrom, start - 1, end)
    except ValueError:
        # Contig not found in tabix index — try without 'chr' prefix or vice versa
        alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
        try:
            yield from tf.fetch(alt_chrom, start - 1, end)
        except ValueError:
            return
