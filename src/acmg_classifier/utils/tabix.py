from __future__ import annotations
import threading
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


class TabixReader:
    """Thread-local, persistent ``pysam.TabixFile`` for one bgzipped file.

    Opening a TabixFile reads and parses the ``.tbi`` index; doing that per
    variant is the dominant cost when annotating hundreds of thousands of
    variants. This holds one open handle **per worker thread** (a pysam handle
    is not safe for concurrent fetch, and the annotation thread pool calls
    fetch() from several threads, so we never share one across threads) and
    reuses it for every lookup. ``fetch`` reproduces ``fetch_region`` exactly,
    including the chr-prefix fallback, so reader-based queries return the same
    lines as the connection-per-call path."""

    def __init__(self, path: Path) -> None:
        self._path = path
        self._local = threading.local()
        self._exists: bool | None = None

    @property
    def path(self) -> Path:
        return self._path

    def exists(self) -> bool:
        """Whether the backing file exists, checked once and memoised.

        query_*_reader calls this on every variant; a reference DB never
        appears or disappears mid-run, so a single ``stat`` is enough. Caching
        it avoids hundreds of thousands of redundant filesystem stats (a real
        cost on network mounts, where it dominated after the open-per-variant
        cost was removed)."""
        if self._exists is None:
            self._exists = self._path.exists()
        return self._exists

    def _handle(self) -> pysam.TabixFile:
        tf = getattr(self._local, "tf", None)
        if tf is None:
            tf = pysam.TabixFile(str(self._path))
            self._local.tf = tf
        return tf

    def fetch(self, chrom: str, start: int, end: int) -> Iterator[str]:
        """Yield lines overlapping [start, end] (1-based inclusive) using the
        thread-local handle, with the identical chr-prefix fallback as
        ``fetch_region``."""
        tf = self._handle()
        try:
            yield from tf.fetch(chrom, start - 1, end)
        except ValueError:
            alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            try:
                yield from tf.fetch(alt_chrom, start - 1, end)
            except ValueError:
                return
