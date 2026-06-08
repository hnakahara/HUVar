"""phyloP conservation lookup for the BP7 "not highly conserved" gate.

Reads per-base phyloP100way scores from a UCSC bigWig (commercial-use OK). The
bigWig is large (~9 GB) and OPTIONAL: when absent the reader reports
"unavailable" and BP7 falls back to its splice-only logic. The reader library
(pyBigWig) is a default dependency, so the gate activates automatically once the
bigWig file is present.

The track is downloaded by default by setup_data.py (pass --skip-phylop to opt
out). Source: hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/
"""
from __future__ import annotations

import threading
from pathlib import Path
from typing import Optional

import structlog

log = structlog.get_logger()


class PhyloPReader:
    """Random-access phyloP100way reader over a UCSC bigWig.

    Construction never raises: a missing file or a missing pyBigWig dependency
    leaves the reader unavailable so callers degrade gracefully.

    The bigWig (~9 GB) is opened LAZILY on the first lookup, not at construction.
    BP7 only consults conservation for synonymous / deep-intronic variants whose
    splice score is benign, so most runs (e.g. an `explain` on a missense or LoF
    variant) never need it — deferring the open avoids the multi-second,
    multi-GB read entirely for those runs."""

    def __init__(self, bigwig_path: Optional[Path]) -> None:
        self._path = bigwig_path
        self._bw = None
        self._opened = False          # whether an open was attempted
        self._lock = threading.Lock()  # _ensure_open runs under the annotation pool
        # "Loadable" = we have a path and pyBigWig imports. The actual (slow)
        # file open is deferred to the first lookup.
        self._loadable = bigwig_path is not None
        if bigwig_path is None:
            return
        try:
            import pyBigWig  # type: ignore # noqa: F401
        except ImportError:
            # pyBigWig is a default dependency, so a missing import means a
            # broken/partial install (or a build failure of its libcurl/zlib
            # extension) rather than an opt-in the user skipped.
            log.warning(
                "phylop_unavailable",
                reason="pyBigWig not importable",
                hint='reinstall to restore the BP7 conservation gate: pip install pyBigWig',
            )
            self._loadable = False

    def _ensure_open(self) -> None:
        """Open the bigWig on first use (thread-safe, runs once)."""
        if self._opened:
            return
        with self._lock:
            if self._opened:
                return
            try:
                if self._loadable:
                    import pyBigWig  # type: ignore
                    from acmg_classifier.utils import progress
                    log.info("phylop_loading", path=str(self._path))
                    with progress.status("Loading phyloP conservation track…"):
                        self._bw = pyBigWig.open(str(self._path))
            except Exception as exc:  # noqa: BLE001 — any open failure → unavailable
                log.warning("phylop_open_failed", path=str(self._path), error=str(exc))
                self._bw = None
                self._loadable = False
            finally:
                self._opened = True

    def is_available(self) -> bool:
        """True if conservation scores can be served — either already open, or
        openable on demand (path present + pyBigWig importable + no prior open
        failure). Does NOT trigger the open, so callers can cheaply gate on it."""
        return self._bw is not None or (self._loadable and not self._opened)

    def value(self, chrom: str, pos: int) -> Optional[float]:
        """phyloP score at a 1-based position, or None if unavailable/missing.

        UCSC bigWig intervals are 0-based half-open, so a 1-based ``pos`` maps to
        the half-open interval [pos-1, pos)."""
        self._ensure_open()
        if self._bw is None:
            return None
        c = chrom if chrom.startswith("chr") else f"chr{chrom}"
        try:
            vals = self._bw.values(c, pos - 1, pos)
        except (RuntimeError, OverflowError):
            return None
        if not vals:
            return None
        v = vals[0]
        # pyBigWig returns NaN for positions with no score.
        return None if v != v else float(v)
