"""phyloP conservation lookup for the BP7 "not highly conserved" gate.

Reads per-base phyloP100way scores from a UCSC bigWig (commercial-use OK). The
bigWig is large (~9 GB) and OPTIONAL: when absent the reader reports
"unavailable" and BP7 falls back to its splice-only logic. The reader library
(pyBigWig) is a default dependency, so the gate activates automatically once the
bigWig file is present.

Download the track: hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/
(or run: setup_data.py --with-phylop)
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import structlog

log = structlog.get_logger()


class PhyloPReader:
    """Random-access phyloP100way reader over a UCSC bigWig.

    Construction never raises: a missing file or a missing pyBigWig dependency
    leaves the reader unavailable so callers degrade gracefully.
    """

    def __init__(self, bigwig_path: Optional[Path]) -> None:
        self._bw = None
        self._chrom_prefix = ""  # UCSC bigWigs use "chr1"; tracked for lookups
        if bigwig_path is None:
            return
        try:
            import pyBigWig  # type: ignore
        except ImportError:
            # pyBigWig is a default dependency, so a missing import means a
            # broken/partial install (or a build failure of its libcurl/zlib
            # extension) rather than an opt-in the user skipped.
            log.warning(
                "phylop_unavailable",
                reason="pyBigWig not importable",
                hint='reinstall to restore the BP7 conservation gate: pip install pyBigWig',
            )
            return
        try:
            self._bw = pyBigWig.open(str(bigwig_path))
        except Exception as exc:  # noqa: BLE001 — any open failure → unavailable
            log.warning("phylop_open_failed", path=str(bigwig_path), error=str(exc))
            self._bw = None

    def is_available(self) -> bool:
        return self._bw is not None

    def value(self, chrom: str, pos: int) -> Optional[float]:
        """phyloP score at a 1-based position, or None if unavailable/missing.

        UCSC bigWig intervals are 0-based half-open, so a 1-based ``pos`` maps to
        the half-open interval [pos-1, pos)."""
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
