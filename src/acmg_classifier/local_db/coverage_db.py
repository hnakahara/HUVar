"""gnomAD per-locus coverage (mean read depth) lookup for the PM2 depth gate.

The ENIGMA BRCA1/2 VCEP requires "the region around the variant must have an
average read depth >= 25" before PM2 (absence) can be applied — so a variant
that is "absent" merely because the region is poorly covered does not earn PM2.

Reads a DuckDB built by ``scripts/build_gnomad_coverage.py`` from the gnomAD
exomes coverage summary (table ``coverage(chrom, pos, mean_dp)``). A missing DB
degrades to "depth unknown" (``None``), so the evaluator skips the gate rather
than blocking — identical graceful-degradation to the phyloP / af_xy paths.
"""
from __future__ import annotations

import os
import threading
from functools import lru_cache
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.utils.chrom import chrom_candidates

log = structlog.get_logger()

# Every DuckDB connection opened by CoverageDB, so they can all be closed right
# before a fork(): a DuckDB database instance owns worker threads, and a child
# that inherited the instance (DuckDB re-uses an open instance for the same
# file) would wait on threads that do not exist in it. Closing them first makes
# the parent re-open lazily and each child open its own instance.
_OPEN: list = []
_OPEN_LOCK = threading.Lock()
_GENERATION = 0


def close_all() -> None:
    """Close every open coverage connection (they re-open on next use)."""
    global _GENERATION
    with _OPEN_LOCK:
        for con in _OPEN:
            try:
                con.close()
            except Exception:  # noqa: BLE001 — closing is best effort
                pass
        _OPEN.clear()
        _GENERATION += 1


def _forget_after_fork() -> None:
    global _OPEN_LOCK, _GENERATION
    _OPEN_LOCK = threading.Lock()
    _OPEN.clear()
    _GENERATION += 1


if hasattr(os, "register_at_fork"):
    os.register_at_fork(before=close_all, after_in_child=_forget_after_fork)


class CoverageDB:
    """Read-only gnomAD coverage lookups over one persistent connection.

    ``mean_depth`` used to open (and close) a DuckDB connection per call; for a
    BRCA1/2 batch every PM2 candidate pays that open, which dominated the
    classification time of large runs. One connection is now kept per thread
    and per process (re-opened after ``fork`` — a DuckDB handle must never be
    shared with a child), and results are memoised per span."""

    def __init__(self, db_path: Path) -> None:
        self._db_path = db_path
        self._exists: bool | None = None
        self._local = threading.local()
        # Per-instance memo (bound method cache would pin ``self``).
        self._memo = lru_cache(maxsize=65536)(self._query)

    @property
    def available(self) -> bool:
        # A reference DB does not appear or vanish mid-run: stat it once.
        if self._exists is None:
            self._exists = self._db_path.exists()
        return self._exists

    def _connection(self):
        con = getattr(self._local, "con", None)
        if con is None or getattr(self._local, "gen", None) != _GENERATION:
            import duckdb
            con = duckdb.connect(str(self._db_path), read_only=True)
            with _OPEN_LOCK:
                _OPEN.append(con)
            self._local.con, self._local.gen = con, _GENERATION
        return con

    def _query(self, c1: str, c2: str, start: int, hi: int) -> Optional[float]:
        try:
            row = self._connection().execute(
                """
                SELECT avg(mean_dp) FROM coverage
                WHERE chrom IN (?, ?) AND pos BETWEEN ? AND ?
                """,
                [c1, c2, start, hi],
            ).fetchone()
        except Exception as exc:
            log.error("coverage_query_error", error=str(exc))
            return None
        return float(row[0]) if row and row[0] is not None else None

    def mean_depth(self, chrom: str, start: int, end: Optional[int] = None) -> Optional[float]:
        """Average gnomAD mean read depth over ``[start, end]`` (end defaults to
        start) on *chrom*, or ``None`` when the DB is absent or no covered locus
        falls in the span (depth unknown → caller skips the gate)."""
        if not self.available:
            return None
        c1, c2 = chrom_candidates(chrom)
        hi = end if end is not None else start
        return self._memo(c1, c2, start, hi)
