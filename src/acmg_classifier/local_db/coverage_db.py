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

from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.utils.chrom import chrom_candidates

log = structlog.get_logger()


class CoverageDB:
    def __init__(self, db_path: Path) -> None:
        self._db_path = db_path

    @property
    def available(self) -> bool:
        return self._db_path.exists()

    def mean_depth(self, chrom: str, start: int, end: Optional[int] = None) -> Optional[float]:
        """Average gnomAD mean read depth over ``[start, end]`` (end defaults to
        start) on *chrom*, or ``None`` when the DB is absent or no covered locus
        falls in the span (depth unknown → caller skips the gate)."""
        if not self._db_path.exists():
            return None
        c1, c2 = chrom_candidates(chrom)
        hi = end if end is not None else start
        try:
            import duckdb
            con = duckdb.connect(str(self._db_path), read_only=True)
            row = con.execute(
                """
                SELECT avg(mean_dp) FROM coverage
                WHERE chrom IN (?, ?) AND pos BETWEEN ? AND ?
                """,
                [c1, c2, start, hi],
            ).fetchone()
            con.close()
        except Exception as exc:
            log.error("coverage_query_error", error=str(exc))
            return None
        return float(row[0]) if row and row[0] is not None else None
