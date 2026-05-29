"""SQUIRLS splice predictor (Apache 2.0, default). Queries precomputed SQLite DB."""
from __future__ import annotations
import sqlite3
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.local_db.splice.base import SplicePredictor
from acmg_classifier.models.annotation import SpliceScore
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.utils.chrom import chrom_candidates

log = structlog.get_logger()


class SquirlsPredictor(SplicePredictor):
    """
    Queries the SQUIRLS precomputed SQLite database.

    SQUIRLS DB schema (squirls-2309-hg38/hg19):
      Table: squirls_scores
      Columns: chrom, pos, ref, alt, squirls_score, ...

    Thresholds are APPROXIMATE (not Walker 2023 calibrated).
    Output reports note this limitation.
    """

    def __init__(self, db_dir: Path) -> None:
        self._db_dir = db_dir
        self._db_path: Optional[Path] = self._find_db()

    def _find_db(self) -> Optional[Path]:
        """Locate the SQUIRLS SQLite file inside the configured directory.

        SQUIRLS releases ship under different filenames per build year
        (e.g. squirls-2309-hg38.db). Globbing for *.db then *.sqlite lets
        the same code work across releases without per-version config."""
        if not self._db_dir.exists():
            return None
        for candidate in self._db_dir.glob("*.db"):
            return candidate
        for candidate in self._db_dir.glob("*.sqlite"):
            return candidate
        return None

    def is_available(self) -> bool:
        return self._db_path is not None and self._db_path.exists()

    def predict(self, variant: VariantRecord) -> SpliceScore:
        if not self.is_available():
            return SpliceScore(tool="squirls", is_available=False)

        c1, c2 = chrom_candidates(variant.chrom)
        try:
            con = sqlite3.connect(str(self._db_path))
            row = con.execute(
                "SELECT squirls_score FROM squirls_scores "
                "WHERE chrom IN (?, ?) AND pos = ? AND ref = ? AND alt = ? LIMIT 1",
                (c1, c2, variant.pos, variant.ref, variant.alt),
            ).fetchone()
            con.close()
            if row is None:
                return SpliceScore(tool="squirls", is_available=True, raw_score=None)
            return SpliceScore(tool="squirls", is_available=True, raw_score=float(row[0]))
        except Exception as exc:
            log.error("squirls_error", error=str(exc))
            return SpliceScore(tool="squirls", is_available=False)
