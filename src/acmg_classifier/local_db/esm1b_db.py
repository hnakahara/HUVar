"""SQLite query for ESM1b precomputed LLR scores (Brandes et al. 2023).

Schema produced by `scripts/setup_data.py` from
`ALL_hum_isoforms_ESM1b_LLR.zip`:

    CREATE TABLE scores (
        transcript_id TEXT NOT NULL,
        aa_pos        INTEGER NOT NULL,
        alt_aa        TEXT NOT NULL,
        llr           REAL NOT NULL,
        PRIMARY KEY (transcript_id, aa_pos, alt_aa)
    );
    CREATE INDEX idx_scores_tx_pos ON scores(transcript_id, aa_pos);

`transcript_id` is the Ensembl transcript ID (ENST...) with the version
suffix stripped (e.g. ENST00000357654). The builder maps Brandes's
UniProt-keyed scores to ENST via the bundled `isoform_list.csv`.
"""
from __future__ import annotations
import sqlite3
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.models.annotation import ESM1bData

log = structlog.get_logger()


def _strip_version(transcript_id: str) -> str:
    """ENST00000357654.4 → ENST00000357654."""
    return transcript_id.split(".", 1)[0]


def query_esm1b(
    sqlite_path: Path,
    transcript_id: str,
    aa_pos: int,
    alt_aa: str,
) -> Optional[ESM1bData]:
    """Look up ESM1b LLR for a missense variant.

    Returns None if the database is missing, the transcript/position is not
    covered, or the lookup raises. A None result causes PP3/BP4 to skip the
    ESM1b branch — never raise upstream.
    """
    if not sqlite_path.exists():
        log.warning("esm1b_missing", path=str(sqlite_path))
        return None
    if not transcript_id or aa_pos is None or not alt_aa:
        return None

    tx = _strip_version(transcript_id)
    try:
        conn = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
        try:
            cur = conn.execute(
                "SELECT llr FROM scores "
                "WHERE transcript_id = ? AND aa_pos = ? AND alt_aa = ? LIMIT 1",
                (tx, int(aa_pos), alt_aa),
            )
            row = cur.fetchone()
        finally:
            conn.close()
    except sqlite3.Error as exc:
        log.error("esm1b_sqlite_error", error=str(exc), tx=tx, aa_pos=aa_pos, alt_aa=alt_aa)
        return None

    if row is None:
        return None
    return ESM1bData(llr=float(row[0]))
