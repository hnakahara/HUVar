"""SQLite query for ESM1b precomputed LLR scores (Brandes et al. 2023).

Schema produced by `scripts/setup_data.py` from
`ALL_hum_isoforms_ESM1b_LLR.zip`:

    CREATE TABLE scores (
        uniprot_id TEXT NOT NULL,
        aa_pos     INTEGER NOT NULL,
        alt_aa     TEXT NOT NULL,
        llr        REAL NOT NULL,
        PRIMARY KEY (uniprot_id, aa_pos, alt_aa)
    );
    CREATE INDEX idx_scores_uni_pos ON scores(uniprot_id, aa_pos);
    CREATE TABLE aliases (entry_name TEXT PRIMARY KEY, uniprot_id TEXT NOT NULL);

`uniprot_id` is the SwissProt/TrEMBL accession that the Brandes archive
uses as its file-name key (e.g. `P38398`, or `P38398-2` for alternative
isoforms). VEP --uniprot usually emits the same accession (the version
suffix "P38398.4" is stripped upstream in vep_runner._parse_transcript),
but some caches (observed: GRCh37 release 111) emit the UniProt entry name
/ mnemonic instead (e.g. `PK3CD_HUMAN`). `_resolve_accession` normalises
those back to the accession via the `aliases` table before querying.
"""
from __future__ import annotations
import sqlite3
import threading
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.models.annotation import ESM1bData

log = structlog.get_logger()


def _select_llr(
    conn: sqlite3.Connection, uniprot_id: str, aa_pos: int, alt_aa: str
) -> Optional[ESM1bData]:
    """Resolve the accession and point-select the LLR on an open connection.

    Shared by the connection-per-call query_esm1b() and the persistent-
    connection Esm1bDB.lookup() so both do the identical accession resolution
    and lookup — only the connection lifecycle differs."""
    accession = _resolve_accession(conn, uniprot_id)
    row = conn.execute(
        "SELECT llr FROM scores "
        "WHERE uniprot_id = ? AND aa_pos = ? AND alt_aa = ? LIMIT 1",
        (accession, int(aa_pos), alt_aa),
    ).fetchone()
    if row is None:
        return None
    return ESM1bData(llr=float(row[0]))


def _resolve_accession(conn: sqlite3.Connection, uniprot_id: str) -> str:
    """Normalise a VEP UniProt token to the accession the scores table uses.

    Most caches emit the accession directly (e.g. ``O00329``) and are returned
    unchanged. Some caches (observed: GRCh37 release 111) emit the UniProt
    *entry name* / mnemonic instead:
      - SwissProt: ``PK3CD_HUMAN`` → resolved via the ``aliases`` table.
      - TrEMBL: ``B7ZM44_HUMAN`` → the accession is the prefix, so stripping
        ``_HUMAN`` recovers it when no alias row exists.
    Falls back to stripping ``_HUMAN`` if the ``aliases`` table is absent
    (older DB build) so the lookup degrades gracefully rather than raising.
    """
    if not uniprot_id.endswith("_HUMAN"):
        return uniprot_id
    try:
        row = conn.execute(
            "SELECT uniprot_id FROM aliases WHERE entry_name = ? LIMIT 1",
            (uniprot_id,),
        ).fetchone()
    except sqlite3.Error:
        row = None
    if row:
        return row[0]
    return uniprot_id[: -len("_HUMAN")]


def query_esm1b(
    sqlite_path: Path,
    uniprot_id: str,
    aa_pos: int,
    alt_aa: str,
) -> Optional[ESM1bData]:
    """Look up ESM1b LLR for a missense variant.

    Returns None if the database is missing, the protein/position is not
    covered, or the lookup raises. A None result causes PP3/BP4 to skip the
    ESM1b branch — never raise upstream.
    """
    if not sqlite_path.exists():
        log.warning("esm1b_missing", path=str(sqlite_path))
        return None
    if not uniprot_id or aa_pos is None or not alt_aa:
        return None

    try:
        # mode=ro opens the database read-only via URI so multiple worker
        # processes can share the same file without write-lock contention.
        # The variant pipeline only reads ESM1b, so RO is both safer and
        # measurably faster under concurrent access.
        conn = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
        try:
            return _select_llr(conn, uniprot_id, aa_pos, alt_aa)
        finally:
            conn.close()
    except sqlite3.Error as exc:
        log.error(
            "esm1b_sqlite_error",
            error=str(exc),
            uniprot=uniprot_id, aa_pos=aa_pos, alt_aa=alt_aa,
        )
        return None


class Esm1bDB:
    """Persistent-connection ESM1b reader for batch annotation.

    Replaces the per-variant ``sqlite3.connect()`` in query_esm1b() with a
    thread-local read-only connection reused across the whole batch (the
    annotate_batch thread pool calls lookup() from several worker threads, and a
    SQLite connection may not cross threads — so each thread opens its own once).
    The accession resolution + point SELECT are identical to query_esm1b via the
    shared _select_llr helper; only the connection lifecycle differs."""

    def __init__(self, sqlite_path: Path) -> None:
        self._path = sqlite_path
        self._local = threading.local()

    def _conn(self) -> Optional[sqlite3.Connection]:
        conn = getattr(self._local, "conn", None)
        if conn is None:
            # RO via URI mirrors query_esm1b: read-only avoids write-lock
            # contention. Kept per-thread (threading.local) so the default
            # check_same_thread guard is satisfied without disabling it.
            conn = sqlite3.connect(f"file:{self._path}?mode=ro", uri=True)
            self._local.conn = conn
        return conn

    def lookup(
        self, uniprot_id: str, aa_pos: int, alt_aa: str
    ) -> Optional[ESM1bData]:
        """Look up ESM1b LLR reusing a thread-local connection.

        Same contract as query_esm1b: returns None when the DB is missing, the
        inputs are incomplete, the protein/position is uncovered, or the lookup
        raises — never propagates upstream."""
        if not self._path.exists():
            log.warning("esm1b_missing", path=str(self._path))
            return None
        if not uniprot_id or aa_pos is None or not alt_aa:
            return None
        try:
            return _select_llr(self._conn(), uniprot_id, aa_pos, alt_aa)
        except sqlite3.Error as exc:
            log.error(
                "esm1b_sqlite_error",
                error=str(exc),
                uniprot=uniprot_id, aa_pos=aa_pos, alt_aa=alt_aa,
            )
            return None
