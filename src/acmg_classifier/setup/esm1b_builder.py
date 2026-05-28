"""Build the ESM1b SQLite from Brandes 2023 distribution.

Inputs (downloaded by `scripts/setup_data.py`):
  - ALL_hum_isoforms_ESM1b_LLR.zip
      One file per isoform inside, keyed by UniProt ID. Each file is a CSV
      matrix where rows are amino-acid positions (1-based) and columns are
      the 20 standard amino acids; cells are the LLR for the substitution.
      The WT position (diagonal) contains 0 or NaN and is skipped.
  - isoform_list.csv
      UniProt-to-Ensembl-transcript mapping with columns
        (uniprot_id, ensembl_transcript_id, ...)
      One UniProt isoform maps to exactly one Ensembl transcript.

Output:
  - SQLite at the supplied path with schema
      scores(transcript_id TEXT, aa_pos INTEGER, alt_aa TEXT, llr REAL)
      indexed on (transcript_id, aa_pos).

Note: the upstream archive ships a matrix-style CSV. If a future release
ships long-form (mutation,score) instead, set `long_form=True` to parse it.
"""
from __future__ import annotations

import csv
import io
import sqlite3
import zipfile
from pathlib import Path
from typing import Iterator

_AA_COLS = list("ACDEFGHIKLMNPQRSTVWY")


def _build_isoform_map(isoform_csv: Path) -> dict[str, str]:
    """Return {uniprot_id (stripped of -N suffix): ensembl_transcript_id}."""
    mapping: dict[str, str] = {}
    with isoform_csv.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            uni = (row.get("uniprot_id") or row.get("uniprot") or "").strip()
            ens = (
                row.get("ensembl_transcript_id")
                or row.get("transcript_id")
                or row.get("ensembl")
                or ""
            ).strip()
            if not uni or not ens:
                continue
            mapping[uni] = ens.split(".", 1)[0]
    return mapping


def _iter_matrix_rows(
    transcript_id: str,
    text: str,
) -> Iterator[tuple[str, int, str, float]]:
    """Yield (transcript_id, aa_pos, alt_aa, llr) from a matrix-style CSV."""
    reader = csv.reader(io.StringIO(text))
    header = next(reader, None)
    if not header:
        return
    # header[0] is the position column; remaining cells are AA labels.
    aa_labels = [h.strip() for h in header[1:]]
    for row in reader:
        if not row or len(row) < 2:
            continue
        try:
            pos = int(row[0])
        except ValueError:
            continue
        for label, cell in zip(aa_labels, row[1:]):
            if not label or not cell:
                continue
            try:
                llr = float(cell)
            except ValueError:
                continue
            if llr == 0.0:
                # WT diagonal — Brandes archives encode WT cells as 0.
                continue
            yield transcript_id, pos, label, llr


def build_esm1b_sqlite(
    zip_path: Path,
    isoform_csv: Path,
    dest: Path,
    *,
    batch_size: int = 50_000,
) -> None:
    """Build the ESM1b SQLite from the Brandes 2023 zip + isoform map."""
    uni2tx = _build_isoform_map(isoform_csv)
    if not uni2tx:
        raise RuntimeError(f"Empty UniProt→ENST map from {isoform_csv}")

    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        dest.unlink()

    conn = sqlite3.connect(str(dest))
    try:
        conn.execute("PRAGMA journal_mode = OFF")
        conn.execute("PRAGMA synchronous = OFF")
        conn.execute(
            "CREATE TABLE scores ("
            " transcript_id TEXT NOT NULL,"
            " aa_pos INTEGER NOT NULL,"
            " alt_aa TEXT NOT NULL,"
            " llr REAL NOT NULL,"
            " PRIMARY KEY (transcript_id, aa_pos, alt_aa)"
            ") WITHOUT ROWID"
        )

        buf: list[tuple[str, int, str, float]] = []
        n_isoforms = 0
        n_rows = 0
        with zipfile.ZipFile(zip_path) as zf:
            for name in zf.namelist():
                if not name.endswith(".csv"):
                    continue
                # File name e.g. "Q14524_LLR.csv" or "Q14524-2_LLR.csv".
                stem = Path(name).stem  # "Q14524_LLR"
                uni_id = stem.replace("_LLR", "").strip()
                tx = uni2tx.get(uni_id)
                if tx is None:
                    continue
                with zf.open(name) as fh:
                    text = io.TextIOWrapper(fh, encoding="utf-8").read()
                for tup in _iter_matrix_rows(tx, text):
                    buf.append(tup)
                    if len(buf) >= batch_size:
                        conn.executemany(
                            "INSERT OR IGNORE INTO scores VALUES (?, ?, ?, ?)",
                            buf,
                        )
                        n_rows += len(buf)
                        buf.clear()
                n_isoforms += 1
        if buf:
            conn.executemany(
                "INSERT OR IGNORE INTO scores VALUES (?, ?, ?, ?)",
                buf,
            )
            n_rows += len(buf)

        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_scores_tx_pos "
            "ON scores(transcript_id, aa_pos)"
        )
        conn.commit()
        print(f"  ESM1b: {n_isoforms} isoforms, {n_rows:,} rows → {dest.name}")
    finally:
        conn.close()
