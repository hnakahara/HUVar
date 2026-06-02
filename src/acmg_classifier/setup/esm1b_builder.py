"""Build the ESM1b SQLite from the Brandes 2023 distribution.

Source archive (downloaded by `scripts/setup_data.py`):
  - ALL_hum_isoforms_ESM1b_LLR.zip
      One CSV per protein isoform, file name `<UniProt>_LLR.csv`
      (e.g. `P38398_LLR.csv`, `P38398-2_LLR.csv`).

CSV layout (one isoform):

    ,M 1,A 2,A 3,E 4,L 5,...
    K,-11.210,-6.968,-6.102,-4.795,-4.373,...
    R,-12.401,-5.842,...,...
    ...

  - First row, first cell is empty; remaining header cells are "<WT_AA> <POS>"
    pairs (note the single space). The position is 1-based and the WT amino
    acid is the residue at that position.
  - Each subsequent row's first cell is the *alt* amino acid; remaining cells
    are LLR values per position. A cell is the LLR for substituting the
    column's WT with the row's alt AA. Cells where the alt AA matches the
    column's WT (the diagonal) are 0.000 — skipped.

Output SQLite (no isoform mapping needed):

    CREATE TABLE scores (
        uniprot_id TEXT NOT NULL,
        aa_pos     INTEGER NOT NULL,
        alt_aa     TEXT NOT NULL,
        llr        REAL NOT NULL,
        PRIMARY KEY (uniprot_id, aa_pos, alt_aa)
    ) WITHOUT ROWID;
    CREATE INDEX idx_scores_uni_pos ON scores(uniprot_id, aa_pos);

    -- UniProt entry-name (mnemonic) → accession, built from isoform_list.csv.
    -- Lets the lookup normalise tokens like 'PK3CD_HUMAN' that some VEP caches
    -- emit instead of the accession 'O00329' the scores table is keyed on.
    CREATE TABLE aliases (
        entry_name TEXT PRIMARY KEY,
        uniprot_id TEXT NOT NULL
    ) WITHOUT ROWID;
"""
from __future__ import annotations

import csv
import io
import re
import sqlite3
import zipfile
from pathlib import Path
from typing import Iterator

from acmg_classifier.utils.progress import progress_bar

_HEADER_CELL = re.compile(r"^([A-Z*])\s+(\d+)$")
_MNEMONIC = re.compile(r"\(([^)]+)\)")


def _iter_aliases(zip_path: Path) -> Iterator[tuple[str, str]]:
    """Yield (entry_name, uniprot_id) pairs from the Brandes isoform_list.csv.

    VEP --uniprot on some caches (observed: GRCh37 release 111) reports the
    UniProt *entry name* / mnemonic (e.g. ``PK3CD_HUMAN``) instead of the
    accession (``O00329``) that the score table is keyed on. isoform_list.csv
    maps each accession to a display string of the form
    ``GENE (MNEMONIC) | ACCESSION``, so the parenthesised mnemonic + ``_HUMAN``
    reconstructs the entry name VEP emits, letting the query normalise it back
    to the accession.

    The CSV is looked up next to the zip first (the Brandes download ships it
    alongside the archive), then inside the zip as a fallback.
    """
    text: str | None = None
    sidecar = zip_path.parent / "isoform_list.csv"
    if sidecar.exists():
        text = sidecar.read_text(encoding="utf-8")
    else:
        try:
            with zipfile.ZipFile(zip_path) as zf:
                for name in zf.namelist():
                    if Path(name).name == "isoform_list.csv":
                        with zf.open(name) as fh:
                            text = io.TextIOWrapper(fh, encoding="utf-8").read()
                        break
        except (zipfile.BadZipFile, OSError):
            return
    if not text:
        return

    reader = csv.reader(io.StringIO(text))
    next(reader, None)  # header: id,txt
    for row in reader:
        if len(row) < 2:
            continue
        accession = row[0].strip()
        if not accession:
            continue
        # The display column can itself contain commas — rejoin the remainder.
        display = ",".join(row[1:]).split("|", 1)[0]
        mnemonics = _MNEMONIC.findall(display)
        if not mnemonics:
            continue
        entry_name = mnemonics[-1].strip() + "_HUMAN"
        if entry_name != "_HUMAN":
            yield entry_name, accession


def _parse_header(header: list[str]) -> list[tuple[str, int] | None]:
    """Map each header column index (excluding the empty leading cell) to
    a (wt_aa, position) tuple. Returns None for cells that do not match
    the "<WT_AA> <POS>" pattern so the corresponding data columns can be
    skipped without aborting the whole file.
    """
    parsed: list[tuple[str, int] | None] = []
    for cell in header[1:]:
        m = _HEADER_CELL.match(cell.strip())
        if m:
            parsed.append((m.group(1), int(m.group(2))))
        else:
            parsed.append(None)
    return parsed


def _iter_matrix_rows(
    uniprot_id: str,
    text: str,
) -> Iterator[tuple[str, int, str, float]]:
    """Yield (uniprot_id, aa_pos, alt_aa, llr) from one isoform CSV."""
    reader = csv.reader(io.StringIO(text))
    header = next(reader, None)
    if not header or len(header) < 2:
        return
    pos_info = _parse_header(header)

    for row in reader:
        if not row:
            continue
        alt_aa = row[0].strip()
        if not alt_aa or len(alt_aa) != 1:
            continue
        for col_idx, cell in enumerate(row[1:]):
            if col_idx >= len(pos_info):
                break
            info = pos_info[col_idx]
            if info is None:
                continue
            wt_aa, pos = info
            if alt_aa == wt_aa:
                # WT-to-self diagonal — Brandes encodes these as 0.000.
                continue
            cell = cell.strip()
            if not cell:
                continue
            try:
                llr = float(cell)
            except ValueError:
                continue
            yield uniprot_id, pos, alt_aa, llr


def _uniprot_from_name(name: str) -> str | None:
    """`P38398_LLR.csv` or `subdir/P38398-2_LLR.csv` → `P38398` / `P38398-2`."""
    stem = Path(name).stem  # `P38398_LLR`
    if not stem.endswith("_LLR"):
        return None
    uni = stem[: -len("_LLR")]
    return uni or None


def build_esm1b_sqlite(
    zip_path: Path,
    dest: Path,
    *,
    batch_size: int = 100_000,
) -> None:
    """Build the ESM1b SQLite from the Brandes 2023 zip."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        dest.unlink()

    conn = sqlite3.connect(str(dest))
    try:
        conn.execute("PRAGMA journal_mode = OFF")
        conn.execute("PRAGMA synchronous = OFF")
        conn.execute(
            "CREATE TABLE scores ("
            " uniprot_id TEXT NOT NULL,"
            " aa_pos INTEGER NOT NULL,"
            " alt_aa TEXT NOT NULL,"
            " llr REAL NOT NULL,"
            " PRIMARY KEY (uniprot_id, aa_pos, alt_aa)"
            ") WITHOUT ROWID"
        )
        # Entry-name → accession map so the lookup can normalise UniProt
        # mnemonics (e.g. PK3CD_HUMAN) that some VEP caches emit instead of
        # the accession (O00329) the score table is keyed on.
        conn.execute(
            "CREATE TABLE aliases ("
            " entry_name TEXT PRIMARY KEY,"
            " uniprot_id TEXT NOT NULL"
            ") WITHOUT ROWID"
        )
        alias_rows = list(_iter_aliases(zip_path))
        if alias_rows:
            conn.executemany(
                "INSERT OR IGNORE INTO aliases VALUES (?, ?)", alias_rows
            )

        buf: list[tuple[str, int, str, float]] = []
        n_isoforms = 0
        n_rows = 0
        with zipfile.ZipFile(zip_path) as zf:
            # Pre-count CSV entries so the progress bar has an accurate total
            # — the namelist is already in memory once the zip is open, so
            # this is essentially free.
            csv_names = [n for n in zf.namelist() if n.endswith(".csv")]
            with progress_bar("Building ESM1b SQLite", total=len(csv_names)) as advance:
                for name in csv_names:
                    uni_id = _uniprot_from_name(name)
                    if uni_id is None:
                        advance()
                        continue
                    with zf.open(name) as fh:
                        text = io.TextIOWrapper(fh, encoding="utf-8").read()
                    for tup in _iter_matrix_rows(uni_id, text):
                        buf.append(tup)
                        if len(buf) >= batch_size:
                            conn.executemany(
                                "INSERT OR IGNORE INTO scores VALUES (?, ?, ?, ?)",
                                buf,
                            )
                            n_rows += len(buf)
                            buf.clear()
                    n_isoforms += 1
                    advance()
        if buf:
            conn.executemany(
                "INSERT OR IGNORE INTO scores VALUES (?, ?, ?, ?)",
                buf,
            )
            n_rows += len(buf)

        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_scores_uni_pos "
            "ON scores(uniprot_id, aa_pos)"
        )
        conn.commit()
        print(
            f"  ESM1b: {n_isoforms} isoforms, {n_rows:,} rows, "
            f"{len(alias_rows):,} aliases → {dest.name}"
        )
    finally:
        conn.close()
