"""tabix query for REVEL precomputed scores (PP3/BP4)."""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional

import structlog

from acmg_classifier.models.annotation import RevelData
from acmg_classifier.utils.tabix import TabixReader, open_tabix, fetch_region

log = structlog.get_logger()


def _base(tx: str) -> str:
    return tx.split(".", 1)[0].strip()


def _match_revel(
    lines: Iterable[str],
    pos: int,
    ref: str,
    alt: str,
    aa: Optional[tuple[str, str]] = None,
    transcripts: Optional[Iterable[str]] = None,
) -> Optional[RevelData]:
    """Pick the REVEL score from tabix region lines matching (pos, ref, alt).

    REVEL v1.3 carries one row per (pos, ref, alt, aaref, aaalt): the same SNV
    can have different scores on different transcripts (e.g. IDUA c.530T>G is
    F>C = 0.961 on the MANE transcript but S>A = 0.400 on ENST00000504568).
    The 8-column file written by scripts/setup_data.py keeps aaref, aaalt and
    the Ensembl transcript IDs, and the row is chosen as

      1. the row whose transcript list contains a MANE Select ENST
         (``transcripts``),
      2. else the row whose amino-acid change equals the MANE consequence
         (``aa`` = (ref_aa, alt_aa), one-letter),
      3. else the highest score among the matching rows.

    A legacy 5-column file (chrom, pos, ref, alt, REVEL) has no transcript
    information, so its first matching row is returned as before.

    Shared by the connection-per-call query_revel and the persistent-handle
    query_revel_reader so both parse identically."""
    tx_set = {_base(t) for t in (transcripts or ()) if t}
    rows: list[tuple[float, str, str, set[str]]] = []
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        f_pos, f_ref, f_alt = fields[1], fields[2], fields[3]
        if f_ref != ref or f_alt != alt:
            continue
        try:
            if int(f_pos) != pos:
                continue
            score = float(fields[4])
        except ValueError:
            # A single malformed row must not break annotation.
            continue
        if len(fields) < 7:
            return RevelData(score=score, match="legacy")
        txs = {_base(t) for t in (fields[7].split(";") if len(fields) > 7 else []) if t and t != "."}
        rows.append((score, fields[5], fields[6], txs))
    if not rows:
        return None
    if tx_set:
        hit = [r for r in rows if r[3] & tx_set]
        if hit:
            return RevelData(score=max(r[0] for r in hit), match="transcript")
    if aa is not None:
        hit = [r for r in rows if (r[1], r[2]) == aa]
        if hit:
            return RevelData(score=max(r[0] for r in hit), match="amino_acid")
    return RevelData(score=max(r[0] for r in rows), match="max" if len(rows) > 1 else "single")


def query_revel(
    tsv_gz_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    aa: Optional[tuple[str, str]] = None,
    transcripts: Optional[Iterable[str]] = None,
) -> Optional[RevelData]:
    """Query the REVEL TSV (tabix-indexed) for one SNV.

    scripts/setup_data.py normalises REVEL into an 8-column, tab-separated,
    bgzipped file indexed on the assembly's position column:

        chrom  pos  ref  alt  REVEL  aaref  aaalt  Ensembl_transcriptids

    REVEL covers only single-nucleotide missense substitutions, so rows are
    matched on (pos, ref, alt) and the MANE row is chosen by ``_match_revel``. Returns None when the file is absent, the
    variant is uncovered, or the lookup raises — PP3/BP4 then simply skip the
    REVEL branch rather than aborting annotation.
    """
    if not tsv_gz_path.exists():
        log.warning("revel_missing", path=str(tsv_gz_path))
        return None

    try:
        with open_tabix(tsv_gz_path) as tf:
            return _match_revel(fetch_region(tf, chrom, pos, pos), pos, ref, alt,
                                aa, transcripts)
    except Exception as exc:
        log.error("revel_error", error=str(exc))
    return None


def query_revel_reader(
    reader: TabixReader,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    aa: Optional[tuple[str, str]] = None,
    transcripts: Optional[Iterable[str]] = None,
) -> Optional[RevelData]:
    """Batch variant of query_revel using a persistent thread-local handle."""
    if not reader.exists():
        log.warning("revel_missing", path=str(reader.path))
        return None
    try:
        return _match_revel(reader.fetch(chrom, pos, pos), pos, ref, alt, aa, transcripts)
    except Exception as exc:
        log.error("revel_error", error=str(exc))
    return None
