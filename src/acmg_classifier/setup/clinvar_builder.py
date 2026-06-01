"""Build ClinVar SQLite from XML for PS1/PM5/PS3/PP1/PS4 lookups."""
from __future__ import annotations
import gzip
import multiprocessing as mp
import re
import sqlite3
import xml.etree.ElementTree as ET
from collections.abc import Iterator
from pathlib import Path

import structlog

from acmg_classifier.utils.progress import progress_bar

log = structlog.get_logger()

# Each parse task carries this many <ClinVarSet> records. Batching keeps the
# inter-process pickling/IPC overhead negligible relative to the parse+regex
# cost that the workers actually do.
_BATCH_SIZE = 200

_CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS variants (
    variation_id TEXT,
    chrom TEXT,
    pos INTEGER,
    ref TEXT,
    alt TEXT,
    gene_symbol TEXT,
    hgvs_c TEXT,
    hgvs_p TEXT,
    amino_acid_change TEXT,
    codon_position INTEGER,
    clinical_significance TEXT,
    review_status TEXT,
    star_rating INTEGER,
    last_evaluated TEXT,
    affected_cases INTEGER DEFAULT 0,
    functional_evidence INTEGER DEFAULT 0,
    segregation_evidence INTEGER DEFAULT 0
);
"""

# Indexes are created AFTER the bulk load. Building them once at the end is far
# faster than maintaining three indexes on every INSERT during the load.
_CREATE_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_gene_codon ON variants (gene_symbol, codon_position);
CREATE INDEX IF NOT EXISTS idx_chrom_pos ON variants (chrom, pos, ref, alt);
CREATE INDEX IF NOT EXISTS idx_gene_aa ON variants (gene_symbol, amino_acid_change);
"""

# --- SCV free-text evidence mining (PS3 functional / PP1 cosegregation) ---
# Positive indicators of a damaging functional assay (PS3).
_FUNCTIONAL_POS = re.compile(
    r"functional stud|enzymatic activit|enzyme activit|experimental stud|"
    r"in vitro|reduced .{0,30}activit|decreased .{0,30}activit|undetectable|"
    r"abolish|loss of function|affects? .{0,20}function|\bPS3\b",
    re.IGNORECASE,
)
# Phrases that negate a damaging functional result — skip the SCV if present.
# Robust to adverbs between the negator and verb ("does not substantially affect")
# and to both word orders ("normal protein function" / "protein function was normal").
_NEG_ADJ = r"(normal|unaffected|intact|unchanged|preserved|retained|comparable|similar|wild[- ]?type)"
_FUNCTIONAL_NEG = re.compile(
    rf"{_NEG_ADJ}\s+(enzyme |enzymatic |protein |splic\w* )?(activit|function|level|express)|"
    rf"(activit\w*|function|express\w*|level)\s+(\w+\s+){{0,2}}"
        rf"(was|were|is|are|appear\w*|remain\w*|seem\w*)\s+(\w+\s+){{0,1}}{_NEG_ADJ}|"
    r"\bno\b[^.]{0,40}?(effect|impact|damaging|deleterious|functional consequence|"
        r"change in (activit|function|express))|"
    r"(does|did|do|was|were|is|are|has|have|could|would|appear\w*)\s*n[o']?t\b[^.]{0,45}?"
        r"(affect|alter|impair|reduc\w+|abolish|disrupt|damag\w+|chang\w+|impact|influence)|"
    r"\bnot\b[^.]{0,30}?(affect|damaging|deleterious|pathogenic|functional)|"
    r"without[^.]{0,30}?(affect|alter|impair|effect|chang)|"
    r"functionally (neutral|benign|tolerated|silent)|\btolerated\b|"
    r"no effect on (protein |enzyme |splic\w* )?function|"
    # In-silico / computational predictions are NOT functional assays. PS3
    # requires wet-lab evidence, so exclude SCVs whose damaging signal comes
    # from prediction tools (e.g. "tools predict the variant abolishes a
    # splicing donor site"), and those explicitly stating the effect is
    # unconfirmed by functional studies or lacks experimental evidence.
    r"\bin[\s-]?silico\b|"
    r"computational (tool|method|analys|predict|algorithm|approach|program|software)|"
    r"predict\w*[^.]{0,60}?(abolish|disrupt|damag|impact|affect|effect|splic|donor|acceptor|function)|"
    r"yet to be (confirmed|validated|established|proven)|"
    r"no (experimental|functional)[^.]{0,40}?(evidence|stud|data|assay)",
    re.IGNORECASE,
)
# Positive indicators of cosegregation (PP1).
_SEGREGATION_POS = re.compile(r"segregat|\bPP1\b", re.IGNORECASE)
# Phrases that negate cosegregation.
_SEGREGATION_NEG = re.compile(
    r"(did not|does not|didn't|doesn't|failed to|no|lack of|without)\s+\w*\s*segregat",
    re.IGNORECASE,
)


def _star_rating(review_status: str) -> int:
    rs = review_status.lower()
    if "practice guideline" in rs:
        return 4
    if "reviewed by expert panel" in rs:
        return 3
    if "multiple submitters, no conflicts" in rs:
        return 2
    if "single submitter" in rs:
        return 1
    return 0


def _parse_aa_change(hgvs_p: str | None) -> tuple[str | None, int | None]:
    """Extract amino_acid_change ('R175H') and codon_position (175) from HGVS p."""
    import re
    if not hgvs_p:
        return None, None
    m = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)", hgvs_p)
    if m:
        aa1 = m.group(1)[:1].upper()
        pos = int(m.group(2))
        aa2 = m.group(3)[:1].upper()
        return aa1 + str(pos) + aa2, pos
    return None, None


def _iter_clinvarset_chunks(fh, read_size: int = 1 << 20) -> Iterator[str]:
    """Yield each ``<ClinVarSet>...</ClinVarSet>`` block as a self-contained string.

    This is a cheap text scan — no XML parsing — so the (single-threaded,
    gzip-bound) main process spends almost nothing here and the expensive
    ElementTree parse + regex mining happens in the worker pool. Each
    ClinVarSet element is independent, so a sliced block parses standalone
    via ``ET.fromstring``.
    """
    buf = ""
    start_tag = "<ClinVarSet"
    end_tag = "</ClinVarSet>"
    while True:
        data = fh.read(read_size)
        if not data:
            break
        buf += data
        while True:
            end = buf.find(end_tag)
            if end == -1:
                break  # need more data for a complete block
            end += len(end_tag)
            start = buf.find(start_tag)
            if start == -1 or start > end:
                # Stray closing tag with no matching open in the buffer — drop
                # the consumed prefix and keep scanning. (Should not happen on
                # well-formed ClinVar XML.)
                buf = buf[end:]
                continue
            yield buf[start:end]
            buf = buf[end:]  # trim so the buffer stays ~one record long


def _batched(it: Iterator[str], n: int) -> Iterator[list[str]]:
    batch: list[str] = []
    for item in it:
        batch.append(item)
        if len(batch) >= n:
            yield batch
            batch = []
    if batch:
        yield batch


def _parse_chunk_batch(args: tuple[list[str], str]) -> tuple[int, list[tuple]]:
    """Worker: parse a batch of ClinVarSet strings into row tuples.

    Returns ``(n_chunks, rows)`` — n_chunks drives the progress bar (rows is
    already filtered, so it under-counts processed records).
    """
    chunks, assembly = args
    rows: list[tuple] = []
    for chunk in chunks:
        try:
            elem = ET.fromstring(chunk)
        except ET.ParseError:
            continue
        row = _parse_clinvarset(elem, assembly)
        if row:
            rows.append(row)
    return len(chunks), rows


def build_clinvar_sqlite(
    xml_gz_path: Path,
    output_db: Path,
    assembly: str,
    workers: int | None = None,
) -> None:
    """Parse ClinVar XML and populate SQLite for PS1/PM5 queries.

    Parsing is parallelized across a process pool: the main process streams the
    gzip and splits it into per-record strings (cheap), workers do the costly
    ElementTree parse + free-text regex mining, and the main process inserts the
    returned rows. This keeps all cores busy instead of bottlenecking on a
    single-threaded ``iterparse``.

    ``workers`` sets the number of parse processes. ``None`` (default) uses 4.
    Values are clamped to the range 1..24.
    """
    log.info("building_clinvar_sqlite", output=str(output_db))
    output_db.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(output_db))
    # Bulk-load PRAGMAs. This DB is built from scratch, so durability is not needed
    # (just rebuild on failure); these cut the per-commit fsync / journal overhead.
    con.execute("PRAGMA journal_mode = OFF")
    con.execute("PRAGMA synchronous = OFF")
    con.execute("PRAGMA temp_store = MEMORY")
    con.execute("PRAGMA cache_size = -262144")  # ~256 MB page cache
    con.executescript(_CREATE_TABLE)  # table only; indexes are built after the load

    insert_sql = "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    n_rows = 0
    n_seen = 0
    # Default to 4 parse processes; cap at 24 so the main process (gzip
    # decompression + SQLite inserts) isn't starved and IPC stays manageable.
    n_workers = min(24, max(1, workers if workers is not None else 4))
    # Read decompressed text directly; ClinVar XML is declared UTF-8.
    if str(xml_gz_path).endswith(".gz"):
        fh = gzip.open(str(xml_gz_path), "rt", encoding="utf-8")
    else:
        fh = open(str(xml_gz_path), "rt", encoding="utf-8")

    # ClinVar XML doesn't expose a record count, so we run an indeterminate
    # progress bar (total=None) — rich shows a spinner + the running tick
    # count, which is sufficient feedback for a multi-minute build.
    rows: list[tuple] = []
    with fh, progress_bar("Parsing ClinVar XML", total=None) as advance, \
            mp.Pool(processes=n_workers) as pool:  # type: ignore[arg-type]
        batches = ((b, assembly) for b in _batched(_iter_clinvarset_chunks(fh), _BATCH_SIZE))
        # imap_unordered: results stream back as workers finish, overlapping
        # parsing with the main thread's SQLite inserts.
        for n_chunks, batch_rows in pool.imap_unordered(_parse_chunk_batch, batches):
            n_seen += n_chunks
            advance(n_chunks)
            if batch_rows:
                rows.extend(batch_rows)
                if len(rows) >= 5000:
                    con.executemany(insert_sql, rows)
                    n_rows += len(rows)
                    rows = []

    if rows:
        con.executemany(insert_sql, rows)
        n_rows += len(rows)
    con.commit()

    # Build indexes once, on the fully loaded table.
    log.info("clinvar_building_indexes", rows=n_rows, seen=n_seen)
    con.executescript(_CREATE_INDEXES)
    con.commit()
    con.close()
    log.info("clinvar_sqlite_built", path=str(output_db), rows=n_rows)


def _parse_clinvarset(elem: ET.Element, assembly: str):
    try:
        ref_assert = elem.find(".//ReferenceClinVarAssertion")
        if ref_assert is None:
            return None

        var_id = ref_assert.get("ID", "")
        clinsig_elem = ref_assert.find(".//ClinicalSignificance/Description")
        clinsig = clinsig_elem.text if clinsig_elem is not None else ""

        rev_elem = ref_assert.find(".//ClinicalSignificance/ReviewStatus")
        rev_status = rev_elem.text if rev_elem is not None else ""
        stars = _star_rating(rev_status)

        last_elem = ref_assert.find(".//ClinicalSignificance")
        last_eval = last_elem.get("DateLastEvaluated", "") if last_elem is not None else ""

        gene_elem = ref_assert.find(".//MeasureSet/Measure/MeasureRelationship/Symbol/ElementValue")
        gene = gene_elem.text if gene_elem is not None else ""

        hgvs_c, hgvs_p = None, None
        for attr in ref_assert.findall(".//MeasureSet/Measure/AttributeSet/Attribute"):
            atype = attr.get("Type", "")
            if atype == "HGVS, coding, RefSeq" and attr.text:
                hgvs_c = attr.text
            if atype == "HGVS, protein, RefSeq" and attr.text:
                hgvs_p = attr.text

        aa_change, codon_pos = _parse_aa_change(hgvs_p)

        loc = ref_assert.find(
            f".//MeasureSet/Measure/SequenceLocation[@Assembly='{assembly}']"
        )
        chrom = loc.get("Chr") if loc is not None else None
        pos = int(loc.get("positionVCF", 0)) if loc is not None else None
        ref = loc.get("referenceAlleleVCF") if loc is not None else None
        alt = loc.get("alternateAlleleVCF") if loc is not None else None

        # Per-SCV evidence mining. Each <ClinVarAssertion> is one SCV submission.
        #   affected     -> PS4 (AffectedStatus="yes")
        #   functional   -> PS3 (free-text comment describing a damaging assay)
        #   segregation  -> PP1 (free-text comment describing cosegregation)
        affected = 0
        functional = 0
        segregation = 0
        for scv in elem.findall(".//ClinVarAssertion"):
            # PS4 counts affected probands ONLY from P/LP submissions. An affected
            # individual reported by a Benign/VUS submitter is an incidental finding,
            # not evidence that the variant causes disease, so it must not be counted.
            scv_sig_elem = scv.find(".//ClinicalSignificance/Description")
            scv_sig = (scv_sig_elem.text or "").strip().lower() if scv_sig_elem is not None else ""
            scv_is_plp = scv_sig.startswith("pathogenic") or scv_sig.startswith("likely pathogenic")
            if scv_is_plp:
                for status_elem in scv.findall(".//ObservedIn/Sample/AffectedStatus"):
                    if (status_elem.text or "").strip().lower() == "yes":
                        affected += 1
                        break  # one affected observation per SCV is enough

            texts: list[str] = []
            for comment in scv.findall(".//Comment"):
                if comment.text:
                    texts.append(comment.text)
            for desc in scv.findall(".//ObservedData/Attribute[@Type='Description']"):
                if desc.text:
                    texts.append(desc.text)
            combined = " ".join(texts)
            if combined:
                if _FUNCTIONAL_POS.search(combined) and not _FUNCTIONAL_NEG.search(combined):
                    functional += 1
                if _SEGREGATION_POS.search(combined) and not _SEGREGATION_NEG.search(combined):
                    segregation += 1

        return (var_id, chrom, pos, ref, alt, gene, hgvs_c, hgvs_p,
                aa_change, codon_pos, clinsig, rev_status, stars, last_eval,
                affected, functional, segregation)
    except Exception:
        return None
