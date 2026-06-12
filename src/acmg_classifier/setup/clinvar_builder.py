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

from acmg_classifier.criteria.grantham import AA3_TO_AA1
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
    segregation_evidence INTEGER DEFAULT 0,
    bs2_evidence INTEGER DEFAULT 0,
    bs2_strength TEXT
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

# --- BS2 mining (expert-panel benign "observed in healthy adult") ---
# Several VCEPs (CDH1, TP53, SERPINC1, ...) apply BS2 from *internal* laboratory
# cohorts that gnomAD cannot supply, so our gnomAD-based BS2 is forced
# not_applicable for these genes. When an expert panel (>=3 stars) has already
# applied BS2 to the specific variant, its review cites the criterion code
# explicitly ("...without a diagnosis... (BS2)"), so we can harvest that
# variant-level judgement. Because the code is named verbatim this is far more
# reliable than the PS3/PP1 free-text heuristics above.
#
# The strength modifier ("BS2_Supporting" / "BS2_Moderate" / "BS2_Strong" /
# "BS2_VeryStrong") is captured when present; a bare "BS2" applies at its ACMG
# Strong default.
_BS2_POS = re.compile(
    r"\bBS2(?:_?(VeryStrong|Very[\s_]Strong|Strong|Moderate|Supporting))?\b",
    re.IGNORECASE,
)
# Skip when BS2 is cited as *not* applied/applicable (so a "BS2 not applicable"
# statement never reads as positive BS2 evidence). Targets BS2 specifically so a
# neighbouring "BA1 not applicable" does not suppress a genuine BS2.
_BS2_NEG = re.compile(
    r"BS2[\s_]*(?:is\s+|was\s+|were\s+|are\s+)?(?:not\s+(?:applicable|met|applied)|n/?a\b)|"
    r"(?:not\s+applicable|does\s+not\s+meet|did\s+not\s+meet|cannot\s+(?:be\s+)?"
    r"(?:applied|met)|unable\s+to\s+(?:apply|meet))[^.]{0,30}?BS2",
    re.IGNORECASE,
)

# Benign-direction strength rank, strongest wins when several BS2 mentions differ.
_BS2_STRENGTH_RANK = {"Supporting": 1, "Moderate": 2, "Strong": 3, "VeryStrong": 4}


def _norm_bs2_strength(modifier: str | None) -> str:
    """Canonical strength label from a captured BS2 modifier; bare BS2 -> Strong
    (its ACMG default)."""
    if not modifier:
        return "Strong"
    key = re.sub(r"[\s_]", "", modifier).lower()
    return {
        "verystrong": "VeryStrong",
        "strong": "Strong",
        "moderate": "Moderate",
        "supporting": "Supporting",
    }.get(key, "Strong")


def _mine_bs2(text: str) -> tuple[int, str | None]:
    """(hit, strength) for an expert-panel BS2 citation in *text*.

    Returns ``(0, None)`` when BS2 is absent or cited as not-applicable; else
    ``(1, strength)`` with the strongest cited strength ("Strong" for a bare
    "BS2"). The caller gates on star_rating >= 3, so this is only mined from
    expert-panel / practice-guideline reviews."""
    if not text or _BS2_NEG.search(text):
        return 0, None
    best: str | None = None
    best_rank = -1
    for m in _BS2_POS.finditer(text):
        strength = _norm_bs2_strength(m.group(1))
        rank = _BS2_STRENGTH_RANK[strength]
        if rank > best_rank:
            best, best_rank = strength, rank
    if best is None:
        return 0, None
    return 1, best


# Non-coding / uncharacterised locus prefixes. When a variant overlaps such a
# locus AND a real gene, the real gene is the functionally relevant one.
_NONCODING_GENE_PREFIX = ("LOC", "LINC", "MIR", "LNC", "SNORD", "SNORA")


def _gene_symbol(ref_assert: ET.Element) -> str:
    """Best gene symbol for the variant, robust to ClinVar overlapping loci.

    ClinVar tags each gene with a ``MeasureRelationship Type``: ``variant in
    gene`` / ``within single gene`` is the functional gene, ``within multiple
    genes by overlap`` is an incidental neighbour (often a LOC/LINC locus). The
    legacy code took the FIRST ``Symbol`` in document order, so for a gene
    overlapped by such a locus (e.g. PAH vs LOC126861615) it picked the wrong
    symbol and silently broke PS1/PM5 gene matching (the comparator was stored
    under LOC…, so ``WHERE gene_symbol='PAH'`` found nothing). Prefer the
    functional relationship, then a non-LOC symbol."""
    best = ""
    best_rank = (99, 99)
    for rel in ref_assert.findall(".//MeasureSet/Measure/MeasureRelationship"):
        sym = rel.find("./Symbol/ElementValue[@Type='Preferred']")
        if sym is None:
            sym = rel.find("./Symbol/ElementValue")
        if sym is None or not sym.text:
            continue
        name = sym.text.strip()
        t = (rel.get("Type") or "").lower()
        type_score = (
            0 if ("single gene" in t or "variant in gene" in t)
            else 2 if "overlap" in t
            else 1
        )
        loc_penalty = 1 if name.upper().startswith(_NONCODING_GENE_PREFIX) else 0
        rank = (type_score, loc_penalty)
        if rank < best_rank:
            best, best_rank = name, rank
    return best


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


def _rcv_classification(ref_assert: ET.Element) -> tuple[str, str, str]:
    """(significance, review_status, date_last_evaluated) for an RCV record.

    Supports both ClinVar XML schemas. The new RCV release (ClinVar_RCV_2.3.xsd,
    ``RCV_release/ClinVarRCVRelease_*``) nests the aggregate germline call under
    ``Classifications/GermlineClassification`` (ReviewStatus + Description, the
    latter carrying DateLastEvaluated); the legacy ``RCV_xml_old_format`` release
    used a flat ``ClinicalSignificance`` element. We try the new path first and
    fall back to the legacy one so either download builds correctly."""
    gc = ref_assert.find(".//Classifications/GermlineClassification")
    if gc is not None:
        desc = gc.find("Description")
        rev = gc.find("ReviewStatus")
        sig = desc.text if desc is not None and desc.text else ""
        rs = rev.text if rev is not None and rev.text else ""
        dle = desc.get("DateLastEvaluated", "") if desc is not None else ""
        return sig, rs, dle
    cs = ref_assert.find(".//ClinicalSignificance")
    if cs is not None:
        desc = cs.find("Description")
        rev = cs.find("ReviewStatus")
        sig = desc.text if desc is not None and desc.text else ""
        rs = rev.text if rev is not None and rev.text else ""
        return sig, rs, cs.get("DateLastEvaluated", "")
    return "", "", ""


def _scv_significance(scv: ET.Element) -> str:
    """Lower-cased SCV germline significance, both schemas.

    New schema: ``Classification/GermlineClassification`` is the significance
    text itself (a sibling of ReviewStatus). Legacy schema:
    ``ClinicalSignificance/Description``."""
    gc = scv.find(".//Classification/GermlineClassification")
    if gc is not None and gc.text:
        return gc.text.strip().lower()
    d = scv.find(".//ClinicalSignificance/Description")
    return (d.text or "").strip().lower() if d is not None else ""


def _parse_aa_change(hgvs_p: str | None) -> tuple[str | None, int | None]:
    """Extract amino_acid_change ('R175H') and codon_position (175) from HGVS p.

    Tolerates ClinVar's predicted-protein parenthesis notation
    ``NP_000537.3:p.(Arg248Gln)`` — the leading "(" otherwise made the regex
    miss the change, leaving codon_position NULL and silently breaking the
    PS1/PM5 same-codon lookup for every variant stored that way.

    The three-letter residue codes are mapped with the full AA3_TO_AA1 table.
    The previous ``code[:1]`` shortcut collapsed every residue to its first
    letter, so Arg/Asn/Asp all became "A" and Gln/Glu both "G" (likewise
    Lys→L, Phe→P, Trp/Tyr→T). That corrupted amino_acid_change for those
    residues and made the PS1 same-AA lookup match the wrong substitutions."""
    if not hgvs_p:
        return None, None
    m = re.search(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)", hgvs_p)
    if not m:
        return None, None
    pos = int(m.group(2))
    aa1 = AA3_TO_AA1.get(m.group(1))
    raw2 = m.group(3)
    aa2 = "*" if raw2 in ("Ter", "*") else AA3_TO_AA1.get(raw2)
    if aa1 is None or aa2 is None:
        # Unknown residue code (e.g. non-standard 'Sec'/'Xaa') — keep the codon
        # number for the same-codon lookup but leave amino_acid_change NULL.
        return None, pos
    return aa1 + str(pos) + aa2, pos


# A bare protein change ("p.Arg248Trp"), tolerant of a transcript/accession
# prefix and ClinVar's predicted-protein parentheses ("...:p.(Arg248Trp)").
_PROT_CHANGE_RE = re.compile(
    r"p\.\(?([A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|Ter|\*|=))\)?"
)
_PREF_CODING_RE = re.compile(r"(N[MR]_\d+\.\d+)")


def _norm_protein(text: str | None) -> str | None:
    """Canonical bare protein change ('p.Arg248Trp') from any HGVS-ish string,
    stripping the accession prefix and predicted-protein parentheses so the
    Preferred-name change and a per-transcript attribute compare equal."""
    if not text:
        return None
    m = _PROT_CHANGE_RE.search(text)
    return "p." + m.group(1) if m else None


def _preferred_protein_change(ref_assert: ET.Element) -> str | None:
    """The protein change from ClinVar's Preferred variant name, which uses the
    canonical (MANE) transcript — e.g. ``NM_000546.6(TP53):c.742C>T
    (p.Arg248Trp)`` -> ``p.Arg248Trp``. This is the numbering VEP's MANE Select
    annotation uses, so anchoring on it keeps codon_position consistent between
    the candidate and the stored comparators."""
    for path in (
        ".//MeasureSet/Measure/Name/ElementValue[@Type='Preferred']",
        ".//MeasureSet/Name/ElementValue[@Type='Preferred']",
    ):
        for nm in ref_assert.findall(path):
            change = _norm_protein(nm.text)
            if change:
                return change
    return None


def _preferred_coding_accession(ref_assert: ET.Element) -> str | None:
    """The canonical coding transcript accession (e.g. ``NM_000546.6``) from the
    Preferred name, used to pair hgvs_c to the same transcript as hgvs_p."""
    for path in (
        ".//MeasureSet/Measure/Name/ElementValue[@Type='Preferred']",
        ".//MeasureSet/Name/ElementValue[@Type='Preferred']",
    ):
        for nm in ref_assert.findall(path):
            if nm.text:
                m = _PREF_CODING_RE.search(nm.text)
                if m:
                    return m.group(1)
    return None


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

    insert_sql = "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
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
        # Aggregate germline classification — schema-agnostic (new RCV release
        # nests it under Classifications/GermlineClassification, legacy used a
        # flat ClinicalSignificance). See _rcv_classification.
        clinsig, rev_status, last_eval = _rcv_classification(ref_assert)
        stars = _star_rating(rev_status)

        # Functional gene symbol, preferring 'variant in gene' over an
        # incidental 'within multiple genes by overlap' locus (see _gene_symbol).
        gene = _gene_symbol(ref_assert)

        # ClinVar lists a coding+protein HGVS for EVERY transcript of the gene.
        # The old code kept the last-seen protein attribute — an arbitrary
        # transcript — so TP53 p.Arg248Trp was stored under the NP_001263628.1
        # isoform as p.Arg89Trp (codon 89). VEP annotates the candidate on MANE
        # Select (p.Arg248Trp, codon 248), so codon_position disagreed and the
        # PS1/PM5 same-codon lookup silently found no comparator. Anchor on
        # ClinVar's Preferred name (the canonical/MANE transcript) and pick the
        # attribute that matches it, so the stored numbering aligns with VEP.
        hgvs_c, hgvs_p = None, None
        prot_attrs: list[str] = []
        cod_attrs: list[str] = []
        for attr in ref_assert.findall(".//MeasureSet/Measure/AttributeSet/Attribute"):
            atype = attr.get("Type", "")
            if atype == "HGVS, coding, RefSeq" and attr.text:
                hgvs_c = attr.text  # legacy last-seen fallback
                cod_attrs.append(attr.text)
            if atype == "HGVS, protein, RefSeq" and attr.text:
                hgvs_p = attr.text  # legacy last-seen fallback
                prot_attrs.append(attr.text)

        pref_change = _preferred_protein_change(ref_assert)
        if pref_change:
            canon_p = next(
                (p for p in prot_attrs if _norm_protein(p) == pref_change), None
            )
            if canon_p:
                hgvs_p = canon_p
            aa_change, codon_pos = _parse_aa_change(canon_p or pref_change)
            # Pair hgvs_c to the same canonical transcript when identifiable.
            pref_nm = _preferred_coding_accession(ref_assert)
            if pref_nm:
                canon_c = next(
                    (c for c in cod_attrs if c.startswith(pref_nm)), None
                )
                if canon_c:
                    hgvs_c = canon_c
        else:
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
            scv_sig = _scv_significance(scv)
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

        # BS2 from an expert-panel (>=3 star) review: harvest the criterion code
        # the VCEP cited verbatim for this variant. Scanned across the WHOLE
        # ClinVarSet (RCV + SCV comments/descriptions) since the interpretation
        # summary may sit at either level. Lower-star records are ignored — an
        # internal-cohort BS2 is only trustworthy from the expert panel that owns
        # the data.
        bs2_evidence = 0
        bs2_strength: str | None = None
        if stars >= 3:
            bs2_texts: list[str] = [
                c.text for c in elem.findall(".//Comment") if c.text
            ]
            bs2_texts += [
                d.text for d in elem.findall(".//Attribute[@Type='Description']")
                if d.text
            ]
            bs2_evidence, bs2_strength = _mine_bs2(" ".join(bs2_texts))

        return (var_id, chrom, pos, ref, alt, gene, hgvs_c, hgvs_p,
                aa_change, codon_pos, clinsig, rev_status, stars, last_eval,
                affected, functional, segregation, bs2_evidence, bs2_strength)
    except Exception:
        return None
