"""tabix query for FuncVEP-CTI precomputed scores (PP3/BP4).

FuncVEP (Kayaalp et al., Nat Genet 2026) publishes precomputed scores for every
possible missense variant. Two tabix-indexed layouts are read, chosen by the
file's extension (see :pyattr:`Config.funcvep_scores`):

* ``.tsv.gz`` — the 6-column file built by ``scripts/setup_data.py --with-funcvep``:

      chrom  pos  ref  alt  FuncVEP_CTI  gene

* ``.vcf.gz`` — a VCF carrying ``CTI=`` and ``GENE=`` in INFO, i.e. the file that
  VEP ``--custom`` consumes. Reading it directly lets an existing copy of the
  scores be reused instead of keeping a second multi-GB conversion of them.

Licence: PolyForm Strict 1.0.0 — non-commercial use only. The orchestrator opens
this file only when ``insilico_tool`` is FUNCVEP, so a commercial-use run never
reads it (mirroring the BayesDel/CADD licence gate).
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional

import structlog

from acmg_classifier.models.annotation import FuncVEPData
from acmg_classifier.utils.tabix import TabixReader, open_tabix, fetch_region

log = structlog.get_logger()


def is_vcf(path) -> bool:
    """Whether the score file is a VCF (INFO CTI=/GENE=) rather than the 6-column TSV."""
    name = str(path).lower()
    return name.endswith(".vcf.gz") or name.endswith(".vcf.bgz") or name.endswith(".vcf")


def _info_value(info: str, key: str) -> Optional[str]:
    """Value of ``key=`` in a VCF INFO field ("CTI=0.93;GENE=BRCA1")."""
    for item in info.split(";"):
        name, sep, value = item.partition("=")
        if sep and name == key:
            return value
    return None


def _vcf_row(fields: list[str]) -> Optional[tuple[str, str, str, float, str]]:
    """(pos, ref, alt, CTI, gene) from a VCF record, or None when unusable.

    A multi-allelic ALT ("A,T") is skipped rather than guessed at: the scores are
    written one record per alternate allele, so a comma here means the file was
    not produced the way this reader expects."""
    if len(fields) < 8:
        return None
    pos, ref, alt, info = fields[1], fields[3], fields[4], fields[7]
    if "," in alt:
        return None
    raw = _info_value(info, "CTI")
    if raw is None:
        return None
    try:
        score = float(raw)
    except ValueError:
        return None
    return pos, ref, alt, score, (_info_value(info, "GENE") or "")


def _match_funcvep(
    lines: Iterable[str], pos: int, ref: str, alt: str, gene: Optional[str] = None,
    vcf: bool = False,
) -> Optional[FuncVEPData]:
    """Pick the CTI score from tabix region lines matching (pos, ref, alt).

    FuncVEP scores a variant per gene, so a locus inside overlapping genes
    carries several rows with different scores. We prefer the row whose gene
    matches the evaluated transcript's symbol; failing that a single row (or
    several agreeing rows) is used, and genuinely conflicting rows yield no
    score rather than an arbitrary one (same rule as the VAS implementation).
    """
    rows: list[tuple[float, str]] = []
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if vcf:
            parsed = _vcf_row(fields)
            if parsed is None:
                continue
            f_pos, f_ref, f_alt, score, row_gene = parsed
        else:
            if len(fields) < 5:
                continue
            f_pos, f_ref, f_alt = fields[1], fields[2], fields[3]
            try:
                score = float(fields[4])
            except ValueError:
                # A single malformed row must not break annotation.
                continue
            row_gene = fields[5].strip() if len(fields) > 5 else ""
        if f_ref != ref or f_alt != alt:
            continue
        try:
            if int(f_pos) != pos:
                continue
        except ValueError:
            continue
        if score != score:  # NaN
            continue
        rows.append((score, row_gene))

    if not rows:
        return None
    if gene:
        for score, row_gene in rows:
            if row_gene and row_gene == gene:
                return FuncVEPData(score=score, gene=row_gene, match="gene")
    if len(rows) == 1:
        return FuncVEPData(score=rows[0][0], gene=rows[0][1] or None, match="single")
    if len({score for score, _ in rows}) == 1:
        return FuncVEPData(score=rows[0][0], gene=rows[0][1] or None, match="uniform")
    # Several genes disagree and none is the evaluated one: report the ambiguity
    # instead of a score, so PP3/BP4 abstain rather than pick a side.
    return FuncVEPData(score=None, gene=None, match="ambiguous")


def query_funcvep(
    tsv_gz_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene: Optional[str] = None,
) -> Optional[FuncVEPData]:
    """Query the FuncVEP-CTI TSV (tabix-indexed) for one SNV.

    Returns None when the file is absent, the variant is uncovered, or the
    lookup raises — PP3/BP4 then report "no in-silico score available" instead
    of aborting annotation.
    """
    if not tsv_gz_path.exists():
        log.warning("funcvep_missing", path=str(tsv_gz_path))
        return None

    try:
        with open_tabix(tsv_gz_path) as tf:
            return _match_funcvep(fetch_region(tf, chrom, pos, pos), pos, ref, alt, gene,
                                  vcf=is_vcf(tsv_gz_path))
    except Exception as exc:
        log.error("funcvep_error", error=str(exc))
    return None


def query_funcvep_reader(
    reader: TabixReader,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene: Optional[str] = None,
) -> Optional[FuncVEPData]:
    """Batch variant of query_funcvep using a persistent thread-local handle."""
    if not reader.exists():
        log.warning("funcvep_missing", path=str(reader.path))
        return None
    try:
        return _match_funcvep(reader.fetch(chrom, pos, pos), pos, ref, alt, gene,
                              vcf=is_vcf(reader.path))
    except Exception as exc:
        log.error("funcvep_error", error=str(exc))
    return None
