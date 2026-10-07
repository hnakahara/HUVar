"""ClinGen gene-level curations: Dosage Sensitivity (haploinsufficiency) and
Gene-Disease Validity.

* **Dosage Sensitivity** — ``ClinGen_gene_curation_list_GRCh38.tsv`` from
  https://ftp.clinicalgenome.org/ (tab-separated, ``#``-prefixed preamble; the
  header line contains ``Gene Symbol`` and ``Haploinsufficiency Score``). A
  haploinsufficiency (HI) score of 3 ("sufficient evidence") establishes loss of
  function as a disease mechanism for PVS1 (Abou Tayoun et al. 2018).
* **Gene-Disease Validity** — the CSV export of
  https://search.clinicalgenome.org/kb/gene-validity/download (a short preamble,
  then a header with ``GENE SYMBOL`` and ``CLASSIFICATION``). Genes whose
  curations are *all* Limited / Disputed / Refuted / No Known Disease
  Relationship are not classified above VUS (Strande et al. 2017,
  PMID 31732716 / ClinGen guidance).

Both files are optional: when absent, the corresponding checks are skipped.
"""
from __future__ import annotations

import csv
import io
from functools import lru_cache
from pathlib import Path
from typing import Optional

import structlog

log = structlog.get_logger()

# Gene-disease validity classes that DO support clinical-grade classification.
SUPPORTED_GDV = {"definitive", "strong", "moderate"}
# Classes for which P/LP classification is not supported.
UNSUPPORTED_GDV = {
    "limited", "disputed", "disputing", "refuted", "refuting",
    "no known disease relationship", "no reported evidence",
}


def _read_text(path: Path) -> str:
    raw = path.read_bytes()
    for enc in ("utf-8-sig", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


@lru_cache(maxsize=4)
def load_hi_scores(path: Path) -> dict[str, str]:
    """``{gene_symbol: haploinsufficiency score}`` (scores kept as strings: "3", "2",
    "1", "0", "30", "40"). Empty dict when the file is absent or unparseable."""
    if not path or not path.exists():
        return {}
    lines = _read_text(path).splitlines()
    header_idx = next(
        (i for i, ln in enumerate(lines)
         if "Gene Symbol" in ln and "Haploinsufficiency Score" in ln), None)
    if header_idx is None:
        log.warning("clingen_dosage_bad_header", path=str(path))
        return {}
    header = [h.lstrip("#").strip() for h in lines[header_idx].split("\t")]
    gi, hi = header.index("Gene Symbol"), header.index("Haploinsufficiency Score")
    out: dict[str, str] = {}
    for ln in lines[header_idx + 1:]:
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) > max(gi, hi) and f[gi].strip():
            out[f[gi].strip()] = f[hi].strip()
    return out


@lru_cache(maxsize=4)
def load_gene_validity(path: Path) -> dict[str, list[str]]:
    """``{gene_symbol: [classification, ...]}`` (lower-case), one entry per curation."""
    if not path or not path.exists():
        return {}
    text = _read_text(path)
    lines = text.splitlines()
    header_idx = next(
        (i for i, ln in enumerate(lines)
         if "GENE SYMBOL" in ln.upper() and "CLASSIFICATION" in ln.upper()), None)
    if header_idx is None:
        log.warning("clingen_validity_bad_header", path=str(path))
        return {}
    reader = csv.reader(io.StringIO("\n".join(lines[header_idx:])))
    header = [h.strip().upper() for h in next(reader)]
    gi = header.index("GENE SYMBOL")
    ci = header.index("CLASSIFICATION")
    out: dict[str, list[str]] = {}
    for row in reader:
        if len(row) <= max(gi, ci):
            continue
        gene, cls = row[gi].strip(), row[ci].strip().lower()
        if not gene or set(gene) <= {"+", "-", "="}:
            continue
        out.setdefault(gene, []).append(cls)
    return out


def hi_score(path: Optional[Path], gene: Optional[str]) -> Optional[str]:
    if not gene or not isinstance(path, Path):
        return None
    return load_hi_scores(path).get(gene)


def gene_validity_status(path: Path, gene: Optional[str]) -> tuple[str, list[str]]:
    """``(status, classifications)`` where status is

    * ``"supported"``   — at least one Definitive/Strong/Moderate curation;
    * ``"unsupported"`` — curated, but every curation is Limited or below;
    * ``"not_curated"`` — no ClinGen gene-disease validity curation (or no file).
    """
    if not gene or not isinstance(path, Path):
        return "not_curated", []
    cls = load_gene_validity(path).get(gene, [])
    if not cls:
        return "not_curated", []
    if any(c in SUPPORTED_GDV for c in cls):
        return "supported", cls
    if all(c in UNSUPPORTED_GDV for c in cls):
        return "unsupported", cls
    return "not_curated", cls
