"""Gene -> inheritance map for inheritance-aware PM2 thresholds.

Loads a simple TSV (``data/shared/gene_inheritance.tsv`` at runtime; the
version-controlled source is generated to ``resources/shared/``):

    gene<TAB>inheritance
    MVK<TAB>AR
    HBB<TAB>AR
    RHAG<TAB>AD/AR
    G6PD<TAB>XL

Inheritance codes follow CGD conventions (AD, AR, AD/AR, XL, XLR, XLD, ...).
A gene is treated as "recessive" (eligible for a relaxed PM2 frequency
threshold) when its code contains AR or XL — i.e. a phenotype in which
unaffected carriers/hemizygote-tolerant frequencies are expected.
"""
from __future__ import annotations
import csv
from functools import lru_cache
from pathlib import Path

import structlog

log = structlog.get_logger()

_GENE_COLS = ("gene", "gene_symbol", "symbol")
_INH_COLS = ("inheritance", "inh", "moi")


@lru_cache(maxsize=8)
def load_inheritance_map(tsv_path: Path) -> dict[str, str]:
    """Return {gene_symbol: inheritance_code}. Empty dict if the file is absent."""
    if not tsv_path.exists():
        log.warning("inheritance_map_missing", path=str(tsv_path))
        return {}

    result: dict[str, str] = {}
    try:
        with tsv_path.open(encoding="utf-8") as fh:
            reader = csv.DictReader(
                (line for line in fh if not line.lstrip().startswith("#")),
                delimiter="\t",
            )
            field_lower = {(c or "").strip().lower(): c for c in (reader.fieldnames or [])}
            gene_col = next((field_lower[c] for c in _GENE_COLS if c in field_lower), None)
            inh_col = next((field_lower[c] for c in _INH_COLS if c in field_lower), None)
            if gene_col is None or inh_col is None:
                log.error("inheritance_map_bad_header", fields=reader.fieldnames)
                return {}
            for row in reader:
                gene = (row.get(gene_col) or "").strip()
                inh = (row.get(inh_col) or "").strip()
                if gene and inh:
                    result[gene] = inh
    except Exception as exc:
        log.error("inheritance_map_error", error=str(exc))
        return {}

    log.info("inheritance_map_loaded", genes=len(result), path=str(tsv_path))
    return result


def is_recessive(inheritance: str | None) -> bool:
    """True if the inheritance code implies a *purely* recessive (or X-linked)
    phenotype — i.e. eligible for the relaxed PM2 frequency threshold.

    A dominant component (AD) disqualifies the gene: combined codes such as
    ``AD/AR`` are treated as dominant (strict threshold), because a dominant
    pathogenic allele must stay rare even when the gene also causes a recessive
    disease. Only ``AR`` / ``XL`` (X-linked, incl. XLR/XLD) without any AD
    component qualify for the relaxed cutoff.
    """
    if not inheritance:
        return False
    code = inheritance.upper()
    if "AD" in code:  # dominant component present -> strict threshold
        return False
    return "AR" in code or "XL" in code
