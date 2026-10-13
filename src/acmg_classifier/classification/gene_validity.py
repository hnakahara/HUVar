"""Gene-disease validity constraint on the final classification.

Following ClinGen guidance (Strande et al. 2017; PMID 31732716 discussion of
gene-disease validity in variant classification), variants in genes whose
ClinGen gene-disease validity curations are *all* Limited, Disputed, Refuted or
No Known Disease Relationship are not classified above VUS. Genes without any
ClinGen curation are flagged but not capped.

Genes with an implemented VCEP specification are not capped: a VCEP is formed
for genes whose relationship with the panel's disease is established, while the
ClinGen gene-disease validity curations of such a gene may concern other diseases
only (e.g. PIK3CA: hereditary breast carcinoma, Refuted; POLG: Leigh syndrome,
Limited). HUVar does not know the patient's disease, so the curations are shown
as a note instead (Reviewer 1, Genome Medicine revision).
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
from typing import Optional

from acmg_classifier.models.enums import Pathogenicity

_PATHOGENIC_SIDE = {Pathogenicity.PATHOGENIC, Pathogenicity.LIKELY_PATHOGENIC}


@lru_cache(maxsize=4)
def _vcep_genes(path: Optional[Path]) -> dict[str, str]:
    """gene -> VCEP name, from the per-gene VCEP table (disease_prevalence.tsv)."""
    out: dict[str, str] = {}
    try:
        if path is None or not Path(path).exists():
            return out
        with open(path, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                g = (row.get("gene_symbol") or "").strip()
                if g:
                    out[g] = (row.get("source_vcep") or "").strip() or "VCEP"
    except Exception:  # noqa: BLE001 — a bad table must not break classification
        return {}
    return out


def apply_gene_validity(
    gene: Optional[str],
    cls_2015: Pathogenicity,
    cls_bayes: Pathogenicity,
    cfg,
) -> tuple[Pathogenicity, Pathogenicity, list[str]]:
    """Return ``(cls_2015, cls_bayes, warnings)`` after the gene-disease validity cap."""
    path = getattr(cfg, "clingen_gene_validity_csv", None)
    if path is None or not gene or not path.exists():
        return cls_2015, cls_bayes, []
    from acmg_classifier.local_db.clingen_gene_db import gene_validity_status
    status, classes = gene_validity_status(path, gene)
    vcep = _vcep_genes(getattr(cfg, "disease_prevalence_tsv", None)).get(gene)
    if status == "unsupported" and vcep:
        return cls_2015, cls_bayes, [
            f"GENE_VALIDITY_NOTE: {gene} ClinGen gene-disease validity curations "
            f"({', '.join(sorted(set(classes)))}) are all Limited or weaker, but the gene has a "
            f"VCEP specification ({vcep}); not capped — confirm that the disease under "
            "consideration is the one covered by the VCEP"]
    if status == "unsupported":
        capped = []
        if cls_2015 in _PATHOGENIC_SIDE:
            capped.append(f"2015 {cls_2015.value}")
            cls_2015 = Pathogenicity.VUS
        if cls_bayes in _PATHOGENIC_SIDE:
            capped.append(f"Bayesian {cls_bayes.value}")
            cls_bayes = Pathogenicity.VUS
        msg = (f"GENE_VALIDITY_LIMITED: {gene} ClinGen gene-disease validity "
               f"({', '.join(sorted(set(classes)))}) does not support P/LP classification")
        if capped:
            msg += f"; capped at VUS (was {', '.join(capped)})"
        return cls_2015, cls_bayes, [msg]
    if status == "not_curated" and (cls_2015 in _PATHOGENIC_SIDE or cls_bayes in _PATHOGENIC_SIDE):
        return cls_2015, cls_bayes, [
            f"GENE_VALIDITY_NOT_CURATED: {gene} has no ClinGen gene-disease validity "
            "curation at Moderate or above; confirm the gene-disease relationship"
        ]
    return cls_2015, cls_bayes, []
