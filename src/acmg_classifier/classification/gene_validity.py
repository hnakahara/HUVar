"""Gene-disease validity constraint on the final classification.

Following ClinGen guidance (Strande et al. 2017; PMID 31732716 discussion of
gene-disease validity in variant classification), variants in genes whose
ClinGen gene-disease validity curations are *all* Limited, Disputed, Refuted or
No Known Disease Relationship are not classified above VUS. Genes without any
ClinGen curation are flagged but not capped.
"""
from __future__ import annotations

from typing import Optional

from acmg_classifier.models.enums import Pathogenicity

_PATHOGENIC_SIDE = {Pathogenicity.PATHOGENIC, Pathogenicity.LIKELY_PATHOGENIC}


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
