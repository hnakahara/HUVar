"""Per-gene PVS1 applicability from ClinGen VCEP specs (``disease_prevalence.tsv``).

* ``pvs1`` — ``applicable`` / ``not_applicable`` / "". A VCEP declares PVS1
  ``not_applicable`` when loss-of-function is NOT the disease mechanism for the
  gene — i.e. gain-of-function or dominant-negative disorders, where a null
  allele is not pathogenic by the established mechanism. Examples: MYOC
  (misfolded-myocilin POAG), the RASopathy panel (BRAF, KRAS, PTPN11, …), the
  cardiomyopathy genes (MYH7, TNNT2, …), the activating PIK3 genes (PIK3CA, …),
  RYR1, VWF. PVS1 must never fire on a null variant in such a gene.

The PVS1 counterpart to :class:`~acmg_classifier.criteria.ps1_genes.PS1Spec`.
"""
from __future__ import annotations

import csv
from pathlib import Path


class PVS1Applicability:
    """VCEP PVS1 applicability per gene, loaded once from the TSV. A missing file
    or column degrades to "no VCEP data" (every gene resolves to "")."""

    def __init__(self, tsv_path: Path) -> None:
        self._not_applicable: set[str] = set()
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                if (row.get("pvs1") or "").strip().lower() == "not_applicable":
                    self._not_applicable.add(gene)

    def is_not_applicable(self, gene: str | None) -> bool:
        """True if the gene's VCEP declined PVS1 entirely (LoF not the mechanism)."""
        return bool(gene) and gene in self._not_applicable
