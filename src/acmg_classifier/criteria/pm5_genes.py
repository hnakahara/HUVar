"""Per-gene PM5 Grantham-distance gating from ClinGen VCEP specs.

A subset of VCEPs require PM5 to clear a Grantham-distance test against the
known same-codon pathogenic/likely-pathogenic comparator: the variant under
evaluation must be chemically *as different or more different* from the
wild-type residue than the comparator. The ``pm5_grantham`` column of
``disease_prevalence.tsv`` records the comparison operator per gene:

* ``ge`` — candidate distance must be >= comparator distance (most VCEPs;
  "equal or greater / equal or worse" wording, e.g. PIK3CD, VHL, HNF1A).
* ``gt`` — candidate distance must be strictly greater (PIK3R1: "higher
  Grantham score"; RYR1: comparator "must be less than" the candidate).
* (blank/absent) — no Grantham gate; the PM5 evaluator uses its plain
  same-codon-different-AA rule.

This is the PM5 counterpart to :class:`~acmg_classifier.criteria.pp2_genes.PP2Applicability`.
"""
from __future__ import annotations

import csv
from pathlib import Path

GE = "ge"  # candidate Grantham distance >= comparator
GT = "gt"  # candidate Grantham distance strictly > comparator
_VALID = {GE, GT}


class PM5Grantham:
    """VCEP PM5 Grantham-gating operator per gene, loaded once from the TSV.

    A missing file or column degrades to "no gating" (every gene resolves to
    ""), so minimal setups and pre-Grantham TSVs keep the plain PM5 behaviour.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._by_gene: dict[str, str] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                op = (row.get("pm5_grantham") or "").strip().lower()
                if op in _VALID:
                    self._by_gene[gene] = op

    def operator(self, gene: str | None) -> str:
        """Grantham operator for *gene*: ``ge`` / ``gt`` / "" (no gate)."""
        if not gene:
            return ""
        return self._by_gene.get(gene, "")
