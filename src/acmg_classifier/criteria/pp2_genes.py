"""Per-gene PP2 applicability from ClinGen VCEP specs (``disease_prevalence.tsv``).

The ``pp2`` column records each VCEP's gene-level decision on PP2:

* ``applicable``     — a VCEP designates PP2 usable for the gene.
* ``not_applicable`` — a VCEP carries a PP2 code but declined it for the gene
  (e.g. KCNQ1: benign missense throughout, z-score 1.83).
* (blank/absent)     — no VCEP covers the gene; the PP2 evaluator falls back to
  its statistical heuristic.

This authoritative list is what curbs PP2 over-assignment: the statistical
heuristic alone qualified ~4x as many genes as the eRepo truth set.
"""
from __future__ import annotations

import csv
from pathlib import Path

from acmg_classifier.models.enums import ACMGCriterion

APPLICABLE = "applicable"
NOT_APPLICABLE = "not_applicable"


class PP2Applicability:
    """VCEP PP2 applicability per gene, loaded once from the prevalence TSV.

    A missing file or missing column is treated as "no VCEP data" (every gene
    resolves to ""), so minimal setups and older TSVs degrade gracefully.

    The ``pp2_requires`` column carries gene-specific co-requirements — other
    ACMG criteria that must also be met for PP2 to apply (BMPR2: PM2 + PP3).
    These are enforced post-hoc in the registry, not here.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._by_gene: dict[str, str] = {}
        self._requires: dict[str, list[ACMGCriterion]] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                value = (row.get("pp2") or "").strip().lower()
                if value:
                    self._by_gene[gene] = value
                req = self._parse_requires(row.get("pp2_requires") or "")
                if req:
                    self._requires[gene] = req

    @staticmethod
    def _parse_requires(raw: str) -> list[ACMGCriterion]:
        out: list[ACMGCriterion] = []
        for code in raw.split(","):
            code = code.strip().upper()
            if not code:
                continue
            try:
                crit = ACMGCriterion(code)
            except ValueError:
                continue  # unknown code in the table — skip defensively
            if crit not in out:
                out.append(crit)
        return out

    def get(self, gene: str | None) -> str:
        """VCEP PP2 status for *gene*: ``applicable`` / ``not_applicable`` / ""."""
        if not gene:
            return ""
        return self._by_gene.get(gene, "")

    def requires(self, gene: str | None) -> list[ACMGCriterion]:
        """Criteria that must ALSO be met for PP2 to apply to *gene* (may be
        empty). Enforced as a post-hoc suppression pass in the registry."""
        if not gene:
            return []
        return self._requires.get(gene, [])
