"""Per-gene BS2 applicability and inheritance mode from ClinGen VCEP specs.

BS2 ("observed in a healthy adult") is only meaningful when the VCEP allows
general-population (gnomAD) data for the gene, and the count that matters depends
on the inheritance mode:

* recessive (AR)  → homozygotes
* X-linked (XL)   → hemizygotes
* dominant (AD)   → heterozygous carriers

Both signals are read from ``disease_prevalence.tsv``:

* ``bs2``         — ``applicable`` / ``not_applicable`` (a VCEP that bars
  population data, e.g. RASopathy GN004, resolves to ``not_applicable``); blank
  when no VCEP covers the gene.
* ``inheritance`` — ``AD`` / ``AR`` / ``XL`` (comma-joined when several), used to
  pick which gnomAD count the BS2 evaluator tests.

The BS2 counterpart to :class:`~acmg_classifier.criteria.pp2_genes.PP2Applicability`.
"""
from __future__ import annotations

import csv
from pathlib import Path

APPLICABLE = "applicable"
NOT_APPLICABLE = "not_applicable"


class BS2Applicability:
    """VCEP BS2 applicability and inheritance mode per gene, loaded from the TSV.

    A missing file or column degrades to "no VCEP data" (every gene resolves to
    "" / no modes), so the evaluator falls back to its mode-agnostic heuristic.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._status: dict[str, str] = {}
        self._modes: dict[str, frozenset[str]] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                status = (row.get("bs2") or "").strip().lower()
                if status in (APPLICABLE, NOT_APPLICABLE):
                    self._status[gene] = status
                modes = frozenset(
                    m.strip().upper()
                    for m in (row.get("inheritance") or "").split(",")
                    if m.strip().upper() in ("AD", "AR", "XL")
                )
                if modes:
                    self._modes[gene] = modes

    def status(self, gene: str | None) -> str:
        """VCEP BS2 status for *gene*: ``applicable`` / ``not_applicable`` / ""."""
        if not gene:
            return ""
        return self._status.get(gene, "")

    def modes(self, gene: str | None) -> frozenset[str]:
        """Inheritance modes for *gene* (subset of {AD, AR, XL}); empty if unknown."""
        if not gene:
            return frozenset()
        return self._modes.get(gene, frozenset())
