"""Per-gene PS1 specifications from ClinGen VCEP specs (``disease_prevalence.tsv``).

* ``ps1``        — ``applicable`` / ``not_applicable``. A VCEP may decline PS1
  for its gene (CDH1: "Not applicable for CDH1"); blank when no VCEP covers it.
* ``ps1_splice`` — the gene's PS1 splice-extension state:

  * ``""``            — no splice extension; PS1 is missense-only (original ACMG;
    e.g. GAA, the HCM genes). A splice/intronic variant must NOT receive PS1.
  * ``canonical``     — splice extension covering canonical sites too (HNF1A,
    GCK, HNF4A, PIK3R1, SLC6A8).
  * ``noncanonical``  — splice extension limited to non-canonical positions
    (InSiGHT MMR, BMPR2, RS1, RUNX1, …); a canonical ±1/±2 change is PVS1
    territory and must not also get PS1.
"""
from __future__ import annotations

import csv
from pathlib import Path


class PS1Spec:
    """VCEP PS1 specifications per gene, loaded once from the TSV. A missing file
    or column degrades to "no restriction" for every gene."""

    def __init__(self, tsv_path: Path) -> None:
        self._not_applicable: set[str] = set()
        self._splice_mode: dict[str, str] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                if (row.get("ps1") or "").strip().lower() == "not_applicable":
                    self._not_applicable.add(gene)
                mode = (row.get("ps1_splice") or "").strip().lower()
                if mode in ("canonical", "noncanonical"):
                    self._splice_mode[gene] = mode

    def is_not_applicable(self, gene: str | None) -> bool:
        """True if the gene's VCEP declined PS1 entirely (e.g. CDH1)."""
        return bool(gene) and gene in self._not_applicable

    def splice_mode(self, gene: str | None) -> str:
        """PS1 splice-extension state for *gene*: "" (missense-only, no splice
        PS1) / ``canonical`` / ``noncanonical``."""
        if not gene:
            return ""
        return self._splice_mode.get(gene, "")
