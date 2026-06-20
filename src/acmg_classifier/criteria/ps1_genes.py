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
* ``ps1_max``    — per-gene PS1 strength ceiling: ``Supporting`` (RMRP, whose
  VCEP downgrades PS1 to Supporting) or ``Moderate``; blank = no cap (PS1 keeps
  its Strong/Moderate comparator-derived strength).
* ``ps1_paralog_group`` — sibling genes whose analogous (same-numbered) residue a
  PS1 same-AA comparator may also come from (RASopathy "highly analogous
  groupings" HRAS/NRAS/KRAS, MAP2K1/MAP2K2, SOS1/SOS2; HBA2↔HBA1). Comma-joined,
  the OTHER genes of the group.
* ``ps1_paralog_strength`` — fixed strength for a paralogue-only PS1 hit:
  ``Moderate`` (HBA2: "PS1_Moderate ... in a paralogue gene"); blank = use the
  comparator-derived strength (RASopathy grants the full Strong/Moderate).
"""
from __future__ import annotations

import csv
from pathlib import Path

from acmg_classifier.models.enums import CriterionStrength

_CAP = {
    "supporting": CriterionStrength.SUPPORTING,
    "moderate": CriterionStrength.MODERATE,
    "strong": CriterionStrength.STRONG,
}


class PS1Spec:
    """VCEP PS1 specifications per gene, loaded once from the TSV. A missing file
    or column degrades to "no restriction" for every gene."""

    def __init__(self, tsv_path: Path) -> None:
        self._not_applicable: set[str] = set()
        self._splice_mode: dict[str, str] = {}
        self._max: dict[str, CriterionStrength] = {}
        self._paralog_group: dict[str, tuple[str, ...]] = {}
        self._paralog_strength: dict[str, CriterionStrength] = {}
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
                cap = (row.get("ps1_max") or "").strip().lower()
                if cap in _CAP:
                    self._max[gene] = _CAP[cap]
                group = tuple(
                    s.strip() for s in (row.get("ps1_paralog_group") or "").split(",")
                    if s.strip()
                )
                if group:
                    self._paralog_group[gene] = group
                pstr = (row.get("ps1_paralog_strength") or "").strip().lower()
                if pstr in _CAP:
                    self._paralog_strength[gene] = _CAP[pstr]

    def is_not_applicable(self, gene: str | None) -> bool:
        """True if the gene's VCEP declined PS1 entirely (e.g. CDH1)."""
        return bool(gene) and gene in self._not_applicable

    def splice_mode(self, gene: str | None) -> str:
        """PS1 splice-extension state for *gene*: "" (missense-only, no splice
        PS1) / ``canonical`` / ``noncanonical``."""
        if not gene:
            return ""
        return self._splice_mode.get(gene, "")

    def max_strength(self, gene: str | None) -> CriterionStrength | None:
        """PS1 strength ceiling for *gene* (``Supporting`` / ``Moderate``), or
        ``None`` when the VCEP does not cap PS1 below its Strong default."""
        if not gene:
            return None
        return self._max.get(gene)

    def paralog_group(self, gene: str | None) -> tuple[str, ...]:
        """Sibling genes whose analogous (same-numbered) residue a PS1 comparator
        may also come from, or empty when the gene has no paralogue PS1 rule."""
        if not gene:
            return ()
        return self._paralog_group.get(gene, ())

    def paralog_strength(self, gene: str | None) -> CriterionStrength | None:
        """Fixed strength for a paralogue-only PS1 hit (HBA2 → Moderate), or
        ``None`` to use the comparator-derived strength (RASopathy)."""
        if not gene:
            return None
        return self._paralog_strength.get(gene)
