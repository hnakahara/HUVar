"""Per-gene REVEL PP3/BP4 thresholds from ClinGen VCEP specs (``disease_prevalence.tsv``).

Many VCEPs state a gene-specific REVEL cutoff for PP3/BP4 instead of the
genome-wide ClinGen SVI / Pejaver 2022 defaults. Those cutoffs are mined into
``disease_prevalence.tsv`` (see ``scripts/build_disease_thresholds.py``) as up
to three tiers per direction:

* ``revel_pp3_supporting`` / ``revel_pp3_moderate`` / ``revel_pp3_strong`` —
  pathogenic-side cutoffs; PP3 fires at the tier whose ``REVEL >= cutoff`` holds.
* ``revel_bp4_supporting`` / ``revel_bp4_moderate`` / ``revel_bp4_strong`` —
  benign-side cutoffs; BP4 fires at the tier whose ``REVEL <= cutoff`` holds.

A VCEP that only grants PP3/BP4 at Supporting fills just the ``*_supporting``
column; the evaluator then caps the gene at Supporting (it will NOT promote to
Moderate/Strong on a high/low score, matching the VCEP). A gene with no REVEL
columns at all keeps the global Pejaver tiers.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from acmg_classifier.models.enums import CriterionStrength


@dataclass(frozen=True)
class RevelRule:
    # strength -> cutoff. PP3 cutoffs are minimums (score >= cutoff); BP4
    # cutoffs are maximums (score <= cutoff). Empty dict on a side means the
    # VCEP did not specialise that direction (use the global default there).
    pp3: dict[CriterionStrength, float]
    bp4: dict[CriterionStrength, float]


_PP3_COLS = {
    "revel_pp3_supporting": CriterionStrength.SUPPORTING,
    "revel_pp3_moderate": CriterionStrength.MODERATE,
    "revel_pp3_strong": CriterionStrength.STRONG,
}
_BP4_COLS = {
    "revel_bp4_supporting": CriterionStrength.SUPPORTING,
    "revel_bp4_moderate": CriterionStrength.MODERATE,
    "revel_bp4_strong": CriterionStrength.STRONG,
}


class RevelSpec:
    """VCEP gene-specific REVEL cutoffs, loaded once from the TSV.

    A missing file/column degrades to "no per-gene rule" for every gene, so the
    PP3/BP4 evaluators simply keep their global Pejaver thresholds.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._by_gene: dict[str, RevelRule] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                pp3 = self._parse_tiers(row, _PP3_COLS)
                bp4 = self._parse_tiers(row, _BP4_COLS)
                if pp3 or bp4:
                    self._by_gene[gene] = RevelRule(pp3, bp4)

    @staticmethod
    def _parse_tiers(
        row: dict[str, str], cols: dict[str, CriterionStrength]
    ) -> dict[CriterionStrength, float]:
        out: dict[CriterionStrength, float] = {}
        for col, strength in cols.items():
            raw = (row.get(col) or "").strip()
            if not raw:
                continue
            try:
                out[strength] = float(raw)
            except ValueError:
                continue
        return out

    def get(self, gene: Optional[str]) -> Optional[RevelRule]:
        if not gene:
            return None
        return self._by_gene.get(gene)
