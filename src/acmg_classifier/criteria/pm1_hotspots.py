"""Per-gene PM1 hotspot regions from ClinGen VCEP specs (``pm1_hotspots.tsv``).

PM1 ("mutational hotspot / critical functional domain") is defined per gene by
each VCEP. ``scripts/build_pm1_hotspots.py`` mines those definitions into a table
of residue ranges and explicit residue positions, with the VCEP strength, plus
genes for which the VCEP declared PM1 ``not_applicable``.

The PM1 counterpart to :class:`~acmg_classifier.criteria.pp2_genes.PP2Applicability`.
"""
from __future__ import annotations

import csv
from pathlib import Path

from acmg_classifier.models.enums import CriterionStrength

_STRENGTH = {
    "supporting": CriterionStrength.SUPPORTING,
    "moderate": CriterionStrength.MODERATE,
    "strong": CriterionStrength.STRONG,
}
# Strongest first — the evaluator awards the highest strength whose region
# contains the residue.
_RANK = {
    CriterionStrength.STRONG: 3,
    CriterionStrength.MODERATE: 2,
    CriterionStrength.SUPPORTING: 1,
}


class PM1Hotspots:
    """VCEP PM1 hotspot regions per gene, loaded once from the TSV.

    A missing file degrades to "no curated data" (every gene resolves to no
    hotspot / not not-applicable), so the evaluator falls back to its statistical
    hotspot heuristic for every gene.
    """

    def __init__(self, tsv_path: Path) -> None:
        # gene -> list of (strength, ranges:list[(a,b)], residues:frozenset[int])
        self._by_gene: dict[str, list[tuple[CriterionStrength, list, frozenset]]] = {}
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
                strength_raw = (row.get("strength") or "").strip().lower()
                if strength_raw == "not_applicable":
                    self._not_applicable.add(gene)
                    continue
                strength = _STRENGTH.get(strength_raw)
                if strength is None:
                    continue
                ranges = _parse_ranges(row.get("regions") or "")
                residues = _parse_residues(row.get("residues") or "")
                if not ranges and not residues:
                    continue
                self._by_gene.setdefault(gene, []).append((strength, ranges, residues))

    def is_not_applicable(self, gene: str | None) -> bool:
        """True if the gene's VCEP declared PM1 not applicable (e.g. ABCA4, ATM,
        the RASopathy genes — benign variation throughout / no defined hotspot)."""
        return bool(gene) and gene in self._not_applicable

    def has_gene(self, gene: str | None) -> bool:
        """True if the gene has curated hotspot rows (so the statistical fallback
        should be skipped)."""
        return bool(gene) and gene in self._by_gene

    def lookup(self, gene: str | None, position: int | None) -> CriterionStrength | None:
        """Strongest PM1 strength whose hotspot contains *position*, or None."""
        if not gene or position is None:
            return None
        best: CriterionStrength | None = None
        for strength, ranges, residues in self._by_gene.get(gene, ()):
            if position in residues or any(a <= position <= b for a, b in ranges):
                if best is None or _RANK[strength] > _RANK[best]:
                    best = strength
        return best


def _parse_ranges(raw: str) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for part in raw.split(";"):
        part = part.strip()
        if not part or "-" not in part:
            continue
        a, _, b = part.partition("-")
        try:
            out.append((int(a), int(b)))
        except ValueError:
            continue
    return out


def _parse_residues(raw: str) -> frozenset[int]:
    out: set[int] = set()
    for tok in raw.split(","):
        tok = tok.strip()
        if tok:
            try:
                out.add(int(tok))
            except ValueError:
                continue
    return frozenset(out)
