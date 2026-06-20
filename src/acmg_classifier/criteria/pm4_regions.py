"""Per-gene PM4 region / strength rules from ClinGen VCEPs (``pm4_regions.tsv``).

Loads the table produced by ``scripts/build_pm4_regions.py``: positive
strength regions/residues (an in-frame indel impacting one gets that strength),
deny regions (withheld), a per-gene ``region_default`` (the strength for an
in-frame indel matching no region), and a stop-loss strength. The PM4 evaluator
consults this before its flat default.
"""
from __future__ import annotations

import csv
from pathlib import Path

from acmg_classifier.models.enums import CriterionStrength

_STRENGTH = {
    "supporting": CriterionStrength.SUPPORTING,
    "moderate": CriterionStrength.MODERATE,
    "strong": CriterionStrength.STRONG,
    "very_strong": CriterionStrength.VERY_STRONG,
}
# Strongest first — the evaluator awards the highest strength whose region/residue
# contains the indel position.
_RANK = {
    CriterionStrength.VERY_STRONG: 4,
    CriterionStrength.STRONG: 3,
    CriterionStrength.MODERATE: 2,
    CriterionStrength.SUPPORTING: 1,
}


class _GeneRule:
    __slots__ = ("tiers", "deny", "region_default", "stoploss",
                 "conserved_phylop", "deletion_content", "excludes")

    def __init__(self) -> None:
        # list of (strength, ranges:list[(a,b)], residues:frozenset[int])
        self.tiers: list[tuple[CriterionStrength, list, frozenset]] = []
        self.deny: list[tuple[int, int]] = []
        self.region_default: CriterionStrength | None | str = None  # str sentinel "not_met"
        self.stoploss: CriterionStrength | None | str = None        # "not_applicable" sentinel
        self.conserved_phylop: float | None = None
        self.deletion_content: bool = False
        self.excludes: tuple[str, ...] = ()


class PM4Regions:
    """VCEP PM4 region/strength rules per gene, loaded once from the TSV. A
    missing file degrades to "no rule" (every gene → ``has_gene`` False), so the
    evaluator keeps its flat default PM4 behaviour."""

    def __init__(self, tsv_path: Path) -> None:
        self._by_gene: dict[str, _GeneRule] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        try:
            if not tsv_path.exists():
                return
            with tsv_path.open(encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        except (OSError, TypeError, AttributeError):
            # Missing file, or a non-Path (mocked) cfg attribute → no rules.
            return
        for row in rows:
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                kind = (row.get("strength") or "").strip().lower()
                rule = self._by_gene.setdefault(gene, _GeneRule())
                if kind == "region_default":
                    rule.region_default = _default_value(row.get("regions") or "")
                elif kind == "stoploss":
                    rule.stoploss = _stoploss_value(row.get("regions") or "")
                elif kind == "conserved_phylop":
                    try:
                        rule.conserved_phylop = float((row.get("regions") or "").strip())
                    except ValueError:
                        pass
                elif kind == "deletion_content":
                    rule.deletion_content = (row.get("regions") or "").strip().lower() in (
                        "yes", "1", "true",
                    )
                elif kind == "excludes":
                    rule.excludes = tuple(
                        c.strip().upper() for c in (row.get("regions") or "").split(",")
                        if c.strip()
                    )
                elif kind == "deny":
                    rule.deny.extend(_parse_ranges(row.get("regions") or ""))
                elif kind in _STRENGTH:
                    ranges = _parse_ranges(row.get("regions") or "")
                    residues = _parse_residues(row.get("residues") or "")
                    if ranges or residues:
                        rule.tiers.append((_STRENGTH[kind], ranges, frozenset(residues)))

    def has_gene(self, gene: str | None) -> bool:
        return bool(gene) and gene in self._by_gene

    def indel_strength(
        self, gene: str | None, position: int | None
    ) -> CriterionStrength | None | str:
        """PM4 strength for an in-frame indel at *position* in *gene*:

        * a :class:`CriterionStrength` (the strongest matching tier, or the
          region default), or
        * the string ``"not_met"`` (denied region, or region default N/A), or
        * ``None`` when the gene has no PM4 region rule (caller uses the flat
          default)."""
        rule = self._by_gene.get(gene or "")
        if rule is None:
            return None
        if position is not None:
            if any(a <= position <= b for a, b in rule.deny):
                return "not_met"
            best: CriterionStrength | None = None
            for strength, ranges, residues in rule.tiers:
                if position in residues or any(a <= position <= b for a, b in ranges):
                    if best is None or _RANK[strength] > _RANK[best]:
                        best = strength
            if best is not None:
                return best
        # No tier matched → the gene's region default.
        if rule.region_default == "not_met":
            return "not_met"
        return rule.region_default  # CriterionStrength or None (→ Moderate default)

    def stoploss_strength(
        self, gene: str | None
    ) -> CriterionStrength | None | str:
        """Stop-loss PM4 strength for *gene*: a strength, ``"not_applicable"``,
        or ``None`` (no override → flat default)."""
        rule = self._by_gene.get(gene or "")
        return rule.stoploss if rule else None

    def conserved_phylop(self, gene: str | None) -> float | None:
        """PhyloP cutoff above which an in-frame indel is "conserved" and PM4 may
        fire (RPE65/CTLA4/PIK3R1), or ``None`` (no conservation gate)."""
        rule = self._by_gene.get(gene or "")
        return rule.conserved_phylop if rule else None

    def requires_deletion_content(self, gene: str | None) -> bool:
        """True if an in-frame DELETION must contain a known ClinVar P/LP or VUS
        within the deleted range to earn PM4 (the SCID panel rule)."""
        rule = self._by_gene.get(gene or "")
        return bool(rule and rule.deletion_content)

    def excludes(self, gene: str | None) -> tuple[str, ...]:
        """ACMG codes PM4 is mutually exclusive with for *gene* (e.g. PVS1,PP3)."""
        rule = self._by_gene.get(gene or "")
        return rule.excludes if rule else ()


def _default_value(raw: str) -> CriterionStrength | str | None:
    raw = raw.strip().lower()
    if raw == "not_met":
        return "not_met"
    return _STRENGTH.get(raw)


def _stoploss_value(raw: str) -> CriterionStrength | str | None:
    raw = raw.strip().lower()
    if raw == "not_applicable":
        return "not_applicable"
    return _STRENGTH.get(raw)


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


def _parse_residues(raw: str) -> set[int]:
    out: set[int] = set()
    for tok in raw.split(","):
        tok = tok.strip()
        if tok:
            try:
                out.add(int(tok))
            except ValueError:
                continue
    return out
