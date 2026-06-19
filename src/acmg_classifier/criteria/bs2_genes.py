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

from acmg_classifier.models.enums import CriterionStrength

APPLICABLE = "applicable"
NOT_APPLICABLE = "not_applicable"

# Strength label (cspec) -> ACMG strength, used by the ``bs2_strength`` count→
# strength tiers (e.g. "Strong:6,Supporting:3").
_STRENGTH_LABELS = {
    "verystrong": CriterionStrength.VERY_STRONG,
    "strong": CriterionStrength.STRONG,
    "moderate": CriterionStrength.MODERATE,
    "supporting": CriterionStrength.SUPPORTING,
}
_STRENGTH_RANK = {
    CriterionStrength.VERY_STRONG: 4,
    CriterionStrength.STRONG: 3,
    CriterionStrength.MODERATE: 2,
    CriterionStrength.SUPPORTING: 1,
}


class BS2Applicability:
    """VCEP BS2 applicability and inheritance mode per gene, loaded from the TSV.

    A missing file or column degrades to "no VCEP data" (every gene resolves to
    "" / no modes), so the evaluator falls back to its mode-agnostic heuristic.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._status: dict[str, str] = {}
        self._modes: dict[str, frozenset[str]] = {}
        self._count: dict[str, int] = {}
        self._tiers: dict[str, tuple[tuple[CriterionStrength, int], ...]] = {}
        self._female_only: set[str] = set()
        self._hom_only: set[str] = set()
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
                raw_count = (row.get("bs2_count") or "").strip()
                if raw_count:
                    try:
                        self._count[gene] = int(raw_count)
                    except ValueError:
                        pass
                tiers = self._parse_tiers(row.get("bs2_strength") or "")
                if tiers:
                    self._tiers[gene] = tiers
                # A VCEP whose BS2 counts only females (e.g. TP53: ">=8 unrelated
                # females ... without cancer"). The evaluator then counts female
                # carriers (gnomAD AC_XX) instead of all-sex carriers.
                if (row.get("bs2_female_only") or "").strip() in ("1", "true", "yes"):
                    self._female_only.add(gene)
                # A dominant gene with incomplete penetrance (BMPR2, PIK3R2)
                # whose VCEP scores BS2 on homozygotes only — healthy
                # heterozygous carriers do not count toward benign evidence.
                if (row.get("bs2_hom_only") or "").strip() in ("1", "true", "yes"):
                    self._hom_only.add(gene)

    @staticmethod
    def _parse_tiers(raw: str) -> tuple[tuple[CriterionStrength, int], ...]:
        """Parse a ``bs2_strength`` cell ("Strong:6,Supporting:3") into
        (strength, min_count) pairs sorted strongest-first. Malformed pairs are
        skipped defensively."""
        out: list[tuple[CriterionStrength, int]] = []
        for part in raw.split(","):
            part = part.strip()
            if not part or ":" not in part:
                continue
            label, _, num = part.partition(":")
            strength = _STRENGTH_LABELS.get(label.strip().lower())
            if strength is None:
                continue
            try:
                count = int(num.strip())
            except ValueError:
                continue
            out.append((strength, count))
        out.sort(key=lambda t: _STRENGTH_RANK[t[0]], reverse=True)
        return tuple(out)

    def tiers(self, gene: str | None) -> tuple[tuple[CriterionStrength, int], ...]:
        """Count→strength BS2 tiers for *gene*, strongest-first (e.g.
        ``((STRONG, 6), (SUPPORTING, 3))``); empty when the gene has no tiering
        (the evaluator then uses the single ``bs2_count`` threshold at Strong)."""
        if not gene:
            return ()
        return self._tiers.get(gene, ())

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

    def count(self, gene: str | None) -> int | None:
        """The VCEP's minimum BS2 observation count for *gene*, or None (use the
        global default). When set it overrides the mode-specific count threshold
        — important for cancer panels (CDH1 >=10, TP53 >=8) where the lower
        global default would FALSELY fire BS2 on a pathogenic variant."""
        if not gene:
            return None
        return self._count.get(gene)

    def female_only(self, gene: str | None) -> bool:
        """True if the gene's VCEP counts only females toward BS2 (e.g. TP53).
        The evaluator then counts female carriers (gnomAD AC_XX - nhomalt_XX)
        rather than all-sex carriers."""
        if not gene:
            return False
        return gene in self._female_only

    def hom_only(self, gene: str | None) -> bool:
        """True if the gene's VCEP scores BS2 on homozygotes only (BMPR2,
        PIK3R2). For these incomplete-penetrance dominant genes the evaluator
        counts gnomAD homozygotes (nhomalt) instead of heterozygous carriers."""
        if not gene:
            return False
        return gene in self._hom_only
