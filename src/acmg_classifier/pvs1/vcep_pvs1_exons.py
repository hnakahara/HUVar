"""Optional exon-aware refinement of the VCEP PVS1 canonical-splice strength.

The flat per-gene splice strength in :mod:`acmg_classifier.pvs1.vcep_pvs1`
(``VERY_STRONG`` for most genes) over-calls canonical ±1,2 splice variants whose
skipped exon is in-frame and/or non-critical — several VCEP decision trees score
those at Strong or Moderate (e.g. DICER1 exon 10 → Strong; exons 5/15/18/22 →
Moderate; CDKL5 exon 17 → Moderate / exon 18 → Strong).

This module lets a reviewer encode those per-(gene, skipped-exon) strengths in a
small TSV that *overrides* the flat default. It is deliberately OPT-IN and
data-driven: with no TSV (or no matching row) the caller keeps the flat default,
so behaviour is unchanged until a reviewer supplies a verified table.

The exon a canonical splice variant skips is derived from VEP's ``intron = n/N``
field per the VCEP convention (also stated verbatim in the ACADVL/GAMT trees):

    * a DONOR (+1/+2) variant in intron *n* skips the exon 5' of it  → exon *n*
    * an ACCEPTOR (-1/-2) variant in intron *n* skips the exon 3' of it → exon *n+1*

Reviewers should populate the override TSV using the coordinate table emitted by
``scripts/build_vcep_pvs1_exons.py`` (which lists, per coding exon, whether
skipping is in-frame and what % of the protein it removes) so the exon numbering
matches the RefSeq transcript VEP annotates against — avoiding the off-by-one
trap from genes with a non-coding exon 1 (HNF4A, CDKL5, MECP2, DICER1, …).

Override TSV columns: ``gene``, ``skipped_exon``, ``strength`` (one of
very_strong / strong / moderate / supporting / na), optional ``note``.
"""
from __future__ import annotations

import csv
from pathlib import Path

from acmg_classifier.models.enums import ConsequenceType, CriterionStrength

_STR_BY_NAME = {
    "very_strong": CriterionStrength.VERY_STRONG,
    "pvs1": CriterionStrength.VERY_STRONG,
    "strong": CriterionStrength.STRONG,
    "pvs1_strong": CriterionStrength.STRONG,
    "moderate": CriterionStrength.MODERATE,
    "pvs1_moderate": CriterionStrength.MODERATE,
    "supporting": CriterionStrength.SUPPORTING,
    "pvs1_supporting": CriterionStrength.SUPPORTING,
    "na": CriterionStrength.NOT_MET,
    "n/a": CriterionStrength.NOT_MET,
    "not_met": CriterionStrength.NOT_MET,
}


def skipped_exon(consequence: ConsequenceType, intron: str | None) -> int | None:
    """Exon predicted to be skipped by a canonical ±1,2 splice variant, from the
    VEP ``intron = n/N`` field. Donor → exon n; acceptor → exon n+1. Returns
    None when the intron index cannot be parsed."""
    if not intron:
        return None
    try:
        n = int(intron.split("/")[0])
    except (ValueError, IndexError):
        return None
    if consequence == ConsequenceType.SPLICE_DONOR:
        return n
    if consequence == ConsequenceType.SPLICE_ACCEPTOR:
        return n + 1
    return None


class SpliceExonOverrides:
    """Per-(gene, skipped-exon) PVS1 splice-strength overrides loaded from a TSV.
    A missing file leaves the table empty, so every lookup misses and the caller
    falls back to the flat per-gene default (behaviour unchanged)."""

    def __init__(self, tsv_path: Path | None) -> None:
        self._table: dict[tuple[str, int], CriterionStrength] = {}
        if tsv_path is not None and tsv_path.exists():
            self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        with tsv_path.open(encoding="utf-8") as fh:
            # Drop leading '#' comment lines so the real column header (not a
            # comment) is what csv.DictReader uses for field names.
            lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
            for row in csv.DictReader(lines, delimiter="\t"):
                gene = (row.get("gene") or "").strip()
                exon_raw = (row.get("skipped_exon") or "").strip()
                name = (row.get("strength") or "").strip().lower()
                if not gene or gene.startswith("#") or not exon_raw:
                    continue
                strength = _STR_BY_NAME.get(name)
                if strength is None:
                    continue
                try:
                    exon = int(exon_raw)
                except ValueError:
                    continue
                self._table[(gene, exon)] = strength

    def __bool__(self) -> bool:
        return bool(self._table)

    def lookup(self, gene: str, exon: int | None) -> CriterionStrength | None:
        if exon is None:
            return None
        return self._table.get((gene, exon))
