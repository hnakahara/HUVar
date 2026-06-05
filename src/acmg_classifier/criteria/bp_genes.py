"""Per-gene BP1 / BP3 applicability from ClinGen VCEP specs (``disease_prevalence.tsv``).

* ``bp1`` / ``bp1_target`` — BP1 ("variant type unlikely to be pathogenic in this
  gene") applicability and the TARGET consequence. Most VCEPs decline BP1; those
  that apply it target ``missense`` (PALB2, APC, BRCA1/2) or — for gain-of-function
  RASopathy genes where loss-of-function is benign — ``truncating``.
* ``bp3`` — BP3 (in-frame indel in a repetitive region) applicability. A VCEP that
  declined BP3 resolves to ``not_applicable`` and suppresses the heuristic.

The BP counterpart to :class:`~acmg_classifier.criteria.pp2_genes.PP2Applicability`.
"""
from __future__ import annotations

import csv
from pathlib import Path

APPLICABLE = "applicable"
NOT_APPLICABLE = "not_applicable"


class BPApplicability:
    """VCEP BP1/BP3 applicability per gene, loaded once from the TSV. A missing
    file or column degrades to "no VCEP data" (every gene resolves to "")."""

    def __init__(self, tsv_path: Path) -> None:
        self._bp1: dict[str, str] = {}
        self._bp1_target: dict[str, str] = {}
        self._bp1_exclude: dict[str, list[tuple[int, int]]] = {}
        self._bp1_strong: set[str] = set()
        self._bp1_no_splice: set[str] = set()
        self._bp3: dict[str, str] = {}
        self._bp3_regions: dict[str, list[tuple[int, int]]] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                bp1 = (row.get("bp1") or "").strip().lower()
                if bp1 in (APPLICABLE, NOT_APPLICABLE):
                    self._bp1[gene] = bp1
                target = (row.get("bp1_target") or "").strip().lower()
                if target in ("missense", "truncating", "broad"):
                    self._bp1_target[gene] = target
                excl = _parse_ranges(row.get("bp1_exclude") or "")
                if excl:
                    self._bp1_exclude[gene] = excl
                if (row.get("bp1_strength") or "").strip().lower() == "strong":
                    self._bp1_strong.add(gene)
                if (row.get("bp1_no_splice") or "").strip().lower() == "yes":
                    self._bp1_no_splice.add(gene)
                bp3 = (row.get("bp3") or "").strip().lower()
                if bp3 in (APPLICABLE, NOT_APPLICABLE):
                    self._bp3[gene] = bp3
                regions = _parse_ranges(row.get("bp3_regions") or "")
                if regions:
                    self._bp3_regions[gene] = regions

    def bp1(self, gene: str | None) -> str:
        """VCEP BP1 status: ``applicable`` / ``not_applicable`` / ""."""
        if not gene:
            return ""
        return self._bp1.get(gene, "")

    def bp1_target(self, gene: str | None) -> str:
        """BP1 target consequence for *gene*: ``missense`` / ``truncating`` / ""."""
        if not gene:
            return ""
        return self._bp1_target.get(gene, "")

    def bp1_excluded(self, gene: str | None, position: int | None) -> bool:
        """True if *position* falls in a region where BP1 must NOT apply (APC's
        β-catenin repeat; the BRCA1/2 clinically-important functional domains)."""
        if not gene or position is None:
            return False
        return any(a <= position <= b for a, b in self._bp1_exclude.get(gene, ()))

    def bp1_is_strong(self, gene: str | None) -> bool:
        """True if BP1 is applied at Strong for *gene* (BRCA1/2)."""
        return bool(gene) and gene in self._bp1_strong

    def bp1_requires_no_splice(self, gene: str | None) -> bool:
        """True if BP1 requires no predicted splice impact (BRCA1/2: SpliceAI<=0.1)."""
        return bool(gene) and gene in self._bp1_no_splice

    def bp3(self, gene: str | None) -> str:
        """VCEP BP3 status: ``applicable`` / ``not_applicable`` / ""."""
        if not gene:
            return ""
        return self._bp3.get(gene, "")

    def bp3_in_region(self, gene: str | None, position: int | None) -> bool | None:
        """For a gene whose BP3 is region-restricted (RPGR ORF15, FOXG1 poly-AA),
        whether *position* is inside an allowed region. ``None`` when the gene has
        no region restriction (the generic repeat heuristic applies)."""
        regions = self._bp3_regions.get(gene or "")
        if not regions:
            return None
        return position is not None and any(a <= position <= b for a, b in regions)


def _parse_ranges(raw: str) -> list[tuple[int, int]]:
    """Parse ";"-joined "a-b" residue ranges (e.g. "1021-1035;2-101")."""
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
