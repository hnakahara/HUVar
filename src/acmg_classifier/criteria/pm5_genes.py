"""Per-gene PM5 specifications from ClinGen VCEPs (``disease_prevalence.tsv``).

Three PM5 columns are loaded here, the PM5 counterpart to
:class:`~acmg_classifier.criteria.pp2_genes.PP2Applicability`:

* ``pm5_grantham`` — Grantham-distance gate operator: ``ge`` (candidate >=
  comparator) or ``gt`` (strictly greater — PIK3R1, RYR1). Blank = no gate.
* ``pm5_excludes`` — criteria PM5 may not be combined with for the gene
  (``PM1`` or ``PM1,PS1``; e.g. RASopathy/Cardiomyopathy genes, RUNX1, DICER1).
* ``pm5_max`` — strength ceiling: ``Supporting`` when the VCEP only allows
  PM5_Supporting (ATM, CDH1, PALB2); blank = default Moderate ceiling.
* ``pm5_lp`` — ``no`` when the VCEP offers no Supporting PM5 strength, so the
  same-codon comparator must reach **Pathogenic**; a Likely-pathogenic-only
  comparator must not trigger PM5 (PTEN, VHL, KCNQ1, RASopathy genes, …).
  Blank = a Likely-pathogenic comparator is accepted (at Supporting).
"""
from __future__ import annotations

import csv
from pathlib import Path

from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength

GE = "ge"  # candidate Grantham distance >= comparator
GT = "gt"  # candidate Grantham distance strictly > comparator
_VALID_OP = {GE, GT}


class PM5Spec:
    """VCEP PM5 specifications per gene, loaded once from the TSV.

    A missing file or column degrades to "no spec" (every gene resolves to its
    empty default), so minimal setups and older TSVs keep plain PM5 behaviour.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._op: dict[str, str] = {}
        self._excludes: dict[str, tuple[ACMGCriterion, ...]] = {}
        self._max: dict[str, CriterionStrength] = {}
        self._require_p: set[str] = set()
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                op = (row.get("pm5_grantham") or "").strip().lower()
                if op in _VALID_OP:
                    self._op[gene] = op
                excl = self._parse_excludes(row.get("pm5_excludes") or "")
                if excl:
                    self._excludes[gene] = excl
                cap = (row.get("pm5_max") or "").strip().lower()
                cap_map = {
                    "supporting": CriterionStrength.SUPPORTING,
                    "moderate": CriterionStrength.MODERATE,
                    "strong": CriterionStrength.STRONG,
                }
                if cap in cap_map:
                    self._max[gene] = cap_map[cap]
                if (row.get("pm5_lp") or "").strip().lower() == "no":
                    self._require_p.add(gene)

    @staticmethod
    def _parse_excludes(raw: str) -> tuple[ACMGCriterion, ...]:
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
        return tuple(out)

    def operator(self, gene: str | None) -> str:
        """Grantham operator for *gene*: ``ge`` / ``gt`` / "" (no gate)."""
        if not gene:
            return ""
        return self._op.get(gene, "")

    def excludes(self, gene: str | None) -> tuple[ACMGCriterion, ...]:
        """Criteria PM5 may not be combined with for *gene* (may be empty)."""
        if not gene:
            return ()
        return self._excludes.get(gene, ())

    def max_strength(self, gene: str | None) -> CriterionStrength | None:
        """PM5 strength ceiling for *gene*: ``Supporting`` / ``Moderate`` /
        ``Strong`` from the VCEP, or ``None`` when no VCEP covers it (the
        evaluator then defaults to a Moderate ceiling — PM5_Strong is granted
        only when a VCEP explicitly allows it)."""
        if not gene:
            return None
        return self._max.get(gene)

    def requires_pathogenic(self, gene: str | None) -> bool:
        """True if *gene*'s VCEP requires a Pathogenic same-codon comparator
        (no Supporting PM5 strength) — a Likely-pathogenic-only comparator must
        not trigger PM5."""
        return bool(gene) and gene in self._require_p
