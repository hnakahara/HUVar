"""Per-gene PM2 specifications from ClinGen VCEP specs (``disease_prevalence.tsv``).

* ``pm2_threshold`` — the VCEP's allele-frequency cutoff for the gene. ``"0"``
  means the variant must be ABSENT (only a truly absent / AC=0 variant
  qualifies). Blank → no VCEP cutoff; the evaluator keeps its global default.
* ``pm2_strength``  — ``Moderate`` for the handful of VCEPs that set PM2 at
  Moderate (GAA, LDLR, ETHE1, PDHA1, POLG, SLC19A3, ITGA2B, ITGB3); blank →
  the SVI default of Supporting.
* ``pm2_basis``     — ``faf`` when the VCEP states the cutoff on the GrpMax
  Filtering Allele Frequency (FAF95); blank → the raw popmax allele frequency.
* ``pm2_subpop``    — highest-subpopulation metric mode that corrects the
  deflated low-AC FAF95: ``point`` (RUNX1 — also require the GrpMax POINT AF
  <= threshold) or ``ci95`` (Cardiomyopathy/HCM — require the UPPER bound of the
  95% CI of the GrpMax AF <= threshold, reconstructed from GrpMax AC/AN); blank
  for genes with no such rule.
* ``pm2_zygosity`` — a homozygote/hemizygote ceiling "<scope>:<max>" (scope =
  hom / hemi / homhemi) PM2 also requires (SLC6A8 ``homhemi:0``, OTC
  ``homhemi:1``, the SCID genes / GATM / GAMT ``hom:0``, ABCD1 ``hemi:0``);
  blank when the VCEP states no such requirement.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from acmg_classifier.models.enums import CriterionStrength


@dataclass(frozen=True)
class PM2Rule:
    threshold: Optional[float]      # None → no per-gene cutoff (use global default)
    strength: CriterionStrength     # Moderate or Supporting
    use_faf: bool                   # compare FAF95 instead of raw popmax AF
    subpop_mode: str = ""           # "" / "point" (RUNX1) / "ci95" (HCM)
    zyg_scope: str = ""             # "" / "hom" / "hemi" / "homhemi"
    zyg_max: int = 0                # highest tolerated homo/hemi count
    subset: str = ""                # "" / "non_cancer" (ENIGMA BRCA1/2: judge
                                    # absence on the gnomAD non-cancer subset)
    min_depth: Optional[float] = None  # ENIGMA BRCA1/2: the region must have a
                                    # gnomAD mean read depth >= this for PM2


class PM2Spec:
    """VCEP PM2 cutoff/strength/basis per gene, loaded once from the TSV.

    A missing file/column degrades to "no per-gene rule" for every gene, so the
    PM2 evaluator simply keeps its global thresholds.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._by_gene: dict[str, PM2Rule] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                raw = (row.get("pm2_threshold") or "").strip()
                threshold: Optional[float] = None
                if raw:
                    try:
                        threshold = float(raw)
                    except ValueError:
                        threshold = None
                strength = (
                    CriterionStrength.MODERATE
                    if (row.get("pm2_strength") or "").strip().lower() == "moderate"
                    else CriterionStrength.SUPPORTING
                )
                use_faf = (row.get("pm2_basis") or "").strip().lower() == "faf"
                subpop_mode = (row.get("pm2_subpop") or "").strip().lower()
                if subpop_mode not in ("point", "ci95"):
                    subpop_mode = ""
                zyg_scope, zyg_max = "", 0
                raw_zyg = (row.get("pm2_zygosity") or "").strip().lower()
                if ":" in raw_zyg:
                    scope, _, mx = raw_zyg.partition(":")
                    if scope in ("hom", "hemi", "homhemi"):
                        try:
                            zyg_scope, zyg_max = scope, int(mx)
                        except ValueError:
                            zyg_scope = ""
                subset = (row.get("pm2_subset") or "").strip().lower()
                if subset != "non_cancer":
                    subset = ""
                min_depth: Optional[float] = None
                raw_depth = (row.get("pm2_min_depth") or "").strip()
                if raw_depth:
                    try:
                        min_depth = float(raw_depth)
                    except ValueError:
                        min_depth = None
                # Only record a rule when the gene carries at least one PM2
                # specialisation; otherwise leave it to the global default.
                if (raw or strength == CriterionStrength.MODERATE or use_faf
                        or subpop_mode or zyg_scope or subset or min_depth is not None):
                    self._by_gene[gene] = PM2Rule(
                        threshold, strength, use_faf, subpop_mode, zyg_scope,
                        zyg_max, subset, min_depth,
                    )

    def get(self, gene: Optional[str]) -> Optional[PM2Rule]:
        if not gene:
            return None
        return self._by_gene.get(gene)
