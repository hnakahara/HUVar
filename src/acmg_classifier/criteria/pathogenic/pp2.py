"""
PP2 -- missense variant in a gene with a low rate of benign missense variation
and in which missense is a common mechanism of disease.

Gene eligibility is decided in priority order:
  1. ClinGen VCEP applicability (``pp2`` column of ``disease_prevalence.tsv``):
     a VCEP's explicit "applicable" / "not applicable" decision is authoritative
     and overrides the statistical heuristic. This is the dominant precision
     lever — most VCEPs declined PP2 for their genes, but the heuristic alone
     ignored that and over-assigned.
  2. For genes no VCEP covers, fall back to ClinVar statistics
     (clinvar_sqlite.query_pp2_eligible): enough P/LP missense and a low benign
     missense fraction (with a gnomAD missense-Z rescue).
Only missense variants in eligible genes receive PP2.
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.pp2_genes import (
    PP2Applicability, APPLICABLE, NOT_APPLICABLE,
)
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


@lru_cache(maxsize=4)
def load_common_missense(path: Path) -> dict[str, int]:
    """``{gene: n gnomAD missense variants meeting the gene's BS1 threshold}``
    from ``pp2_gene_stats.tsv`` (scripts/build_pp2_gene_stats.py). Empty when the
    file is absent, in which case PP2 uses ClinVar B/LB missense only."""
    out: dict[str, int] = {}
    if not isinstance(path, Path) or not path.exists():
        return out
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            try:
                out[r["gene"].strip()] = int(r["common_missense"])
            except (KeyError, ValueError):
                continue
    return out


class PP2Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._vcep = PP2Applicability(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # 1. Manual supplement always wins so curators can assert PP2 for
        #    genes that don't meet the automatic statistical thresholds but
        #    are known by domain experts to be missense-driven.
        for e in (supplement or []):
            if e.criterion == ACMGCriterion.PP2:
                return CriteriaResult.met(ACMGCriterion.PP2, e.strength, e.evidence)

        # PP2 is, by definition, a missense-only criterion.
        pc = annotation.primary_consequence
        if pc is None or pc.consequence != ConsequenceType.MISSENSE:
            return CriteriaResult.not_met(ACMGCriterion.PP2, "Not a missense variant")

        # 2. ClinGen VCEP applicability is authoritative when present: a VCEP
        #    that curated this gene has already decided whether PP2 applies, so
        #    we honour that over the statistical heuristic (which over-fires).
        vcep = self._vcep.get(pc.gene_symbol)
        if vcep == NOT_APPLICABLE:
            return CriteriaResult.not_met(
                ACMGCriterion.PP2,
                f"{pc.gene_symbol}: VCEP designates PP2 not applicable",
            )
        if vcep == APPLICABLE:
            return CriteriaResult.met(
                ACMGCriterion.PP2,
                CriterionStrength.SUPPORTING,
                f"{pc.gene_symbol}: VCEP designates PP2 applicable",
            )

        # 3. No VCEP covers this gene — fall back to ClinVar statistics (low
        #    benign-missense rate relative to P/LP missense) AND/OR the gnomAD
        #    missense Z-score (population-level missense constraint). mis_z may
        #    be None for genes not in the constraint table; the query handles it.
        from acmg_classifier.local_db.clinvar_sqlite import query_pp2_eligible
        mis_z = annotation.gnomad.mis_z if annotation.gnomad else None
        cfg = self._cfg
        stats_path = getattr(cfg, "pp2_gene_stats_tsv", None)
        common = load_common_missense(stats_path).get(pc.gene_symbol, 0) if stats_path else 0
        eligible, evidence = query_pp2_eligible(
            cfg.clinvar_sqlite, pc.gene_symbol, mis_z=mis_z,
            common_missense=common,
            min_path=_num(cfg, "pp2_min_path", 10),
            max_benign_frac=_num(cfg, "pp2_max_benign_frac", 0.05),
            min_mis_z=_num(cfg, "pp2_min_mis_z", 3.09),
            z_max_benign_frac=_num(cfg, "pp2_z_max_benign_frac", 0.15),
        )
        if not eligible:
            return CriteriaResult.not_met(ACMGCriterion.PP2, evidence)
        return CriteriaResult.met(ACMGCriterion.PP2, CriterionStrength.SUPPORTING, evidence)


def _num(cfg, name: str, default):
    """Config value if numeric (tests pass MagicMock configs), else the default."""
    v = getattr(cfg, name, default)
    return v if isinstance(v, (int, float)) and not isinstance(v, bool) else default
