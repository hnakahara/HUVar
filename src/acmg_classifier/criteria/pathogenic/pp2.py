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
        eligible, evidence = query_pp2_eligible(
            self._cfg.clinvar_sqlite, pc.gene_symbol, mis_z=mis_z,
        )
        if not eligible:
            return CriteriaResult.not_met(ACMGCriterion.PP2, evidence)
        return CriteriaResult.met(ACMGCriterion.PP2, CriterionStrength.SUPPORTING, evidence)
