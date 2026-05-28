"""
PP2 -- missense variant in a gene with a low rate of benign missense variation
and in which missense is a common mechanism of disease.

Gene eligibility is derived from ClinVar (clinvar_sqlite.query_pp2_eligible):
a gene qualifies when it has enough P/LP missense records and a low benign
missense fraction. Only missense variants in eligible genes receive PP2.
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class PP2Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # 1. Manual supplement override
        for e in (supplement or []):
            if e.criterion == ACMGCriterion.PP2:
                return CriteriaResult.met(ACMGCriterion.PP2, e.strength, e.evidence)

        pc = annotation.primary_consequence
        if pc is None or pc.consequence != ConsequenceType.MISSENSE:
            return CriteriaResult.not_met(ACMGCriterion.PP2, "Not a missense variant")

        # 2. ClinVar + gnomAD missense constraint (Z-score) gene eligibility
        from acmg_classifier.local_db.clinvar_sqlite import query_pp2_eligible
        mis_z = annotation.gnomad.mis_z if annotation.gnomad else None
        eligible, evidence = query_pp2_eligible(
            self._cfg.clinvar_sqlite, pc.gene_symbol, mis_z=mis_z,
        )
        if not eligible:
            return CriteriaResult.not_met(ACMGCriterion.PP2, evidence)
        return CriteriaResult.met(ACMGCriterion.PP2, CriterionStrength.SUPPORTING, evidence)
