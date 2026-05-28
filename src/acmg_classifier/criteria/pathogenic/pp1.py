"""
PP1 -- cosegregation with disease in multiple affected family members.

Evidence sources (in priority order):
1. Manual supplement entry (curated) -- takes precedence.
2. ClinVar SCV free-text comments describing cosegregation (text-mined).

Strength is kept at Supporting: rigorous PP1 up-weighting (Moderate/Strong) requires
counting informative meioses, which cannot be derived from free-text comments.
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class PP1Evaluator(CriterionEvaluator):
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
            if e.criterion == ACMGCriterion.PP1:
                return CriteriaResult.met(ACMGCriterion.PP1, e.strength, e.evidence)

        # 2. ClinVar text-mined cosegregation evidence
        from acmg_classifier.local_db.clinvar_sqlite import query_segregation_evidence
        n = query_segregation_evidence(
            self._cfg.clinvar_sqlite,
            variant.chrom, variant.pos, variant.ref, variant.alt,
        )
        if n < 1:
            return CriteriaResult.not_met(
                ACMGCriterion.PP1, "No ClinVar SCV describing cosegregation"
            )
        return CriteriaResult.met(
            ACMGCriterion.PP1,
            CriterionStrength.SUPPORTING,
            evidence=f"{n} ClinVar SCV(s) report cosegregation with disease",
        )
