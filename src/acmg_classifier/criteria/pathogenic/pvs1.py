"""PVS1 — null variant in gene where LoF is disease mechanism (ClinGen 2019 decision tree)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_LOF_CONSEQUENCES = {
    ConsequenceType.FRAMESHIFT,
    ConsequenceType.STOP_GAINED,
    ConsequenceType.SPLICE_ACCEPTOR,
    ConsequenceType.SPLICE_DONOR,
    ConsequenceType.START_LOST,
    ConsequenceType.TRANSCRIPT_ABLATION,
}


class PVS1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        pc = annotation.primary_consequence
        if pc is None or pc.consequence not in _LOF_CONSEQUENCES:
            return CriteriaResult.not_met(ACMGCriterion.PVS1, "Not a LoF consequence")

        # Full ClinGen 2019 decision tree lives in pvs1/ sub-package
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        strength, evidence = evaluate_pvs1(variant, annotation, self._cfg)
        if strength == CriterionStrength.NOT_MET:
            return CriteriaResult.not_met(ACMGCriterion.PVS1, evidence)
        return CriteriaResult.met(ACMGCriterion.PVS1, strength, evidence)
