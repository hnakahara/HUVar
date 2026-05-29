"""PVS1 — null variant in gene where LoF is disease mechanism (ClinGen 2019 decision tree)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Gate for the PVS1 decision tree: only these VEP consequences are candidates
# for "predicted null variant" per ACMG 2015. INFRAME_INSERTION/DELETION are
# handled by PM4, not PVS1. SPLICE_REGION (±1-2 outside donor/acceptor core)
# is intentionally excluded — its impact is uncertain and is covered by PP3
# splice predictors instead.
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
        # Cheap pre-filter before invoking the full decision tree: if the
        # primary transcript's consequence is not a candidate LoF type, PVS1
        # cannot fire and we skip the expensive NMD / last-exon / cryptic-
        # splice analysis in pvs1/decision_tree.py.
        pc = annotation.primary_consequence
        if pc is None or pc.consequence not in _LOF_CONSEQUENCES:
            return CriteriaResult.not_met(ACMGCriterion.PVS1, "Not a LoF consequence")

        # The ClinGen SVI 2019 decision tree is large enough (NMD prediction,
        # last-exon position, biologically-relevant transcript check, etc.)
        # that it lives in its own sub-package. Import is local to keep the
        # criteria layer free of the heavy pvs1 dependencies until needed.
        from acmg_classifier.pvs1.decision_tree import evaluate_pvs1
        strength, evidence = evaluate_pvs1(variant, annotation, self._cfg)
        if strength == CriterionStrength.NOT_MET:
            return CriteriaResult.not_met(ACMGCriterion.PVS1, evidence)
        return CriteriaResult.met(ACMGCriterion.PVS1, strength, evidence)
