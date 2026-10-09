"""Manual benign criteria: BS3, BS4, BP2, BP5."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Benign criteria that cannot be automated from public databases:
#   BS3 = functional studies show no damaging effect (requires literature)
#   BS4 = lack of segregation in affected family members
#   BP2 = observed in trans with pathogenic in dominant gene / in cis
#   BP5 = found in case with alternate molecular basis
# Same supplement-based curator workflow as the pathogenic manual.py.
_MANUAL_CRITERIA = (
    ACMGCriterion.BS3,
    ACMGCriterion.BS4,
    ACMGCriterion.BP2,
    ACMGCriterion.BP5,
)


class ManualBenignEvaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        from acmg_classifier.criteria.tp53_functional import TP53Functional
        # Opt-in (non-commercial data terms): Config.use_tp53_functional.
        self._tp53 = TP53Functional(
            cfg.tp53_functional_tsv if getattr(cfg, "use_tp53_functional", False) is True else None)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> list[CriteriaResult]:  # type: ignore[override]
        # Returns a list (one record per criterion) so the registry sees a
        # not_met entry for criteria the curator did not provide — keeps the
        # audit trail complete.
        results = []
        sup = supplement or []
        for criterion in _MANUAL_CRITERIA:
            entries = [e for e in sup if e.criterion == criterion]
            if entries:
                entry = entries[0]
                results.append(CriteriaResult.met(criterion, entry.strength, entry.evidence))
            elif criterion == ACMGCriterion.BS3 and self._tp53_bs3(annotation) is not None:
                # TP53: BS3 from the VCEP functional flowchart on systematic assay data.
                results.append(self._tp53_bs3(annotation))
            else:
                results.append(CriteriaResult.not_met(criterion, "No manual evidence provided"))
        return results

    def _tp53_bs3(self, annotation: AnnotationData) -> CriteriaResult | None:
        pc = annotation.primary_consequence
        if pc is None or pc.gene_symbol != "TP53" or not self._tp53:
            return None
        from acmg_classifier.criteria.pathogenic.ps3 import tp53_functional_result
        return tp53_functional_result(self._tp53, pc, ACMGCriterion.BS3)
