"""Registry that instantiates and runs all criterion evaluators."""
from __future__ import annotations

from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class CriteriaRegistry:
    """Loads and runs all automated criterion evaluators for a given config."""

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._evaluators: list[CriterionEvaluator] = self._build_evaluators()

    def _build_evaluators(self) -> list[CriterionEvaluator]:
        from acmg_classifier.criteria.pathogenic.pvs1 import PVS1Evaluator
        from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator
        from acmg_classifier.criteria.pathogenic.ps3 import PS3Evaluator
        from acmg_classifier.criteria.pathogenic.ps4 import PS4Evaluator
        from acmg_classifier.criteria.pathogenic.pm1 import PM1Evaluator
        from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
        from acmg_classifier.criteria.pathogenic.pm4 import PM4Evaluator
        from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
        from acmg_classifier.criteria.pathogenic.pp1 import PP1Evaluator
        from acmg_classifier.criteria.pathogenic.pp2 import PP2Evaluator
        from acmg_classifier.criteria.pathogenic.pp3 import PP3Evaluator
        # PP5 (reputable-source) intentionally NOT registered — deprecated by
        # ClinGen SVI (Biesecker & Harrison, Genet Med 2018;20:1687-1688). See pp5.py.
        from acmg_classifier.criteria.pathogenic.manual import ManualPathogenicEvaluator
        from acmg_classifier.criteria.benign.ba1 import BA1Evaluator
        from acmg_classifier.criteria.benign.bs1 import BS1Evaluator
        from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
        from acmg_classifier.criteria.benign.bp3 import BP3Evaluator
        from acmg_classifier.criteria.benign.bp4 import BP4Evaluator
        from acmg_classifier.criteria.benign.bp7 import BP7Evaluator
        from acmg_classifier.criteria.benign.manual import ManualBenignEvaluator

        return [
            PVS1Evaluator(self._cfg),
            PS1Evaluator(self._cfg),
            PS3Evaluator(self._cfg),
            PS4Evaluator(self._cfg),
            PM1Evaluator(self._cfg),
            PM2Evaluator(self._cfg),
            PM4Evaluator(self._cfg),
            PM5Evaluator(self._cfg),
            PP1Evaluator(self._cfg),
            PP2Evaluator(self._cfg),
            PP3Evaluator(self._cfg),
            ManualPathogenicEvaluator(self._cfg),
            BA1Evaluator(self._cfg),
            BS1Evaluator(self._cfg),
            BS2Evaluator(self._cfg),
            BP3Evaluator(self._cfg),
            BP4Evaluator(self._cfg),
            BP7Evaluator(self._cfg),
            ManualBenignEvaluator(self._cfg),
        ]

    def evaluate_all(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> list[CriteriaResult]:
        results: list[CriteriaResult] = []
        for evaluator in self._evaluators:
            result = evaluator.evaluate(variant, annotation, supplement)
            if isinstance(result, list):
                results.extend(result)
            else:
                results.append(result)

        # PVS1 ↔ PP3 mutual exclusion (Walker 2023)
        pvs1_triggered = any(
            r.criterion == ACMGCriterion.PVS1 and r.triggered for r in results
        )
        if pvs1_triggered:
            for r in results:
                if r.criterion == ACMGCriterion.PP3 and r.triggered:
                    r.suppressed = True
                    r.evidence = (r.evidence + " [suppressed: PVS1 active]").strip()

        return results
