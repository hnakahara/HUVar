"""Registry that instantiates and runs all criterion evaluators."""
from __future__ import annotations

from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Allele-frequency criteria are mutually exclusive (ClinGen SVI): a variant
# gets at most ONE. Ordered by frequency-evidence strength — BA1 (>5%,
# stand-alone benign) outranks BS1 (benign) which outranks PM2 (rare,
# pathogenic-supporting).
_AF_EXCLUSIVE_PRIORITY = (ACMGCriterion.BA1, ACMGCriterion.BS1, ACMGCriterion.PM2)


def _apply_af_mutual_exclusion(results: list[CriteriaResult]) -> None:
    """Enforce BA1 > BS1 > PM2 mutual exclusivity, in place.

    The highest-priority triggered allele-frequency criterion wins; any
    lower-priority triggered ones are suppressed (retained in the audit trail,
    contributing 0 points) so the same allele-frequency observation is never
    double-counted or counted in both the benign and pathogenic directions.
    """
    active = next(
        (
            c
            for c in _AF_EXCLUSIVE_PRIORITY
            if any(
                r.criterion == c and r.triggered and not r.suppressed
                for r in results
            )
        ),
        None,
    )
    if active is None:
        return
    losers = set(_AF_EXCLUSIVE_PRIORITY[_AF_EXCLUSIVE_PRIORITY.index(active) + 1:])
    for r in results:
        if r.criterion in losers and r.triggered and not r.suppressed:
            r.suppressed = True
            r.evidence = (r.evidence + f" [suppressed: {active.value} active]").strip()


class CriteriaRegistry:
    """Loads and runs all automated criterion evaluators for a given config."""

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._evaluators: list[CriterionEvaluator] = self._build_evaluators()
        # Gene-specific PP2 co-requirements (e.g. BMPR2 needs PM2 + PP3).
        from acmg_classifier.criteria.pp2_genes import PP2Applicability
        self._pp2 = PP2Applicability(cfg.disease_prevalence_tsv)
        # Gene-specific PM5 exclusions (e.g. RUNX1: not with PM1; DICER1: not
        # with PM1/PS1; RASopathy/Cardiomyopathy genes: not with PM1).
        from acmg_classifier.criteria.pm5_genes import PM5Spec
        self._pm5 = PM5Spec(cfg.disease_prevalence_tsv)

    def _build_evaluators(self) -> list[CriterionEvaluator]:
        # Local imports avoid an import cycle: each evaluator module imports
        # config/models, and config in turn references the criteria layer for
        # default thresholds during validation. Deferring imports here also
        # makes the registry cheap to import in tests that stub specific
        # evaluators.
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
        """Run every registered evaluator and apply post-hoc suppression rules.

        Evaluators are independent (each criterion is evaluated in isolation)
        but ACMG/ClinGen forbid certain combinations from being counted
        together. Those interactions are resolved here, after all results are
        collected, rather than inside individual evaluators — otherwise each
        evaluator would need to inspect every other evaluator's output.
        """
        results: list[CriteriaResult] = []
        for evaluator in self._evaluators:
            result = evaluator.evaluate(variant, annotation, supplement)
            # ManualPathogenicEvaluator / ManualBenignEvaluator may emit
            # multiple results (one per supplement row), so accept either a
            # single CriteriaResult or a list.
            if isinstance(result, list):
                results.extend(result)
            else:
                results.append(result)

        # PVS1 ↔ PP3 mutual exclusion (Walker 2023, ClinGen SVI splicing WG):
        # PVS1 already encodes null-variant/splice-disruption evidence at
        # Very Strong, which subsumes in-silico splicing/missense evidence
        # captured by PP3. Letting both fire would double-count the same
        # mechanistic claim, so PP3 is suppressed (kept in the audit trail but
        # contributes 0 points to the Bayesian sum).
        pvs1_triggered = any(
            r.criterion == ACMGCriterion.PVS1 and r.triggered for r in results
        )
        if pvs1_triggered:
            for r in results:
                if r.criterion == ACMGCriterion.PP3 and r.triggered:
                    r.suppressed = True
                    r.evidence = (r.evidence + " [suppressed: PVS1 active]").strip()

        # PP2 gene-specific co-requirements (e.g. BMPR2 / GN125: "PM2_supporting
        # and PP3 must be met"). When a VCEP makes PP2 conditional on other
        # criteria, suppress PP2 unless every required criterion is itself
        # triggered (and not suppressed) for this variant. Run after the PVS1↔PP3
        # pass so a PP3 already suppressed there counts as not-met here.
        self._apply_pp2_co_requirements(results, annotation)

        # Gene-specific PM5 exclusions (RUNX1: not with PM1; DICER1: not with
        # PM1/PS1). Suppress PM5 when an excluded criterion fired for the gene.
        self._apply_pm5_exclusions(results, annotation)

        # BA1 / BS1 / PM2 are mutually exclusive frequency criteria — keep only
        # the highest-priority one (see _apply_af_mutual_exclusion).
        _apply_af_mutual_exclusion(results)

        return results

    def _apply_pm5_exclusions(
        self, results: list[CriteriaResult], annotation: AnnotationData
    ) -> None:
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None
        excluded = self._pm5.excludes(gene)
        if not excluded:
            return
        active = {r.criterion for r in results if r.triggered and not r.suppressed}
        clash = [c for c in excluded if c in active]
        if not clash:
            return
        with_ = "+".join(c.value for c in clash)
        for r in results:
            if r.criterion == ACMGCriterion.PM5 and r.triggered and not r.suppressed:
                r.suppressed = True
                r.evidence = (r.evidence + f" [suppressed: PM5 not with {with_}]").strip()

    def _apply_pp2_co_requirements(
        self, results: list[CriteriaResult], annotation: AnnotationData
    ) -> None:
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None
        required = self._pp2.requires(gene)
        if not required:
            return
        active = {r.criterion for r in results if r.triggered and not r.suppressed}
        missing = [c for c in required if c not in active]
        if not missing:
            return
        need = "+".join(c.value for c in required)
        for r in results:
            if r.criterion == ACMGCriterion.PP2 and r.triggered and not r.suppressed:
                r.suppressed = True
                r.evidence = (r.evidence + f" [suppressed: PP2 requires {need}]").strip()
