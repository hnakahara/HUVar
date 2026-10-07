"""
PS4 -- prevalence in affected individuals significantly increased over controls.

Curated-evidence-only criterion (Genome Medicine revision). ClinVar submissions
cannot be de-duplicated at the level of individuals (the same patient may be
reported by several laboratories, in literature-only submissions or via
GenomeConnect) and do not indicate whether the reported individual is an
affected proband, a carrier of a recessive condition or someone ascertained by
screening. PS4 is therefore applied only from

1. curator-supplied evidence (user curation or the eRepo supplement), or
2. PS4 applied by a ClinGen expert panel to this variant in ClinVar (>=3 stars).

The number of P/LP ClinVar submissions reporting affected individuals is still
reported in the evidence text, for the curator's information only.
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry



class PS4Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        from acmg_classifier.criteria.curated_sources import (
            SRC_INFO, from_expert_panel, from_supplement,
        )
        for r in (from_supplement(ACMGCriterion.PS4, supplement),
                  from_expert_panel(ACMGCriterion.PS4, self._cfg, variant)):
            if r is not None:
                return r
        from acmg_classifier.local_db.clinvar_sqlite import query_affected_cases
        n = query_affected_cases(
            self._cfg.clinvar_sqlite, variant.chrom, variant.pos, variant.ref, variant.alt,
        )
        note = (f"; {SRC_INFO} {n} P/LP ClinVar submission(s) report affected individuals "
                "(not de-duplicated; review case-level data)") if n else ""
        return CriteriaResult.not_met(
            ACMGCriterion.PS4, f"PS4 requires curated case-level evidence{note}"
        )
