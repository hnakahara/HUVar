"""
PP1 -- cosegregation with disease in multiple affected family members.

Curated-evidence-only criterion (Genome Medicine revision). The number of
informative meioses cannot be derived from ClinVar submission text, and several
submissions may describe the same family; automated PP1 from text mining
over-assigned the criterion relative to the eRepo. PP1 is applied only from

1. curator-supplied evidence (user curation or the eRepo supplement), or
2. PP1 applied by a ClinGen expert panel to this variant in ClinVar (>=3 stars).

ClinVar submissions mentioning co-segregation are reported for information only.
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
        from acmg_classifier.criteria.curated_sources import (
            SRC_INFO, from_expert_panel, from_supplement,
        )
        for r in (from_supplement(ACMGCriterion.PP1, supplement),
                  from_expert_panel(ACMGCriterion.PP1, self._cfg, variant)):
            if r is not None:
                return r
        from acmg_classifier.local_db.clinvar_sqlite import query_segregation_evidence
        n = query_segregation_evidence(
            self._cfg.clinvar_sqlite, variant.chrom, variant.pos, variant.ref, variant.alt,
        )
        note = (f"; {SRC_INFO} {n} ClinVar submission(s) mention co-segregation") if n else ""
        return CriteriaResult.not_met(
            ACMGCriterion.PP1, f"PP1 requires curated segregation evidence{note}"
        )
