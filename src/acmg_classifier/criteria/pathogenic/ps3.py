"""
PS3 -- well-established functional studies show a damaging effect.

Evidence sources (in priority order):
1. Manual supplement entry (curated) -- takes precedence, may reach Strong if the
   curator has performed OddsPath calibration.
2. ClinVar SCV free-text comments describing a damaging functional assay (text-mined).

Strength per Brnich et al. 2019 (Genome Med 13073): PS3 Strong requires OddsPath > 18.7
established from documented pathogenic/benign control counts. That calibration cannot be
derived from ClinVar free text, so text-mined PS3 is capped at Moderate. SCV count is used
only as a weak confidence proxy:
   1-2 SCVs  -> Supporting
   >=3 SCVs  -> Moderate (cap)
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


def _functional_strength(n: int) -> CriterionStrength | None:
    # Capped at Moderate: OddsPath calibration for Strong is unavailable from ClinVar text.
    if n >= 3:
        return CriterionStrength.MODERATE
    if n >= 1:
        return CriterionStrength.SUPPORTING
    return None


class PS3Evaluator(CriterionEvaluator):
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
            if e.criterion == ACMGCriterion.PS3:
                return CriteriaResult.met(ACMGCriterion.PS3, e.strength, e.evidence)

        # 2. ClinVar text-mined functional evidence
        from acmg_classifier.local_db.clinvar_sqlite import query_functional_evidence
        n = query_functional_evidence(
            self._cfg.clinvar_sqlite,
            variant.chrom, variant.pos, variant.ref, variant.alt,
        )
        strength = _functional_strength(n)
        if strength is None:
            return CriteriaResult.not_met(
                ACMGCriterion.PS3, "No ClinVar SCV describing damaging functional study"
            )
        return CriteriaResult.met(
            ACMGCriterion.PS3,
            strength,
            evidence=f"{n} ClinVar SCV(s) report damaging functional study ({strength.value})",
        )
