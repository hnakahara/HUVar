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


# Genes whose VCEP does not permit a text-mined PS3 from ClinVar:
#   * no PS3 code in the VCEP spec at all — PALB2, PDHA1, POLG;
#   * PS3 explicitly "not applicable ... for in vitro assays" (only a rare,
#     case-by-case variant-specific animal model qualifies, which a curator must
#     assert via the manual supplement) — CAPN3, ANO5.
# A manual supplement PS3 still applies (it takes precedence below); only the
# free-text ClinVar fallback is suppressed for these genes.
_PS3_NOT_APPLICABLE = frozenset({"PALB2", "PDHA1", "POLG", "CAPN3", "ANO5"})

# Genes whose VCEP caps PS3 at Supporting — the text-mined fallback must not reach
# Moderate (n>=3) for them. Transcribed from the cspec PS3 strength descriptors.
_PS3_MAX_SUPPORTING = frozenset({
    "ABCD1", "AIPL1", "ETHE1", "F8", "F9", "GALT", "GAMT", "GATM", "HBA2",
    "HBB", "RPE65", "RPGR", "SERPINC1", "SLC19A3", "VHL",
})


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
        # 1. Manual supplement takes precedence over the text-mined fallback:
        #    if a curator has reviewed the underlying paper they can assert a
        #    higher strength (including PS3 Strong, which the text-mined path
        #    cannot reach — see module docstring on OddsPath calibration).
        for e in (supplement or []):
            if e.criterion == ACMGCriterion.PS3:
                return CriteriaResult.met(ACMGCriterion.PS3, e.strength, e.evidence)

        # Gene gate: some VCEPs do not allow a text-mined PS3 at all (no PS3 code,
        # or PS3 not applicable for in vitro assays). Withhold the free-text
        # fallback for those genes (a manual supplement above still applies).
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None
        if gene in _PS3_NOT_APPLICABLE:
            return CriteriaResult.not_met(
                ACMGCriterion.PS3,
                f"{gene}: VCEP does not permit a ClinVar-text PS3 (no PS3 / in vitro N/A)",
            )

        # 2. Fallback: count ClinVar SCVs whose free-text comment matches a
        #    damaging-functional-study pattern. Strength is capped at Moderate
        #    because OddsPath cannot be derived from free text alone.
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
        # Per-gene VCEP cap: several VCEPs allow PS3 only at Supporting.
        if gene in _PS3_MAX_SUPPORTING and strength != CriterionStrength.SUPPORTING:
            strength = CriterionStrength.SUPPORTING
        return CriteriaResult.met(
            ACMGCriterion.PS3,
            strength,
            evidence=f"{n} ClinVar SCV(s) report damaging functional study ({strength.value})",
        )
