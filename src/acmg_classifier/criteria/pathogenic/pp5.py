"""PP5 -- classified as pathogenic by reputable source (ClinVar >=2 stars).

DEPRECATED — NOT registered in CriteriaRegistry and never applied automatically.
ClinGen's Sequence Variant Interpretation Working Group recommends discontinuing
PP5 (and BP6): assertions decoupled from the underlying primary data risk
double-counting with evidence-based criteria (here PS1/PS3/PM5 already consume
ClinVar), so they may distort the final classification.

Reference: Biesecker LG & Harrison SM, ClinGen SVI Working Group. "The ACMG/AMP
reputable source criteria for the interpretation of sequence variants."
Genet Med 2018;20:1687-1688. doi:10.1038/gim.2018.42

The evaluator is retained only for backward compatibility (enum / strength map /
manual-supplement override, which remains the lab director's responsibility).
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_P_LP = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}


class PP5Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        for rec in annotation.clinvar_vcf:
            if rec.star_rating >= 2 and rec.clinical_significance in _P_LP:
                return CriteriaResult.met(
                    ACMGCriterion.PP5,
                    evidence=f"ClinVar {rec.star_rating} stars: {rec.clinical_significance}",
                )
        return CriteriaResult.not_met(ACMGCriterion.PP5, "No ClinVar >=2 stars P/LP record")
