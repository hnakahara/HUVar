"""Manual pathogenic criteria: PS2, PM3, PM6, PP4.

PS3 (functional) and PP1 (cosegregation) have dedicated evaluators that mine ClinVar
SCV comments (with manual supplement override), so they are no longer handled here.
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# These four criteria cannot be derived from public databases alone:
#   PS2 = de novo (requires confirmed parental testing)
#   PM3 = in-trans pathogenic allele (requires phasing data)
#   PM6 = assumed de novo (parental testing not confirmed)
#   PP4 = phenotype-genotype match (requires clinician judgement)
# They are emitted as not_met by default and only triggered when a curator
# supplies an entry through the supplement TSV.
_MANUAL_CRITERIA = (
    ACMGCriterion.PS2,
    ACMGCriterion.PM3,
    ACMGCriterion.PM6,
    ACMGCriterion.PP4,
)


class ManualPathogenicEvaluator(CriterionEvaluator):
    """Reads manual evidence from supplement TSV for criteria requiring human curation."""

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> list[CriteriaResult]:  # type: ignore[override]
        # Returns a list (not a single CriteriaResult) so a single evaluator
        # can emit all four manual criteria in one pass. CriteriaRegistry
        # handles the list/scalar branch when collecting results.
        results = []
        sup = supplement or []
        for criterion in _MANUAL_CRITERIA:
            entries = [e for e in sup if e.criterion == criterion]
            if entries:
                # If the curator entered the same criterion more than once,
                # only the first row is honoured — supplement rows are
                # expected to be unique per (variant, criterion).
                entry = entries[0]
                results.append(CriteriaResult.met(criterion, entry.strength, entry.evidence))
            else:
                # Always emit a not_met record so reviewers can see which
                # criteria were explicitly checked but had no curator input.
                results.append(CriteriaResult.not_met(criterion, "No manual evidence provided"))
        return results
