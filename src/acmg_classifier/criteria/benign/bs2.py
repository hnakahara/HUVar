"""BS2 -- observed in healthy adult individual (gnomAD homozygote/hemizygote counts)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Healthy-individual observation thresholds. The 5-count cutoff is a
# common practical threshold: 1-2 may reflect sample contamination or
# data-processing artefacts, while 5+ is much harder to explain by chance.
# gnomAD already excludes severely-affected individuals, so adult
# homozygotes/hemizygotes there support a benign interpretation for
# fully-penetrant recessive/X-linked conditions.
_MIN_HOMALT = 5   # homozygote observations (recessive interpretation)
_MIN_HEMI = 5     # hemizygote observations (X-linked male interpretation)


class BS2Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BS2, "No valid gnomAD record")

        # `or 0` collapses missing fields (None) into 0 so the comparison
        # always works without a separate None check per branch.
        nhomalt = gd.nhomalt or 0
        nhemi = gd.nhemi or 0

        # Homozygote check first because it covers autosomal recessive
        # (the most common BS2 interpretation). Hemi-only triggering covers
        # X-linked variants where homozygotes are biologically impossible.
        if nhomalt >= _MIN_HOMALT:
            return CriteriaResult.met(
                ACMGCriterion.BS2,
                evidence=f"gnomAD nhomalt={nhomalt} >= {_MIN_HOMALT}",
            )
        if nhemi >= _MIN_HEMI:
            return CriteriaResult.met(
                ACMGCriterion.BS2,
                evidence=f"gnomAD nhemi={nhemi} >= {_MIN_HEMI}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BS2,
            f"gnomAD nhomalt={nhomalt}, nhemi={nhemi} (below threshold)",
        )
