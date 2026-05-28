"""PM2 — absent or very rare in population databases (SVI update: Supporting by default)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# SVI recommendation: use FAF95 < 0.0001 for Supporting (dominant default).
_FAF95_ABSENT = 0.0001
# Recessive / X-linked phenotypes tolerate higher carrier frequencies, so a
# pathogenic variant (e.g. a founder allele) can sit well above the dominant
# threshold. Use a relaxed FAF95 cutoff for genes flagged recessive in the map.
_FAF95_RECESSIVE = 0.005


class PM2Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def _threshold(self, annotation: AnnotationData) -> tuple[float, str]:
        """Pick the FAF95 cutoff from the gene's inheritance; fall back to dominant."""
        from acmg_classifier.local_db.inheritance_db import (
            load_inheritance_map,
            is_recessive,
        )
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None
        if gene:
            inh = load_inheritance_map(self._cfg.gene_inheritance_tsv).get(gene)
            if is_recessive(inh):
                return _FAF95_RECESSIVE, f"recessive ({inh})"
        return _FAF95_ABSENT, "dominant/unknown"

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        gd = annotation.gnomad
        if gd is None:
            return CriteriaResult.met(
                ACMGCriterion.PM2,
                CriterionStrength.SUPPORTING,
                "Absent from gnomAD (no record)",
            )
        if not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.PM2, "gnomAD filter failed")

        faf = gd.faf95_popmax
        if faf is None:
            faf = gd.popmax_af or gd.af or 0.0

        if faf == 0.0 or gd.ac == 0:
            return CriteriaResult.met(
                ACMGCriterion.PM2,
                CriterionStrength.SUPPORTING,
                "Absent from gnomAD (AC=0)",
            )

        threshold, basis = self._threshold(annotation)
        if faf < threshold:
            return CriteriaResult.met(
                ACMGCriterion.PM2,
                CriterionStrength.SUPPORTING,
                f"gnomAD FAF95={faf:.2e} < {threshold} [{basis}]",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.PM2,
            f"gnomAD FAF95={faf:.2e} ≥ {threshold} [{basis}]",
        )
