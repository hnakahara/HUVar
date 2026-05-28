"""BS1 -- allele frequency greater than expected for disorder."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_DEFAULT_BS1_THRESHOLD = 0.005


class BS1Evaluator(CriterionEvaluator):
    """
    BS1 threshold is disorder-specific (from disease_prevalence.tsv).
    Falls back to a conservative default of 0.5% if gene not in table.
    """

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._thresholds: dict[str, float] = self._load_thresholds()

    def _load_thresholds(self) -> dict[str, float]:
        tsv = self._cfg.disease_prevalence_tsv
        if not tsv.exists():
            return {}
        thresholds: dict[str, float] = {}
        import csv
        with tsv.open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = row.get("gene_symbol", "").strip()
                threshold = row.get("bs1_threshold", "").strip()
                if gene and threshold:
                    try:
                        thresholds[gene] = float(threshold)
                    except ValueError:
                        pass
        return thresholds

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BS1, "No valid gnomAD record")

        faf = gd.faf95_popmax or gd.popmax_af or gd.af or 0.0

        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else ""
        threshold = self._thresholds.get(gene, _DEFAULT_BS1_THRESHOLD)

        if faf >= threshold:
            return CriteriaResult.met(
                ACMGCriterion.BS1,
                evidence=f"gnomAD FAF95={faf:.4f} >= BS1 threshold {threshold:.4f} for {gene or 'default'}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BS1,
            f"gnomAD FAF95={faf:.4f} < BS1 threshold {threshold:.4f}",
        )
