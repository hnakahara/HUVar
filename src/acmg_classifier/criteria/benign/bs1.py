"""BS1 -- allele frequency greater than expected for disorder."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Default BS1 threshold (0.5%) is used when no disease-specific value is
# supplied. This is lower than BA1 (5%) because BS1 is "more common than
# expected for the disease", which depends on the disease prevalence —
# 0.5% is a conservative catch-all between BA1 and PM2's strict cutoff.
_DEFAULT_BS1_THRESHOLD = 0.005


class BS1Evaluator(CriterionEvaluator):
    """
    BS1 threshold is disorder-specific (from disease_prevalence.tsv).
    Falls back to a conservative default of 0.5% if gene not in table.
    """

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Thresholds are loaded once at instantiation rather than per-variant
        # to avoid re-parsing the TSV across thousands of variant evaluations.
        self._thresholds: dict[str, float] = self._load_thresholds()

    def _load_thresholds(self) -> dict[str, float]:
        """Read per-gene BS1 cutoffs from disease_prevalence.tsv.

        Missing file is treated as "no per-gene overrides" rather than as a
        hard error, so the evaluator works in minimal setups. Lines with
        malformed numbers are silently skipped — the goal is best-effort
        loading, not strict validation."""
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

        # Prefer FAF95 (conservative), fall back to popmax, then global AF.
        faf = gd.faf95_popmax or gd.popmax_af or gd.af or 0.0

        # Per-gene threshold lookup: if the gene appears in the prevalence
        # table, use its calibrated cutoff; otherwise the default applies.
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
