"""PM2 — absent or very rare in population databases (SVI update: Supporting by default)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# PM2 judges on the RAW grpmax allele frequency, not a FAF/CI estimate:
# ClinGen Hearing Loss VCEP specifies that BA1/BS1 use the Filtering Allele
# Frequency (95% CI lower bound) but PM2 uses the observed gnomAD frequency
# directly. SVI dominant default: raw AF < 0.0001 → Supporting.
_RAW_AF_ABSENT = 0.0001
# Recessive / X-linked phenotypes tolerate higher carrier frequencies, so a
# pathogenic variant (e.g. a founder allele) can sit well above the dominant
# threshold. Use a relaxed cutoff for genes flagged recessive in the map.
_RAW_AF_RECESSIVE = 0.005


class PM2Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def _threshold(self, annotation: AnnotationData) -> tuple[float, str]:
        """Select the FAF95 cutoff based on the gene's inheritance pattern.

        Recessive disease genes tolerate much higher carrier frequencies in
        the general population (founder alleles are routinely seen at 0.1-1%),
        so a strict dominant-disease cutoff would over-trigger PM2 against
        legitimate pathogenic recessive variants. Falling back to the
        dominant threshold when inheritance is unknown is the conservative
        choice — it is harder to *meet* PM2, never easier."""
        from acmg_classifier.local_db.inheritance_db import (
            load_inheritance_map,
            is_recessive,
        )
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None
        if gene:
            inh = load_inheritance_map(self._cfg.gene_inheritance_tsv).get(gene)
            if is_recessive(inh):
                return _RAW_AF_RECESSIVE, f"recessive ({inh})"
        return _RAW_AF_ABSENT, "dominant/unknown"

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        gd = annotation.gnomad
        # Total absence from gnomAD is the strongest possible PM2 signal — it
        # implies neither presence nor a failed filter, just no record at all.
        if gd is None:
            return CriteriaResult.met(
                ACMGCriterion.PM2,
                CriterionStrength.SUPPORTING,
                "Absent from gnomAD (no record)",
            )
        # A failed gnomAD QC filter means we cannot trust the AF estimate, so
        # we abstain rather than assert PM2 on dubious data.
        if not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.PM2, "gnomAD filter failed")

        # PM2 uses the RAW grpmax allele frequency (no FAF/CI), per the ClinGen
        # Hearing Loss VCEP. Prefer the grpmax (popmax) raw AF; fall back to the
        # global raw AF only when grpmax is genuinely absent (None). A real 0.0
        # is kept as-is — it means the variant is unobserved in the grpmax
        # group, which is the strongest PM2 signal.
        raw_af = gd.popmax_af
        if raw_af is None:
            raw_af = gd.af if gd.af is not None else 0.0

        # AC=0 means "observed only in samples that failed QC" — effectively
        # absent for ACMG purposes, so we treat it the same as raw AF == 0.
        if raw_af == 0.0 or gd.ac == 0:
            return CriteriaResult.met(
                ACMGCriterion.PM2,
                CriterionStrength.SUPPORTING,
                "Absent from gnomAD (AC=0)",
            )

        threshold, basis = self._threshold(annotation)
        if raw_af < threshold:
            return CriteriaResult.met(
                ACMGCriterion.PM2,
                CriterionStrength.SUPPORTING,
                f"gnomAD AF={raw_af:.2e} < {threshold} [{basis}]",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.PM2,
            f"gnomAD AF={raw_af:.2e} ≥ {threshold} [{basis}]",
        )
