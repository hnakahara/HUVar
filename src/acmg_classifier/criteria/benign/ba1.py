"""BA1 -- allele frequency >5% in gnomAD (stand-alone benign)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.criteria.allele_frequency import DiseaseThresholds


class BA1Evaluator(CriterionEvaluator):
    """BA1 cutoff is disease-specific (Whiffin/Ware: min(0.05, 10 x maxAF)),
    falling back to the ACMG 2015 stand-alone benign 5% when no per-gene
    parameters are available."""

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._thresholds = DiseaseThresholds(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # Absent or filter-failed gnomAD records cannot establish BA1 — we
        # need a *high* observed frequency, which requires a valid record.
        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BA1, "No valid gnomAD record")

        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else ""
        gt = self._thresholds.get(gene)
        threshold = gt.ba1

        # Homozygote/hemizygote-count OR-clause (SLC6A8, OTC: BA1 if >=10 hom or
        # hemizygotes regardless of frequency). Evaluated first since it is an
        # independent path to BA1.
        if gt.ba1_hom_count is not None:
            n_hh = (gd.nhomalt or 0) + (gd.nhemi or 0)
            if n_hh >= gt.ba1_hom_count:
                return CriteriaResult.met(
                    ACMGCriterion.BA1,
                    evidence=(f"{n_hh} gnomAD homo/hemizygotes >= {gt.ba1_hom_count} "
                              f"for {gene or 'default'} (BA1 count rule)"),
                )

        # Same prefer-FAF95 logic as PM2/PS4: FAF gives the most conservative
        # estimate; only fall back to the POINT popmax AF if FAF is unavailable.
        # The fallback is flagged in the metric label because the point AF lacks
        # the 95% CI sparse-data correction and can over-fire.
        faf = gd.faf95_popmax
        metric = "gnomAD FAF95_popmax"
        if faf is None:
            faf = gd.popmax_af if gd.popmax_af is not None else (gd.af or 0.0)
            metric = "gnomAD popmax AF (FAF95 unavailable)"

        # Non-cancer subset (ENIGMA BRCA1/2): judge BA1 on the non-cancer subset's
        # popmax FAF95 (recomputed at build time from the per-group non-cancer
        # AC/AN). Fall back to the overall FAF95 (noting it) when the companion
        # non-cancer value is unavailable.
        if gt.af_subset == "non_cancer":
            if gd.faf95_non_cancer is not None:
                faf = gd.faf95_non_cancer
                metric = "gnomAD FAF95 (non-cancer)"
            else:
                metric += " [non-cancer subset unavailable → overall]"

        # X-linked genes whose VCEP defines the cutoff "in males" (RPGR, RS1,
        # ABCD1, SLC6A8, OTC): compare against the male (XY) allele frequency.
        # Fall back to the overall FAF when AF_XY is unavailable (e.g. a gnomAD
        # DB built before the af_xy column).
        if gt.af_basis == "males" and gd.af_xy is not None:
            faf = gd.af_xy
            metric = "gnomAD AF_XY (males)"

        # Point-estimate basis: a few VCEPs define BA1/BS1 on the grpmax/popmax
        # POINT allele frequency (not FAF95). Honour that for af_basis="popmax"
        # genes unless globally disabled (Config.popmax_af_basis=False → FAF95).
        if (gt.af_basis == "popmax" and self._cfg.popmax_af_basis
                and gd.popmax_af is not None):
            faf = gd.popmax_af
            metric = "gnomAD popmax AF (point)"

        # 3 significant figures rather than 4 fixed decimals — disease-specific
        # BA1 cutoffs can be ~1e-3/1e-5, which ".4f" rounds misleadingly. The
        # comparison uses full float precision regardless of this display.
        if faf >= threshold:
            return CriteriaResult.met(
                ACMGCriterion.BA1,
                evidence=f"{metric}={faf:.3g} >= BA1 threshold {threshold:.3g} for {gene or 'default'}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BA1,
            f"{metric}={faf:.3g} < BA1 threshold {threshold:.3g}",
        )
