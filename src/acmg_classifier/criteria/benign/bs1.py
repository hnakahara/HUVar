"""BS1 -- allele frequency greater than expected for disorder."""
from __future__ import annotations
import re

from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry
from acmg_classifier.criteria.allele_frequency import DiseaseThresholds

# Bare protein change from an HGVS p. string, tolerant of a transcript prefix
# and ClinVar/VEP parentheses; "*" is normalised to "Ter" so p.Gln368* and
# p.Gln368Ter compare equal. Lower-cased for case-insensitive matching.
_PCHANGE_RE = re.compile(r"p\.\(?([A-Za-z]{3}\d+(?:[A-Za-z]{3}|Ter|\*|=))\)?")


def _norm_pchange(hgvs_p: str | None) -> str:
    if not hgvs_p:
        return ""
    m = _PCHANGE_RE.search(hgvs_p)
    return m.group(1).replace("*", "Ter").lower() if m else ""


class BS1Evaluator(CriterionEvaluator):
    """
    BS1 threshold is disorder-specific, derived from disease_prevalence.tsv
    (Whiffin/Ware maximum credible AF, floored at 0.05%). Falls back to a
    conservative flat default if the gene is absent or lacks parameters.
    """

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Thresholds are loaded once at instantiation rather than per-variant
        # to avoid re-parsing the TSV across thousands of variant evaluations.
        self._thresholds = DiseaseThresholds(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BS1, "No valid gnomAD record")

        # Prefer FAF95 (the conservative 95% CI lower bound). Use it as-is when
        # present — including a legitimate 0.0 — and only fall back to raw
        # popmax/global AF when FAF95 is genuinely absent (None). This mirrors
        # BA1; the previous `faf95 or popmax or af` chain collapsed a real
        # FAF95 of 0.0 to the raw AF, over-firing BS1 on wide-CI rare variants.
        faf = gd.faf95_popmax
        metric = "gnomAD FAF95"
        if faf is None:
            faf = gd.popmax_af if gd.popmax_af is not None else (gd.af or 0.0)
            metric = "gnomAD popmax AF (FAF95 unavailable)"

        # Per-gene threshold lookup: if the gene appears in the prevalence
        # table, use its disease-specific cutoff; otherwise the default applies.
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else ""
        gt = self._thresholds.get(gene)
        threshold = gt.bs1

        # Variant-level BS1 exclusion: a recurrent disease allele the VCEP bars
        # from BS1 regardless of its population frequency (e.g. MYOC p.Gln368Ter,
        # whose ~2.6% allelic contribution is disease-driven, not benign).
        if gt.bs1_exclude:
            excluded = {_norm_pchange(x) for x in gt.bs1_exclude.split(",")}
            if _norm_pchange(pc.hgvs_p if pc else None) in excluded:
                return CriteriaResult.not_met(
                    ACMGCriterion.BS1,
                    f"{gene}: VCEP excludes BS1 for {gt.bs1_exclude} "
                    "(recurrent disease allele)",
                )

        # X-linked "in males" genes (RPGR, RS1): the VCEP cutoff is on the male
        # (XY) allele frequency. Fall back to the overall FAF when AF_XY is
        # unavailable (gnomAD DB predating the af_xy column).
        if gt.af_basis == "males" and gd.af_xy is not None:
            faf = gd.af_xy
            metric = "gnomAD AF_XY (males)"

        # Point-estimate basis (mirrors BA1): VCEPs that define the cutoff on the
        # grpmax/popmax POINT AF rather than FAF95 (af_basis="popmax"), unless
        # globally disabled (Config.popmax_af_basis=False).
        if (gt.af_basis == "popmax" and self._cfg.popmax_af_basis
                and gd.popmax_af is not None):
            faf = gd.popmax_af
            metric = "gnomAD popmax AF (point)"

        # Format with 3 significant figures (not fixed 4 decimals): disease-
        # specific cutoffs are often ~1e-4, where ".4f" rounds both the FAF and
        # the threshold to misleadingly coarse values (0.000358 → "0.0004",
        # 0.000316 → "0.0003"). The comparison itself uses full float precision.
        if faf >= threshold:
            return CriteriaResult.met(
                ACMGCriterion.BS1,
                strength=gt.bs1_strength,
                evidence=f"{metric}={faf:.3g} >= BS1 threshold {threshold:.3g} "
                         f"({gt.bs1_strength.value}) for {gene or 'default'}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.BS1,
            f"{metric}={faf:.3g} < BS1 threshold {threshold:.3g}",
        )
