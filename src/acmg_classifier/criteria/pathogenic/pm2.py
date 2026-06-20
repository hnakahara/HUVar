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

# Poisson exact upper one-sided 95% confidence limits, λ_U = 0.5·χ²₀.₉₅(2(k+1)),
# indexed by observed allele count k. Used by the Cardiomyopathy/HCM PM2 rule,
# which requires the UPPER bound of the 95% CI of the GrpMax-subpopulation AF
# (= λ_U / AN_grpmax) to be <= the threshold. gnomAD displays only the FAF (the
# CI LOWER bound), so we reconstruct the upper bound from the GrpMax AC/AN. These
# values reproduce the VCEP's published AC/AN equivalence table for the 0.00004
# cutoff (e.g. k=1: 4.744/120000 = 3.95e-5; k=4: 9.154/230000 = 3.98e-5).
_POISSON_U95 = (
    2.996, 4.744, 6.296, 7.754, 9.154, 10.513, 11.842, 13.148, 14.435,
    15.705, 16.962, 18.208, 19.443, 20.669, 21.886, 23.097, 24.301,
    25.499, 26.692, 27.879, 29.062,
)


def _upper_af_95(ac: int | None, an: int | None) -> float | None:
    """Upper one-sided 95% CI bound of the allele frequency for *ac* observed in
    *an* alleles (Poisson exact for small counts; a conservative normal-ish
    extension for large counts, which are common variants that fail PM2 anyway).
    Returns None when AC/AN are unavailable."""
    if ac is None or an is None or an <= 0 or ac < 0:
        return None
    if ac < len(_POISSON_U95):
        lam_u = _POISSON_U95[ac]
    else:
        import math
        # Wilson/score-style upper extension; only reached for high AC (common
        # variants), where the precise value is immaterial — PM2 won't apply.
        lam_u = ac + 1.645 * math.sqrt(ac) + 1.353
    return lam_u / an


class PM2Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Per-gene VCEP PM2 cutoff / strength (Moderate for a handful) / basis
        # (FAF vs raw AF). Genes absent from the spec use the global defaults.
        from acmg_classifier.criteria.pm2_genes import PM2Spec
        self._spec = PM2Spec(cfg.disease_prevalence_tsv)

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

    def _subpop_block(self, rule, gd, threshold: float) -> str | None:
        """Reason PM2 must be withheld because the highest-subpopulation metric
        exceeds the threshold (the FAF95 lower bound under-states it), or None.

        * ``point`` (RUNX1): the GrpMax POINT AF must be < threshold.
        * ``ci95``  (HCM): the UPPER 95% CI of the GrpMax AF (from GrpMax AC/AN)
          must be < threshold; falls back to the POINT AF when AC/AN are absent
          (a DB built before the grpmax columns)."""
        if rule is None or not rule.subpop_mode:
            return None
        if rule.subpop_mode == "point":
            if gd.popmax_af is not None and gd.popmax_af >= threshold:
                return (f"GrpMax point AF={gd.popmax_af:.2e} ≥ {threshold} "
                        "(subpopulation exceeds)")
            return None
        if rule.subpop_mode == "ci95":
            upper = _upper_af_95(gd.ac_grpmax, gd.an_grpmax)
            if upper is None:  # GrpMax AC/AN unavailable → POINT-AF proxy
                if gd.popmax_af is not None and gd.popmax_af >= threshold:
                    return (f"GrpMax point AF={gd.popmax_af:.2e} ≥ {threshold} "
                            "(95% CI unavailable)")
                return None
            if upper >= threshold:
                return (f"GrpMax 95%CI upper={upper:.2e} ≥ {threshold} "
                        "(subpopulation exceeds)")
        return None

    def _zygosity_block(self, rule, gd) -> str | None:
        """Reason PM2 must be withheld because gnomAD shows more homo-/hemizygotes
        than the VCEP tolerates (SLC6A8 0, OTC <=1, SCID genes 0 homozygotes,
        ABCD1 0 hemizygotes), or None when within the allowed count."""
        if rule is None or not rule.zyg_scope:
            return None
        hom = gd.nhomalt or 0
        hemi = gd.nhemi or 0
        if rule.zyg_scope == "hom":
            n, label = hom, f"{hom} homozygote(s)"
        elif rule.zyg_scope == "hemi":
            n, label = hemi, f"{hemi} hemizygote(s)"
        else:  # homhemi (combined)
            n, label = hom + hemi, f"{hom + hemi} homo-/hemizygote(s)"
        if n > rule.zyg_max:
            return f"{label} in gnomAD (> {rule.zyg_max} allowed)"
        return None

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        gd = annotation.gnomad
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None
        rule = self._spec.get(gene)
        # Strength is Supporting by default (SVI); a few VCEPs set PM2 Moderate.
        strength = rule.strength if rule else CriterionStrength.SUPPORTING

        # Total absence from gnomAD is the strongest possible PM2 signal — it
        # implies neither presence nor a failed filter, just no record at all.
        if gd is None:
            return CriteriaResult.met(
                ACMGCriterion.PM2, strength, "Absent from gnomAD (no record)",
            )
        # A failed gnomAD QC filter (most often AC0 — zero high-quality
        # observations) does NOT preclude PM2: such a variant is effectively
        # absent / extremely rare, exactly what PM2 captures. The eRepo benchmark
        # showed 376 PM2 false negatives from blanket-blocking on filter failure,
        # 98% with AF < 1e-5. We instead judge a filter-failed record on rarity
        # like any other; a genuinely common filter-failed call (e.g. an inflated
        # segmental-duplication AF) still fails the threshold below, so this adds
        # no false positives. The filter status is surfaced in the evidence.
        filter_note = "" if gd.filter_pass else " (gnomAD QC filter not passed)"

        # Homozygote/hemizygote ceiling (SLC6A8, OTC, the SCID genes, ABCD1, …):
        # PM2 is withheld when gnomAD shows more homo-/hemizygotes than the VCEP
        # tolerates, regardless of the allele frequency. Computed once and applied
        # to every "met" path below.
        zyg_block = self._zygosity_block(rule, gd)

        # Comparison metric: most VCEPs (and the global default) judge PM2 on the
        # RAW grpmax allele frequency, but some state the cutoff on the GrpMax
        # Filtering Allele Frequency (FAF95) — e.g. HNF1A/HNF4A/GCK. Prefer the
        # spec's metric; fall back through popmax→global AF when it is absent.
        if rule and rule.use_faf:
            value = gd.faf95_popmax
            if value is None:
                value = gd.popmax_af if gd.popmax_af is not None else (gd.af or 0.0)
            metric = "FAF95"
        else:
            value = gd.popmax_af
            if value is None:
                value = gd.af if gd.af is not None else 0.0
            metric = "AF"

        # Non-cancer subset (ENIGMA BRCA1/2): the VCEP judges absence on gnomAD's
        # non-cancer subset, so a variant present only in cancer cohorts still
        # qualifies. Use the non-cancer AF when available; fall back to the overall
        # value (and note it) when the gnomAD DB predates the column — graceful
        # degradation identical to af_xy/ac_xx.
        absent = value == 0.0 or gd.ac == 0
        if rule and rule.subset == "non_cancer":
            if gd.af_non_cancer is not None:
                value = gd.af_non_cancer
                metric = "AF(non-cancer)"
                absent = value == 0.0
            else:
                metric = "AF [non-cancer subset unavailable → overall]"

        # AC=0 means "observed only in samples that failed QC" — effectively
        # absent for ACMG purposes, so we treat it the same as value == 0.
        if absent:
            if zyg_block:  # rare edge: FAF95≈0 yet homozygotes recorded
                return CriteriaResult.not_met(
                    ACMGCriterion.PM2, f"{zyg_block}{filter_note}",
                )
            using_nc = (rule and rule.subset == "non_cancer"
                        and gd.af_non_cancer is not None)
            detail = "non-cancer subset" if using_nc else "AC=0"
            return CriteriaResult.met(
                ACMGCriterion.PM2, strength,
                f"Absent from gnomAD ({detail}){filter_note}",
            )

        # Threshold: the VCEP's per-gene cutoff when set (threshold 0 means the
        # variant must be ABSENT — only the AC=0 branch above qualifies), else
        # the inheritance-aware global default.
        if rule and rule.threshold is not None:
            threshold = rule.threshold
            basis = f"{gene} VCEP"
        else:
            threshold, basis = self._threshold(annotation)

        if value < threshold:
            # Highest-subpopulation correction for the deflated low-AC FAF95
            # (FAF95 is the CI LOWER bound), then the homo-/hemizygote ceiling —
            # either withholds PM2 despite the low frequency.
            block = self._subpop_block(rule, gd, threshold) or zyg_block
            if block:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM2,
                    f"gnomAD {metric}={value:.2e} < {threshold} but {block} "
                    f"[{basis}]{filter_note}",
                )
            return CriteriaResult.met(
                ACMGCriterion.PM2, strength,
                f"gnomAD {metric}={value:.2e} < {threshold} [{basis}]{filter_note}",
            )
        return CriteriaResult.not_met(
            ACMGCriterion.PM2,
            f"gnomAD {metric}={value:.2e} ≥ {threshold} [{basis}]{filter_note}",
        )
