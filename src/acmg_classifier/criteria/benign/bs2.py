"""BS2 -- observed in a healthy adult individual (gnomAD), inheritance-mode aware.

The count that argues *benign* depends on the gene's inheritance mode:

* recessive (AR)  → homozygotes (`nhomalt`)
* X-linked (XL)   → hemizygotes (`nhemi`)
* dominant (AD)   → heterozygous carriers (`AC - nhomalt`) — a healthy adult
  heterozygote is itself evidence against a (fully-penetrant) dominant disorder.

Whether BS2 may use general-population data at all is a per-gene VCEP decision
(`bs2` column of ``disease_prevalence.tsv``): a VCEP that bars population data
(e.g. RASopathy GN004) resolves to ``not_applicable`` and BS2 is withheld. When
no VCEP covers the gene, the evaluator falls back to the mode-agnostic
homozygote/hemizygote heuristic.
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.bs2_genes import BS2Applicability, NOT_APPLICABLE
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Expert-panel BS2 strength label (as stored by clinvar_builder) -> ACMG strength.
_BS2_STRENGTH = {
    "VeryStrong": CriterionStrength.VERY_STRONG,
    "Strong": CriterionStrength.STRONG,
    "Moderate": CriterionStrength.MODERATE,
    "Supporting": CriterionStrength.SUPPORTING,
}


class BS2Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._vcep = BS2Applicability(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None

        # VCEP gate: a VCEP that declined BS2 or barred population data for the
        # gene withholds the gnomAD-count path (e.g. RASopathy GN004; or CDH1 /
        # TP53 / SERPINC1 whose BS2 needs an internal cohort gnomAD cannot
        # supply). For those genes the only admissible BS2 is a variant-level
        # judgement already made by the expert panel — harvested from ClinVar
        # >=3-star reviews. This path does NOT require a gnomAD record (these
        # internal-cohort variants are often absent from gnomAD).
        if self._vcep.status(gene) == NOT_APPLICABLE:
            return self._clinvar_fallback(variant, gene)

        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BS2, "No valid gnomAD record")

        nhomalt = gd.nhomalt or 0
        nhemi = gd.nhemi or 0
        ac = gd.ac or 0
        het_carriers = max(0, ac - nhomalt)  # individuals carrying >=1 allele (AD)

        # A VCEP whose BS2 counts only females (TP53: ">=8 unrelated females ...
        # without cancer") must not count male carriers. Use gnomAD's female
        # allele/homozygote counts: female carriers = AC_XX - nhomalt_XX. When the
        # gnomAD DB predates these columns (ac_xx is None), we cannot confirm the
        # female count, so BS2 is withheld rather than counting both sexes — the
        # conservative choice that avoids a false benign call on a pathogenic
        # variant (the whole reason these VCEPs restrict to females).
        female_only = self._vcep.female_only(gene)
        if female_only:
            if gd.ac_xx is None:
                return CriteriaResult.not_met(
                    ACMGCriterion.BS2,
                    f"{gene}: VCEP counts only females, but gnomAD female counts "
                    "(AC_XX) unavailable in this DB build",
                )
            het_carriers = max(0, (gd.ac_xx or 0) - (gd.nhomalt_xx or 0))

        # An incomplete-penetrance dominant gene (BMPR2, PIK3R2) scores BS2 on
        # HOMOZYGOTES only — a healthy heterozygote is not evidence against the
        # disease, so counting carriers (the AD default below) would FALSELY fire
        # BS2 on a pathogenic variant. The VCEP states "≥N homozygotes in gnomAD".
        hom_only = self._vcep.hom_only(gene)

        modes = self._vcep.modes(gene)

        # Count→strength tiers (e.g. GUCY2D Strong>=6 / Supporting>=3; BMPR2
        # Strong>=3 / Moderate>=2 / Supporting>=1). When a gene defines tiers the
        # BS2 strength tracks the observed count instead of the flat Strong
        # default — fire at the strongest tier whose count threshold is met.
        tiers = self._vcep.tiers(gene)
        if tiers:
            observed, label = self._observed_count(
                modes, nhomalt, nhemi, het_carriers, hom_only, female_only
            )
            for strength, thr in tiers:  # strongest-first
                if observed >= thr:
                    return self._met(f"{label}={observed} >= {thr}", strength=strength)
            return self._not_met(nhomalt, nhemi, het_carriers)
        # A per-gene VCEP count (e.g. CDH1 >=10, TP53 >=8) overrides ALL mode
        # thresholds; otherwise use the inheritance-mode global defaults.
        vcep_count = self._vcep.count(gene)
        if vcep_count is not None:
            hom_thr = hemi_thr = het_thr = vcep_count
        else:
            hom_thr = self._cfg.bs2_min_homalt
            hemi_thr = self._cfg.bs2_min_hemi
            het_thr = self._cfg.bs2_min_het

        # When the inheritance mode is known, restrict to the count that mode
        # implies (a dominant gene's homozygotes are irrelevant; a recessive
        # gene's heterozygotes are expected carriers). Without mode information
        # fall back to the homozygote/hemizygote heuristic (mode-agnostic).
        if not modes:
            if nhomalt >= hom_thr:
                return self._met(f"nhomalt={nhomalt} >= {hom_thr}")
            if nhemi >= hemi_thr:
                return self._met(f"nhemi={nhemi} >= {hemi_thr}")
            return self._not_met(nhomalt, nhemi, het_carriers)

        if "AR" in modes and nhomalt >= hom_thr:
            return self._met(f"recessive: nhomalt={nhomalt} >= {hom_thr}")
        if "XL" in modes and nhemi >= hemi_thr:
            return self._met(f"X-linked: nhemi={nhemi} >= {hemi_thr}")
        if "AD" in modes:
            # Incomplete-penetrance dominant gene: count homozygotes, not the
            # healthy heterozygous carriers that the standard AD path uses.
            if hom_only:
                if nhomalt >= hom_thr:
                    return self._met(
                        f"dominant (incomplete penetrance): "
                        f"nhomalt={nhomalt} >= {hom_thr}"
                    )
                return self._not_met(nhomalt, nhemi, het_carriers)
            if het_carriers >= het_thr:
                who = "healthy female carriers" if female_only else "healthy carriers"
                return self._met(f"dominant: {who}={het_carriers} >= {het_thr}")
        return self._not_met(nhomalt, nhemi, het_carriers)

    def _clinvar_fallback(
        self, variant: VariantRecord, gene: str | None
    ) -> CriteriaResult:
        """BS2 from an expert-panel (>=3-star) ClinVar review when the VCEP bars
        the gnomAD-count path. The strength is the one the panel cited (a bare
        "BS2" applies at its Strong default)."""
        from acmg_classifier.local_db.clinvar_sqlite import query_bs2_benign_evidence

        has_bs2, strength_label = query_bs2_benign_evidence(
            self._cfg.clinvar_sqlite,
            variant.chrom, variant.pos, variant.ref, variant.alt,
        )
        if not has_bs2:
            return CriteriaResult.not_met(
                ACMGCriterion.BS2,
                f"{gene}: VCEP bars gnomAD-based BS2 and no ClinVar "
                "expert-panel (>=3-star) BS2 found for this variant",
            )
        strength = _BS2_STRENGTH.get(strength_label or "Strong", CriterionStrength.STRONG)
        return CriteriaResult.met(
            ACMGCriterion.BS2,
            strength,
            evidence=(
                f"{gene}: ClinVar expert-panel (>=3-star) applied BS2"
                f"{'_' + strength_label if strength_label and strength_label != 'Strong' else ''} "
                "(VCEP internal-cohort evidence)"
            ),
        )

    @staticmethod
    def _observed_count(
        modes, nhomalt: int, nhemi: int, het_carriers: int,
        hom_only: bool, female_only: bool,
    ) -> tuple[int, str]:
        """The single gnomAD count that the gene's inheritance mode scores BS2
        on, with a human-readable label. Mirrors the mode selection of the
        flat-threshold path: AR→homozygotes, XL→hemizygotes, AD→carriers (or
        homozygotes for incomplete-penetrance hom_only genes); no mode →
        homozygotes."""
        if "AR" in modes:
            return nhomalt, "recessive: nhomalt"
        if "XL" in modes:
            return nhemi, "X-linked: nhemi"
        if "AD" in modes:
            if hom_only:
                return nhomalt, "dominant (incomplete penetrance): nhomalt"
            who = "healthy female carriers" if female_only else "healthy carriers"
            return het_carriers, f"dominant: {who}"
        return nhomalt, "nhomalt"

    def _met(self, detail: str, strength=None) -> CriteriaResult:
        if strength is not None:
            return CriteriaResult.met(
                ACMGCriterion.BS2, strength=strength, evidence=f"gnomAD {detail}"
            )
        return CriteriaResult.met(ACMGCriterion.BS2, evidence=f"gnomAD {detail}")

    def _not_met(self, nhomalt: int, nhemi: int, het: int) -> CriteriaResult:
        return CriteriaResult.not_met(
            ACMGCriterion.BS2,
            f"gnomAD nhomalt={nhomalt}, nhemi={nhemi}, carriers={het} (below threshold)",
        )
