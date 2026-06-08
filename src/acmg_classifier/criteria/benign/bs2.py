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
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


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
        gd = annotation.gnomad
        if gd is None or not gd.filter_pass:
            return CriteriaResult.not_met(ACMGCriterion.BS2, "No valid gnomAD record")

        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None

        # VCEP gate: a VCEP that declined BS2 or barred population data for the
        # gene withholds it (e.g. RASopathy GN004 — variable expressivity).
        if self._vcep.status(gene) == NOT_APPLICABLE:
            return CriteriaResult.not_met(
                ACMGCriterion.BS2, f"{gene}: VCEP designates BS2 not applicable"
            )

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

        modes = self._vcep.modes(gene)
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
        if "AD" in modes and het_carriers >= het_thr:
            who = "healthy female carriers" if female_only else "healthy carriers"
            return self._met(f"dominant: {who}={het_carriers} >= {het_thr}")
        return self._not_met(nhomalt, nhemi, het_carriers)

    def _met(self, detail: str) -> CriteriaResult:
        return CriteriaResult.met(ACMGCriterion.BS2, evidence=f"gnomAD {detail}")

    def _not_met(self, nhomalt: int, nhemi: int, het: int) -> CriteriaResult:
        return CriteriaResult.not_met(
            ACMGCriterion.BS2,
            f"gnomAD nhomalt={nhomalt}, nhemi={nhemi}, carriers={het} (below threshold)",
        )
