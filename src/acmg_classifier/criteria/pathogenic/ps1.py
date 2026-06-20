"""PS1 — same amino acid change as established pathogenic variant."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


_RANK = {
    CriterionStrength.SUPPORTING: 1,
    CriterionStrength.MODERATE: 2,
    CriterionStrength.STRONG: 3,
}


def _ps1_strength(hits) -> CriterionStrength:
    """PS1 strength from the comparator's ClinVar classification (ClinGen SVI).

    Strong when at least one same-change comparator is classified Pathogenic
    (or the P-containing aggregate 'Pathogenic/Likely pathogenic'); Moderate
    when every comparator is only Likely pathogenic. Without this tiering PS1
    always fired at its Strong default, over-weighting LP-only comparators.
    """
    has_pathogenic = any(
        "pathogenic" in (h.clinical_significance or "").lower()
        and (h.clinical_significance or "").lower() != "likely pathogenic"
        for h in hits
    )
    return CriterionStrength.STRONG if has_pathogenic else CriterionStrength.MODERATE


def _cap(strength: CriterionStrength, ceiling: CriterionStrength | None) -> CriterionStrength:
    """Clamp a PS1 strength to the gene's VCEP ceiling (e.g. RMRP → Supporting)."""
    if ceiling is None or _RANK[strength] <= _RANK[ceiling]:
        return strength
    return ceiling


class PS1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Per-gene PS1 spec: some VCEPs restrict the splice extension to
        # non-canonical splice nucleotides (InSiGHT MMR genes).
        from acmg_classifier.criteria.ps1_genes import PS1Spec
        self._spec = PS1Spec(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # PS1 has two routes:
        #  * Missense -> the amino-acid rule: a DIFFERENT nucleotide change
        #    producing the SAME amino-acid substitution as a known pathogenic
        #    variant (matched on hgvs_p).
        #  * Splice / intronic -> the ClinGen SVI splicing extension: these have
        #    no protein change, so the amino-acid rule can never fire. Instead a
        #    DIFFERENT nucleotide change at the SAME splice-site position is
        #    recognised as having the same predicted splicing effect (matched on
        #    genomic position).
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.PS1, "No primary consequence")

        # VCEP gate: a VCEP that declined PS1 for the gene withholds it (CDH1).
        if self._spec.is_not_applicable(pc.gene_symbol):
            return CriteriaResult.not_met(
                ACMGCriterion.PS1, f"{pc.gene_symbol}: VCEP designates PS1 not applicable"
            )

        if pc.consequence == ConsequenceType.MISSENSE:
            return self._evaluate_missense(variant, pc)
        if pc.consequence in (
            ConsequenceType.SPLICE_ACCEPTOR,
            ConsequenceType.SPLICE_DONOR,
            ConsequenceType.SPLICE_REGION,
            ConsequenceType.INTRON,
        ):
            return self._evaluate_splice(variant, pc)
        return CriteriaResult.not_met(ACMGCriterion.PS1, "Not a missense or splice variant")

    def _evaluate_missense(self, variant: VariantRecord, pc) -> CriteriaResult:
        # Same amino-acid change at the same codon via a DIFFERENT nucleotide.
        # Exclude the variant itself by genomic coordinate so it cannot "match
        # itself". min_stars=1: the prior assertion must be a reviewed submission.
        from acmg_classifier.local_db.clinvar_sqlite import query_same_aa_change
        hits = query_same_aa_change(
            self._cfg.clinvar_sqlite,
            pc.gene_symbol,
            pc.hgvs_p,
            exclude_chrom=variant.chrom,
            exclude_pos=variant.pos,
            exclude_ref=variant.ref,
            exclude_alt=variant.alt,
            min_stars=1,
        )
        if hits:
            strength = _cap(_ps1_strength(hits), self._spec.max_strength(pc.gene_symbol))
            evidence = (
                f"ClinVar same AA (different nucleotide), comparator {strength.value}: "
                f"{', '.join(h.variation_id or '' for h in hits[:3])}"
            )
            return CriteriaResult.met(ACMGCriterion.PS1, strength, evidence)

        # Paralogue (analogous-residue) route: a VCEP may grant PS1 from the same
        # amino-acid change at the ANALOGOUS (same-numbered) residue of a sibling
        # gene — RASopathy groups (HRAS/NRAS/KRAS, MAP2K1/MAP2K2, SOS1/SOS2) and
        # HBA2/HBA1. Only consulted when the same-gene rule did not fire.
        siblings = self._spec.paralog_group(pc.gene_symbol)
        if siblings:
            # No self-exclusion / codon-proximity guard for paralogues: the hit is
            # in a DIFFERENT gene (different chromosome), so it can never be the
            # variant itself and the same-codon proximity rule does not apply.
            para = []
            for sib in siblings:
                para.extend(query_same_aa_change(
                    self._cfg.clinvar_sqlite, sib, pc.hgvs_p, min_stars=1,
                ))
            if para:
                fixed = self._spec.paralog_strength(pc.gene_symbol)
                strength = fixed if fixed is not None else _cap(
                    _ps1_strength(para), self._spec.max_strength(pc.gene_symbol)
                )
                evidence = (
                    f"ClinVar same AA in paralogue gene ({', '.join(siblings)}), "
                    f"{strength.value}: {', '.join(h.variation_id or '' for h in para[:3])}"
                )
                return CriteriaResult.met(ACMGCriterion.PS1, strength, evidence)

        return CriteriaResult.not_met(ACMGCriterion.PS1, "No ClinVar >=1 star same-AA hit (excluding self)")

    def _evaluate_splice(self, variant: VariantRecord, pc) -> CriteriaResult:
        # PS1's splice extension is opt-in per VCEP. Genes whose PS1 is the
        # original missense-only ACMG rule (no splice extension — e.g. GAA, the
        # HCM genes) must NOT receive PS1 for a splice/intronic variant.
        mode = self._spec.splice_mode(pc.gene_symbol)
        if not mode:
            return CriteriaResult.not_met(
                ACMGCriterion.PS1,
                f"{pc.gene_symbol}: PS1 is missense-only (no splice extension)",
            )
        # A "noncanonical" extension excludes canonical ±1/±2 sites
        # (SPLICE_DONOR / SPLICE_ACCEPTOR) — those are PVS1 territory.
        if mode == "noncanonical" and pc.consequence in (
            ConsequenceType.SPLICE_ACCEPTOR,
            ConsequenceType.SPLICE_DONOR,
        ):
            return CriteriaResult.not_met(
                ACMGCriterion.PS1,
                f"{pc.gene_symbol}: PS1 splice limited to non-canonical sites "
                "(canonical ±1/±2 is PVS1)",
            )

        # A different nucleotide change at the SAME splice-site position as a
        # known P/LP variant — the splicing counterpart of PS1.
        from acmg_classifier.local_db.clinvar_sqlite import query_same_splice_site
        hits = query_same_splice_site(
            self._cfg.clinvar_sqlite,
            pc.gene_symbol,
            variant.chrom, variant.pos, variant.ref, variant.alt,
            min_stars=1,
        )
        if not hits:
            return CriteriaResult.not_met(
                ACMGCriterion.PS1, "No ClinVar >=1 star same-splice-site P/LP hit"
            )
        strength = _cap(_ps1_strength(hits), self._spec.max_strength(pc.gene_symbol))
        evidence = (
            f"ClinVar same splice-site position (different nucleotide), comparator {strength.value}: "
            + ", ".join(h.variation_id or "" for h in hits[:3])
        )
        return CriteriaResult.met(ACMGCriterion.PS1, strength, evidence)
