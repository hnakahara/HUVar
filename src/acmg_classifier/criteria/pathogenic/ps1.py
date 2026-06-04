"""PS1 — same amino acid change as established pathogenic variant."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class PS1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

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
        if not hits:
            return CriteriaResult.not_met(ACMGCriterion.PS1, "No ClinVar >=1 star same-AA hit (excluding self)")
        evidence = f"ClinVar same AA (different nucleotide): {', '.join(h.variation_id or '' for h in hits[:3])}"
        return CriteriaResult.met(ACMGCriterion.PS1, evidence=evidence)

    def _evaluate_splice(self, variant: VariantRecord, pc) -> CriteriaResult:
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
        evidence = (
            "ClinVar same splice-site position (different nucleotide): "
            + ", ".join(h.variation_id or "" for h in hits[:3])
        )
        return CriteriaResult.met(ACMGCriterion.PS1, evidence=evidence)
