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
        pc = annotation.primary_consequence
        if pc is None or pc.consequence not in (
            ConsequenceType.MISSENSE,
            ConsequenceType.SPLICE_ACCEPTOR,
            ConsequenceType.SPLICE_DONOR,
        ):
            return CriteriaResult.not_met(ACMGCriterion.PS1, "Not a missense or splice variant")

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
