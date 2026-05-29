"""PM5 -- different missense at same codon as established pathogenic variant."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class PM5Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # PM5 is the "different AA, same codon" companion to PS1's "same AA"
        # rule, so only missense is eligible — splice variants don't have a
        # comparable codon-level interpretation, and synonymous changes by
        # definition do not change the amino acid.
        pc = annotation.primary_consequence
        if pc is None or pc.consequence != ConsequenceType.MISSENSE:
            return CriteriaResult.not_met(ACMGCriterion.PM5, "Not a missense variant")

        # Look for ClinVar P/LP variants at the same protein position with a
        # DIFFERENT amino-acid substitution (pc.hgvs_p is passed so the query
        # can exclude exact-match variants — those are PS1, not PM5).
        # min_stars=1 enforces the ACMG requirement for a reviewed assertion.
        from acmg_classifier.local_db.clinvar_sqlite import query_same_codon_different_aa
        hits = query_same_codon_different_aa(
            self._cfg.clinvar_sqlite,
            pc.gene_symbol,
            pc.protein_position,
            pc.hgvs_p,
            min_stars=1,
        )
        if not hits:
            return CriteriaResult.not_met(
                ACMGCriterion.PM5, "No ClinVar >=1 star same-codon different-AA hit"
            )
        # Cap at 3 IDs so the evidence string stays human-readable in TSV/JSON.
        evidence = "ClinVar same codon, diff AA: " + ", ".join(
            h.variation_id or "" for h in hits[:3]
        )
        return CriteriaResult.met(ACMGCriterion.PM5, evidence=evidence)
