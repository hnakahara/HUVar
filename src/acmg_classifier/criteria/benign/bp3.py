"""BP3 -- in-frame indel in repetitive region without known function."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class BP3Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        from acmg_classifier.criteria.bp_genes import BPApplicability
        self._vcep = BPApplicability(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # BP3 is the benign-side counterpart to PM4: in-frame indels in
        # repetitive regions are likely tolerated because length variation
        # in repeats is biologically common. PM4 explicitly hands these
        # variants to BP3 (see pm4.py).
        pc = annotation.primary_consequence
        if pc is None or pc.consequence not in (
            ConsequenceType.INFRAME_INSERTION,
            ConsequenceType.INFRAME_DELETION,
        ):
            return CriteriaResult.not_met(ACMGCriterion.BP3, "Not an in-frame indel")

        # VCEP gate: a VCEP that declined BP3 for the gene suppresses the repeat
        # heuristic (curbs the generic over-assignment). Genes with no VCEP
        # opinion keep the heuristic.
        from acmg_classifier.criteria.bp_genes import NOT_APPLICABLE
        if self._vcep.bp3(pc.gene_symbol) == NOT_APPLICABLE:
            return CriteriaResult.not_met(
                ACMGCriterion.BP3, f"{pc.gene_symbol}: VCEP designates BP3 not applicable"
            )

        # Region-restricted BP3 (RPGR ORF15 aa585-1078; FOXG1 poly-AA tracts):
        # the VCEP enumerates the repetitive regions itself, so being inside one
        # IS the "repetitive region" — award BP3 directly (no Dfam needed); a
        # variant outside them is not BP3.
        in_region = self._vcep.bp3_in_region(pc.gene_symbol, pc.protein_position)
        if in_region is not None:
            if in_region:
                return CriteriaResult.met(
                    ACMGCriterion.BP3,
                    evidence=f"{pc.gene_symbol}: in VCEP BP3 repetitive region "
                    f"(residue {pc.protein_position})",
                )
            return CriteriaResult.not_met(
                ACMGCriterion.BP3, f"{pc.gene_symbol}: outside the VCEP BP3 region"
            )

        # Repeat overlap is determined by Dfam/RepeatMasker (see
        # local_db.repeatmasker_db). The "without known function" qualifier
        # from the ACMG text is not enforced here — we assume the repeat
        # tracks already exclude domains with known function.
        rep = annotation.repeat
        if rep and rep.in_repeat:
            return CriteriaResult.met(
                ACMGCriterion.BP3,
                evidence=f"In-frame indel in repeat element ({rep.repeat_class}: {rep.repeat_name})",
            )
        return CriteriaResult.not_met(ACMGCriterion.BP3, "Not in a repeat region (Dfam)")
