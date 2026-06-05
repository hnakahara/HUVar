"""BP1 -- variant type unlikely to be pathogenic for the gene (VCEP-gated).

Original ACMG BP1: a missense variant in a gene where truncating variants are
the primary disease mechanism. ClinGen VCEPs make BP1 strongly gene-specific:

* Most VCEPs decline BP1 (``not_applicable``) — it is then never applied.
* Genes that apply it target a specific consequence:
  - ``missense`` (PALB2, APC, BRCA1/2): a missense variant is BP1.
  - ``truncating`` (gain-of-function RASopathy genes — BRAF, the RAS/MAPK genes,
    PTPN11, …): loss-of-function is benign, so a TRUNCATING variant is BP1.

Applicability and the target consequence are read from ``disease_prevalence.tsv``
(``bp1`` / ``bp1_target`` columns). A gene with no VCEP BP1 decision does not
receive BP1 (the criterion is too gene-dependent to apply by default).
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.bp_genes import BPApplicability, APPLICABLE
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_TRUNCATING = (
    ConsequenceType.STOP_GAINED,
    ConsequenceType.FRAMESHIFT,
)
# "broad" target (BRCA1/2): silent substitution, missense, or in-frame indel.
_BROAD = (
    ConsequenceType.SYNONYMOUS,
    ConsequenceType.MISSENSE,
    ConsequenceType.INFRAME_INSERTION,
    ConsequenceType.INFRAME_DELETION,
)
# SpliceAI max-delta ceiling for the BRCA "no splicing predicted" condition.
_NO_SPLICE_MAX = 0.10


class BP1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        self._vcep = BPApplicability(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None

        # BP1 is applied only where a VCEP explicitly designates it applicable.
        if self._vcep.bp1(gene) != APPLICABLE:
            return CriteriaResult.not_met(
                ACMGCriterion.BP1, "BP1 not applicable for this gene (no VCEP designation)"
            )
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.BP1, "No primary consequence")

        target = self._vcep.bp1_target(gene)

        # 1. Eligible consequence for the gene's BP1 target.
        if target == "truncating":
            if pc.consequence not in _TRUNCATING:
                return CriteriaResult.not_met(ACMGCriterion.BP1, "Not a truncating variant")
        elif target == "broad":
            if pc.consequence not in _BROAD:
                return CriteriaResult.not_met(
                    ACMGCriterion.BP1, "Not a silent / missense / in-frame variant"
                )
        else:  # missense (default)
            if pc.consequence != ConsequenceType.MISSENSE:
                return CriteriaResult.not_met(ACMGCriterion.BP1, "Not a missense variant")

        # 2. Region exclusion: BP1 does not apply inside the APC β-catenin repeat
        #    or a BRCA1/2 clinically-important functional domain.
        if self._vcep.bp1_excluded(gene, pc.protein_position):
            return CriteriaResult.not_met(
                ACMGCriterion.BP1,
                f"{gene}: in an excluded region / functional domain",
            )

        # 3. "No predicted splice impact" gate (BRCA1/2: SpliceAI <= 0.1). Without
        #    a splice prediction the condition cannot be confirmed, so BP1 is
        #    withheld (a splice-altering variant must not be called benign).
        if self._vcep.bp1_requires_no_splice(gene):
            sp = annotation.splice
            if sp is None or sp.max_delta is None:
                return CriteriaResult.not_met(
                    ACMGCriterion.BP1, f"{gene}: splice prediction unavailable (required)"
                )
            if sp.max_delta > _NO_SPLICE_MAX:
                return CriteriaResult.not_met(
                    ACMGCriterion.BP1,
                    f"{gene}: predicted splice impact (SpliceAI {sp.max_delta:.2f})",
                )

        strength = CriterionStrength.STRONG if self._vcep.bp1_is_strong(gene) else None
        return CriteriaResult.met(
            ACMGCriterion.BP1,
            strength=strength,
            evidence=f"{gene}: {pc.consequence.value} (VCEP BP1{', Strong' if strength else ''})",
        )
