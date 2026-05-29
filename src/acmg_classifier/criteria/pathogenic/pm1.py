"""PM1 — located in mutational hotspot or critical functional domain."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


class PM1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # PM1 only applies to variants that change protein sequence at a
        # specific codon — missense and small in-frame indels. Truncating
        # variants (frameshift, stop_gain) are PVS1 territory, and silent
        # variants do not move an amino-acid position into a hotspot.
        pc = annotation.primary_consequence
        if pc is None or pc.consequence not in (
            ConsequenceType.MISSENSE,
            ConsequenceType.INFRAME_INSERTION,
            ConsequenceType.INFRAME_DELETION,
        ):
            return CriteriaResult.not_met(ACMGCriterion.PM1, "Not a missense/in-frame variant")

        # Operationally we approximate "mutational hotspot or critical
        # functional domain" by counting nearby pathogenic ClinVar entries at
        # the same protein position / cluster window. The cluster definition
        # itself lives in clinvar_sqlite.query_hotspot_cluster.
        from acmg_classifier.local_db.clinvar_sqlite import query_hotspot_cluster
        is_hotspot, evidence = query_hotspot_cluster(
            self._cfg.clinvar_sqlite,
            pc.gene_symbol,
            pc.protein_position,
        )
        if not is_hotspot:
            return CriteriaResult.not_met(ACMGCriterion.PM1, "Not in hotspot cluster")
        return CriteriaResult.met(ACMGCriterion.PM1, evidence=evidence)
