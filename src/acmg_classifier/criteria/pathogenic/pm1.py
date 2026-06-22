"""PM1 — located in mutational hotspot or critical functional domain."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.pm1_hotspots import PM1Hotspots
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


class PM1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Per-gene VCEP hotspot definitions (residue ranges / residues / strength)
        # mined from the cspec summaries. Authoritative where present.
        self._hotspots = PM1Hotspots(cfg.pm1_hotspots_tsv)

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

        gene = pc.gene_symbol

        # 1. VCEP declared PM1 not applicable for the gene (benign variation
        #    throughout / no defined hotspot — e.g. ABCA4, ATM, RASopathy genes).
        if self._hotspots.is_not_applicable(gene):
            return CriteriaResult.not_met(
                ACMGCriterion.PM1, f"{gene}: VCEP designates PM1 not applicable"
            )

        # 2. VCEP-curated hotspot regions are authoritative where defined: a hit
        #    awards PM1 at the VCEP strength; a miss withholds PM1 (do NOT fall
        #    back to the statistical heuristic, which the VCEP regions supersede).
        if self._hotspots.has_gene(gene) or self._hotspots.has_cys_creating(gene):
            strength = self._hotspots.lookup(gene, pc.protein_position)
            evidence = f"{gene} residue {pc.protein_position} in VCEP PM1 hotspot"
            # Cys-creating missense: a substitution introducing a new cysteine in a
            # disulfide-bonded domain (FBN1 EGF/cbEGF/TB/hybrid) earns PM1_Moderate
            # regardless of the curated residue list.
            if (pc.consequence == ConsequenceType.MISSENSE
                    and _alt_aa(pc.amino_acids) == "C"
                    and self._hotspots.in_cys_creating_region(gene, pc.protein_position)):
                if strength is None or _RANK[strength] < _RANK[CriterionStrength.MODERATE]:
                    strength = CriterionStrength.MODERATE
                    evidence = (f"{gene} Cys-creating missense at residue "
                                f"{pc.protein_position} in a disulfide-bonded domain")
            if strength is None:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM1, f"{gene}: not in a VCEP PM1 hotspot region"
                )
            return CriteriaResult.met(
                ACMGCriterion.PM1, strength=strength, evidence=evidence,
            )

        # 3. No curated VCEP data for the gene — fall back to the statistical
        #    hotspot heuristic (nearby pathogenic ClinVar clustering).
        from acmg_classifier.local_db.clinvar_sqlite import query_hotspot_cluster
        is_hotspot, evidence = query_hotspot_cluster(
            self._cfg.clinvar_sqlite,
            gene,
            pc.protein_position,
        )
        if not is_hotspot:
            return CriteriaResult.not_met(ACMGCriterion.PM1, "Not in hotspot cluster")
        return CriteriaResult.met(ACMGCriterion.PM1, evidence=evidence)


def _alt_aa(amino_acids: str | None) -> str | None:
    """The 1-letter alternate residue from a VEP ``amino_acids`` "X/Y" missense
    field, or None when absent / not a single-residue substitution."""
    if not amino_acids or "/" not in amino_acids:
        return None
    alt = amino_acids.partition("/")[2].strip()
    return alt if len(alt) == 1 else None
