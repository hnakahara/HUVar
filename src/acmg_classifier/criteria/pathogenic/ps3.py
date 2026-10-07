"""
PS3 -- well-established functional studies show a damaging effect.

Evidence sources (in priority order):
1. Curator-supplied supplement (user curation or the eRepo supplement) — any strength.
2. PS3 applied by a ClinGen expert panel to this variant in ClinVar (>=3 stars),
   imported at the strength the panel stated (curated evidence, not text mining).
3. Text mining of other (non-expert-panel) ClinVar submissions describing a
   damaging wet-lab assay. Evidence is counted by DISTINCT cited PMIDs (so several
   submissions citing the same study count once) and is capped at Supporting:
   the assay validation required by Brnich et al. 2019 (positive/negative
   controls, OddsPath) cannot be assessed from submission text.
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


# Genes whose VCEP does not permit a text-mined PS3 from ClinVar:
#   * no PS3 code in the VCEP spec at all — PALB2, PDHA1, POLG;
#   * PS3 explicitly "not applicable ... for in vitro assays" (only a rare,
#     case-by-case variant-specific animal model qualifies, which a curator must
#     assert via the manual supplement) — CAPN3, ANO5.
# A manual supplement PS3 still applies (it takes precedence below); only the
# free-text ClinVar fallback is suppressed for these genes.
_PS3_NOT_APPLICABLE = frozenset({"PALB2", "PDHA1", "POLG", "CAPN3", "ANO5"})

# Genes whose VCEP caps PS3 at Supporting — the text-mined fallback must not reach
# Moderate (n>=3) for them. Transcribed from the cspec PS3 strength descriptors.
_PS3_MAX_SUPPORTING = frozenset({
    "ABCD1", "AIPL1", "ETHE1", "F8", "F9", "GALT", "GAMT", "GATM", "HBA2",
    "HBB", "RPE65", "RPGR", "SERPINC1", "SLC19A3", "VHL",
})


def _functional_strength(n: int) -> CriterionStrength | None:
    # Text-mined functional evidence is limited to Supporting regardless of the
    # number of publications (Brnich 2019 validation cannot be read from text).
    return CriterionStrength.SUPPORTING if n >= 1 else None


class PS3Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        from acmg_classifier.criteria.curated_sources import (
            SRC_TEXT_MINING, from_expert_panel, from_supplement,
        )
        # 1. Curator-supplied evidence, 2. expert-panel assertion in ClinVar.
        for r in (from_supplement(ACMGCriterion.PS3, supplement),
                  from_expert_panel(ACMGCriterion.PS3, self._cfg, variant)):
            if r is not None:
                return r

        # Gene gate: some VCEPs do not allow a text-mined PS3 at all (no PS3 code,
        # or PS3 not applicable for in vitro assays). Withhold the free-text
        # fallback for those genes (a manual supplement above still applies).
        pc = annotation.primary_consequence
        gene = pc.gene_symbol if pc else None
        if gene in _PS3_NOT_APPLICABLE:
            return CriteriaResult.not_met(
                ACMGCriterion.PS3,
                f"{gene}: VCEP does not permit a ClinVar-text PS3 (no PS3 / in vitro N/A)",
            )

        # 3. Text mining: distinct PMIDs cited by non-expert-panel SCVs that
        #    describe a damaging functional assay (Supporting cap).
        from acmg_classifier.local_db.clinvar_sqlite import (
            query_functional_evidence, query_functional_pmids,
        )
        pmids = query_functional_pmids(
            self._cfg.clinvar_sqlite, variant.chrom, variant.pos, variant.ref, variant.alt,
        )
        if pmids is None:  # old ClinVar build without PMIDs -> SCV count
            n = query_functional_evidence(
                self._cfg.clinvar_sqlite, variant.chrom, variant.pos, variant.ref, variant.alt,
            )
            unit = "ClinVar SCV(s)"
        else:
            n = len(pmids)
            unit = "distinct cited publication(s)" + (f" (PMID {', '.join(sorted(pmids)[:5])})" if pmids else "")
        strength = _functional_strength(n)
        if strength is None:
            return CriteriaResult.not_met(
                ACMGCriterion.PS3, "No ClinVar submission citing a damaging functional study"
            )
        return CriteriaResult.met(
            ACMGCriterion.PS3,
            strength,
            evidence=f"{SRC_TEXT_MINING} {n} {unit} report a damaging functional study "
                     "(capped at Supporting)",
        )
