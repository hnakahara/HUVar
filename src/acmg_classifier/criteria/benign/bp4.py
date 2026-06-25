"""BP4 -- computational evidence suggesting no impact (Bergquist 2024 thresholds)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion, ConsequenceType, CriterionStrength, InSilicoTool,
)
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry


def _alphamissense_bp4(score: float) -> CriterionStrength | None:
    """AlphaMissense BP4 strength per Bergquist 2024 Table 2.

    No Strong (-4) category exists for AlphaMissense; the strongest BP4 the
    table assigns is ThreePoint at score ≤ 0.070.
    """
    if score <= 0.070:
        return CriterionStrength.THREE_POINT
    if score <= 0.099:
        return CriterionStrength.MODERATE
    if score <= 0.169:
        return CriterionStrength.SUPPORTING
    return None


# REVEL genome-wide BP4 thresholds — Bergquist et al. 2024 Table 2 (same source
# and tier scheme as AlphaMissense/ESM1b above, including the ThreePoint (−3)
# tier). Lower REVEL ⇒ more benign; REVEL reaches Strong but has no Very Strong
# tier. Overridden per-gene when a VCEP states its own cutoff (revel_genes.py).
_REVEL_BP4_DEFAULT: dict[CriterionStrength, float] = {
    CriterionStrength.STRONG: 0.016,
    CriterionStrength.THREE_POINT: 0.052,
    CriterionStrength.MODERATE: 0.183,
    CriterionStrength.SUPPORTING: 0.290,
}


def _revel_bp4(score: float, tiers: dict[CriterionStrength, float] | None) -> CriterionStrength | None:
    """Strongest BP4 strength whose ``REVEL <= cutoff`` holds.

    Iterates tiers by ascending cutoff (the smallest cutoff is the strongest
    benign tier), so the first satisfied tier is the strongest one met.
    ``tiers`` is the gene's VCEP cutoff map when one exists, else None → the
    Pejaver defaults; a VCEP that only grants Supporting caps the gene there."""
    table = tiers if tiers else _REVEL_BP4_DEFAULT
    for strength, cutoff in sorted(table.items(), key=lambda kv: kv[1]):
        if score <= cutoff:
            return strength
    return None


def _esm1b_bp4(llr: float) -> CriterionStrength | None:
    """Bergquist 2024 Table 2 ESM1b BP4 thresholds (higher LLR ⇒ more benign).

    No Strong (-4) category; strongest BP4 is ThreePoint at LLR ≥ 8.8.
    """
    if llr >= 8.8:
        return CriterionStrength.THREE_POINT
    if llr >= -3.2:
        return CriterionStrength.MODERATE
    if llr >= -6.3:
        return CriterionStrength.SUPPORTING
    return None


def _spliceai_bp4(max_delta: float, cutoff: float = 0.10) -> CriterionStrength | None:
    """SpliceAI BP4 per Walker 2023 (no-impact cutoff, default 0.10; some VCEPs
    tighten/loosen it per gene)."""
    if max_delta <= cutoff:
        return CriterionStrength.SUPPORTING
    return None


def _openspliceai_bp4(max_delta: float, cutoff: float = 0.10) -> CriterionStrength | None:
    """OpenSpliceAI BP4. Same 0–1 delta scale as SpliceAI, so the Walker 2023
    <= 0.10 "no impact" cutoff (per-gene overridable) applies; awarded as Supporting."""
    if max_delta <= cutoff:
        return CriterionStrength.SUPPORTING
    return None


def _squirls_bp4(score: float) -> CriterionStrength | None:
    """SQUIRLS BP4: score < 0.50 (i.e. below the PP3 Supporting threshold)."""
    if score < 0.50:
        return CriterionStrength.SUPPORTING
    return None


# MMSplice DISABLED (dependency conflict). Retained, commented out, for later:
# def _mmsplice_bp4(delta_logit_psi: float) -> CriterionStrength | None:
#     """MMSplice BP4: |delta_logit_psi| < 0.5 (minimal predicted splice effect).
#
#     Well below the |delta_logit_psi| ≥ 2 PP3 (splice-impact) threshold, so a
#     near-zero score is conservative evidence of no splice effect."""
#     if abs(delta_logit_psi) < 0.5:
#         return CriterionStrength.SUPPORTING
#     return None


class BP4Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Gene-specific REVEL cutoffs (only consulted when insilico_tool=revel).
        from acmg_classifier.criteria.revel_genes import RevelSpec
        from acmg_classifier.criteria.bp_genes import BPApplicability
        self._revel_spec = RevelSpec(cfg.disease_prevalence_tsv)
        # Per-gene SpliceAI no-impact cutoff for the BP4 splice branch.
        self._bp_spec = BPApplicability(cfg.disease_prevalence_tsv)
        # Per-gene auxiliary BayesDel/CADD/2-of-3 rules (opt-in, licence-gated to
        # REVEL/AlphaMissense — see insilico_genes.combo_active). TP53 uses the
        # VCEP's precomputed PP3/BP4 code table (aGVGD baked in).
        from acmg_classifier.criteria.insilico_genes import InSilicoGeneSpec, TP53Codes
        self._insilico_spec = InSilicoGeneSpec(TP53Codes(cfg.tp53_codes_tsv))

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # BP4 mirrors PP3: same predictors and the same one-tool-only rule
        # for protein-level missense scoring (to prevent dual-counting
        # tools that share training data). The strength functions above use
        # the benign-side Bergquist 2024 thresholds.
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.BP4, "No consequence")

        if pc.consequence == ConsequenceType.MISSENSE:
            sp = annotation.splice
            # A missense with predicted splice impact (SpliceAI ≥ 0.20)
            # cannot be claimed as benign by computation alone — the splice
            # mechanism overrides the protein-level prediction.
            if sp and sp.is_available and sp.tool in ("spliceai", "openspliceai") and sp.max_delta is not None:
                if sp.max_delta >= 0.20:
                    tool_label = "SpliceAI" if sp.tool == "spliceai" else "OpenSpliceAI"
                    return CriteriaResult.not_met(
                        ACMGCriterion.BP4,
                        f"{tool_label} max_delta={sp.max_delta:.3f} — predicted splice impact, BP4 not applicable",
                    )

            # Per-gene auxiliary rule (BayesDel / REVEL∧CADD / 2-of-3). When
            # active it is AUTHORITATIVE for the gene and replaces the single-tool
            # dispatch below — so e.g. CTLA4 cannot meet BP4 on REVEL alone when
            # its VCEP requires REVEL∧CADD agreement. A None outcome means the rule
            # does not govern this consequence → fall through to the default.
            from acmg_classifier.criteria.insilico_genes import build_scores, combo_active
            gene = pc.gene_symbol
            if combo_active(self._insilico_spec, gene, self._cfg):
                outcome = self._insilico_spec.bp4(gene, build_scores(annotation, pc))
                if outcome is not None:
                    strength, note = outcome
                    if strength:
                        return CriteriaResult.met(ACMGCriterion.BP4, strength, note)
                    return CriteriaResult.not_met(ACMGCriterion.BP4, note)

            if self._cfg.insilico_tool == InSilicoTool.REVEL:
                rv = annotation.revel
                if rv and rv.score is not None:
                    rule = self._revel_spec.get(pc.gene_symbol)
                    tiers = rule.bp4 if (rule and rule.bp4) else None
                    strength = _revel_bp4(rv.score, tiers)
                    if strength:
                        src = f" [{pc.gene_symbol} VCEP cutoff]" if tiers else ""
                        return CriteriaResult.met(
                            ACMGCriterion.BP4, strength,
                            f"REVEL={rv.score:.3f} ({strength.value}){src}",
                        )
                    return CriteriaResult.not_met(
                        ACMGCriterion.BP4,
                        f"REVEL={rv.score:.3f} (not in benign range)",
                    )
                return CriteriaResult.not_met(ACMGCriterion.BP4, "No in-silico score available")

            if self._cfg.insilico_tool == InSilicoTool.ESM1B:
                es = annotation.esm1b
                if es and es.llr is not None:
                    strength = _esm1b_bp4(es.llr)
                    if strength:
                        return CriteriaResult.met(
                            ACMGCriterion.BP4, strength,
                            f"ESM1b LLR={es.llr:.3f} ({strength.value})",
                        )
                    return CriteriaResult.not_met(
                        ACMGCriterion.BP4,
                        f"ESM1b LLR={es.llr:.3f} (not in benign range)",
                    )
                return CriteriaResult.not_met(ACMGCriterion.BP4, "No in-silico score available")

            am = annotation.alphamissense
            if am and am.score is not None:
                strength = _alphamissense_bp4(am.score)
                if strength:
                    return CriteriaResult.met(
                        ACMGCriterion.BP4, strength,
                        f"AlphaMissense={am.score:.3f} ({strength.value})",
                    )
                return CriteriaResult.not_met(
                    ACMGCriterion.BP4,
                    f"AlphaMissense={am.score:.3f} (not in benign range)",
                )

        if pc.consequence in (
            ConsequenceType.SPLICE_REGION,
            ConsequenceType.INTRON,
            ConsequenceType.SYNONYMOUS,
        ):
            # Gene auxiliary CADD path for synonymous (e.g. ABCA4, BMPR2). Only a
            # *met* outcome short-circuits — a non-met must not suppress the
            # splice-based BP4 below, so we fall through otherwise.
            from acmg_classifier.criteria.insilico_genes import build_scores, combo_active
            if combo_active(self._insilico_spec, pc.gene_symbol, self._cfg):
                outcome = self._insilico_spec.bp4(pc.gene_symbol, build_scores(annotation, pc))
                if outcome is not None and outcome[0] is not None:
                    return CriteriaResult.met(ACMGCriterion.BP4, outcome[0], outcome[1])

            sp = annotation.splice
            # No splice predictor (e.g. default --splice-tool none) → say so
            # explicitly rather than implying a score was computed and rejected.
            if not (sp and sp.is_available):
                return CriteriaResult.not_met(
                    ACMGCriterion.BP4, "No splice prediction (splice evaluation disabled)",
                )
            cutoff = self._bp_spec.bp4_splice_cutoff(pc.gene_symbol)
            gene_cut = cutoff if cutoff is not None else 0.10
            cut_src = f" [{pc.gene_symbol} VCEP cutoff {gene_cut}]" if cutoff is not None else ""
            if sp.tool == "spliceai" and sp.max_delta is not None:
                strength = _spliceai_bp4(sp.max_delta, gene_cut)
                score_str = f"SpliceAI max_delta={sp.max_delta:.3f}{cut_src}"
            elif sp.tool == "openspliceai" and sp.max_delta is not None:
                strength = _openspliceai_bp4(sp.max_delta, gene_cut)
                score_str = f"OpenSpliceAI max_delta={sp.max_delta:.3f}{cut_src}"
            elif sp.tool == "squirls" and sp.raw_score is not None:
                strength = _squirls_bp4(sp.raw_score)
                score_str = f"SQUIRLS={sp.raw_score:.3f}"
            # MMSplice DISABLED — retained, commented out, for later:
            # elif sp.tool == "mmsplice" and sp.raw_score is not None:
            #     strength = _mmsplice_bp4(sp.raw_score)
            #     score_str = f"MMSplice delta_logit_psi={sp.raw_score:.3f}"
            else:
                return CriteriaResult.not_met(ACMGCriterion.BP4, "Splice score unavailable")

            if strength:
                return CriteriaResult.met(
                    ACMGCriterion.BP4, strength,
                    f"{score_str} ({strength.value})",
                )
            return CriteriaResult.not_met(ACMGCriterion.BP4, "Splice score not in benign range")

        return CriteriaResult.not_met(ACMGCriterion.BP4, "Consequence not applicable for BP4")
