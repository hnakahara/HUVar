"""
ClinGen PVS1 2019 decision tree implementation.

Reference: Abou Tayoun et al. (2018) Hum Mutat 39:1517-1524
Flowchart: https://clinicalgenome.org/site/assets/files/3460/pvs1_decision_tree.pdf

Consequences handled:
  - Frameshift / stop-gained    -> NMD branch
  - Splice donor/acceptor       -> splice branch (integrates SQUIRLS/SpliceAI)
  - Start loss                  -> always Moderate (no NMD)
  - Transcript ablation         -> VeryStrong
"""
from __future__ import annotations

from acmg_classifier.config import Config
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.enums import ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pvs1.nmd_predictor import (
    predicts_nmd,
    is_last_exon,
    is_penultimate_exon,
)
from acmg_classifier.pvs1.transcript_evaluator import (
    has_alternative_transcript_rescue,
    gene_has_lof_mechanism,
)

# ClinGen SVI PVS1 strength caps (Franklin-aligned). When the gene has fewer
# than _MIN_PLP_NULL_FOR_FULL_PVS1 P/LP null variants in ClinVar, PVS1 strength
# is limited to Moderate in either of two situations:
#   - missense-dominant:  >= _MIN_PLP_MISS_FOR_CAP P/LP missense (mechanism
#                         skews to missense, not haploinsufficiency).
#   - undercurated:       total P/LP (null + missense) <
#                         _MIN_PLP_TOTAL_FOR_FULL_PVS1, so there is too little
#                         clinical evidence to establish LoF as the mechanism
#                         even when LOEUF/Z indicate population-level constraint.
_MIN_PLP_NULL_FOR_FULL_PVS1 = 3
_MIN_PLP_MISS_FOR_CAP = 10
_MIN_PLP_TOTAL_FOR_FULL_PVS1 = 5


def evaluate_pvs1(
    variant: VariantRecord,
    annotation: AnnotationData,
    cfg: Config,
) -> tuple[CriterionStrength, str]:
    """
    Run the ClinGen 2019 PVS1 decision tree.

    Returns (strength, evidence_string).
    strength == NOT_MET means PVS1 should not be applied.
    """
    pc = annotation.primary_consequence
    if pc is None:
        return CriterionStrength.NOT_MET, "No primary consequence"

    gd = annotation.gnomad
    loeuf = gd.loeuf if gd else None
    alt_rescue = has_alternative_transcript_rescue(annotation)

    # "Is LoF a known disease mechanism?" — primary signal is the count of P/LP
    # null variants already reported in ClinVar for the gene (Franklin-style);
    # LOEUF is only a secondary constraint hint.
    from acmg_classifier.local_db.clinvar_sqlite import query_pathogenic_null_count
    plp_null = query_pathogenic_null_count(cfg.clinvar_sqlite, pc.gene_symbol)
    lof_mechanism = gene_has_lof_mechanism(pc, loeuf, plp_null)

    # ---- Branch dispatch (compute strength + evidence) ---------------------
    if pc.consequence == ConsequenceType.TRANSCRIPT_ABLATION:
        if lof_mechanism:
            strength, evidence = (
                CriterionStrength.VERY_STRONG,
                "Transcript ablation; LoF is an established disease mechanism",
            )
        else:
            strength, evidence = (
                CriterionStrength.STRONG,
                "Transcript ablation; LoF mechanism uncertain",
            )
    elif pc.consequence == ConsequenceType.START_LOST:
        # Already at Moderate — never subject to cap.
        return CriterionStrength.MODERATE, "Start-loss; assumed partial LoF (no downstream AUG data)"
    elif pc.consequence in (ConsequenceType.SPLICE_DONOR, ConsequenceType.SPLICE_ACCEPTOR):
        strength, evidence = _splice_branch(
            variant, annotation, cfg, lof_mechanism, alt_rescue, pc,
        )
    elif pc.consequence in (ConsequenceType.FRAMESHIFT, ConsequenceType.STOP_GAINED):
        strength, evidence = _nmd_branch(annotation, lof_mechanism, alt_rescue, pc)
    else:
        return (
            CriterionStrength.NOT_MET,
            f"Consequence {pc.consequence.value} not handled by PVS1",
        )

    # ---- ClinGen SVI strength caps -----------------------------------------
    # When the gene has few P/LP null variants, two caps may apply:
    #   (1) missense-dominant: many P/LP missense imply the disease mechanism is
    #       NOT haploinsufficiency, so a new null cannot be VeryStrong/Strong.
    #   (2) undercurated: very few P/LP overall — there is too little clinical
    #       evidence in ClinVar to establish LoF as the disease mechanism even
    #       when LOEUF/Z indicate population-level constraint.
    # Either situation caps PVS1 strength to Moderate (cf. Franklin's PVS1).
    if strength in (CriterionStrength.VERY_STRONG, CriterionStrength.STRONG) \
            and plp_null < _MIN_PLP_NULL_FOR_FULL_PVS1:
        from acmg_classifier.local_db.clinvar_sqlite import query_pathogenic_missense_count
        plp_miss = query_pathogenic_missense_count(cfg.clinvar_sqlite, pc.gene_symbol)
        plp_total = plp_null + plp_miss
        if plp_miss >= _MIN_PLP_MISS_FOR_CAP:
            evidence = (
                f"{evidence} [capped to Moderate: only {plp_null} P/LP null "
                f"but {plp_miss} P/LP missense — missense-dominant gene]"
            )
            strength = CriterionStrength.MODERATE
        elif plp_total < _MIN_PLP_TOTAL_FOR_FULL_PVS1:
            evidence = (
                f"{evidence} [capped to Moderate: only {plp_null} P/LP null "
                f"and {plp_miss} P/LP missense (<{_MIN_PLP_TOTAL_FOR_FULL_PVS1} "
                f"total) — insufficient ClinVar evidence for LoF mechanism]"
            )
            strength = CriterionStrength.MODERATE
    return strength, evidence


# ---------------------------------------------------------------------------
# NMD branch (frameshift / stop-gained)
# ---------------------------------------------------------------------------

def _nmd_branch(
    annotation: AnnotationData,
    lof_mechanism: bool,
    alt_rescue: bool,
    pc,
) -> tuple[CriterionStrength, str]:
    nmd = predicts_nmd(pc)

    if not lof_mechanism:
        return CriterionStrength.NOT_MET, "Gene LoF mechanism not established"

    if nmd:
        # NMD predicted
        if not alt_rescue:
            return CriterionStrength.VERY_STRONG, f"{pc.consequence.value}; NMD predicted; no rescue transcript"
        else:
            return CriterionStrength.STRONG, f"{pc.consequence.value}; NMD predicted; alt transcript may rescue"
    else:
        # NMD NOT predicted (last exon or penultimate)
        last = is_last_exon(pc)
        penult = is_penultimate_exon(pc)
        note = "last exon" if last else ("penultimate exon" if penult else "NMD escape")

        if last or penult:
            # Check if truncated region is critical
            domains = pc.domains or []
            has_domain = bool(domains)
            if has_domain:
                return CriterionStrength.STRONG, f"{pc.consequence.value}; {note}; truncated region contains functional domain"
            else:
                return CriterionStrength.MODERATE, f"{pc.consequence.value}; {note}; no critical domain data"
        else:
            return CriterionStrength.SUPPORTING, f"{pc.consequence.value}; NMD not predicted; uncertain impact"


# ---------------------------------------------------------------------------
# Splice branch
# ---------------------------------------------------------------------------

def _splice_branch(
    variant: VariantRecord,
    annotation: AnnotationData,
    cfg: Config,
    lof_mechanism: bool,
    alt_rescue: bool,
    pc,
) -> tuple[CriterionStrength, str]:
    if not lof_mechanism:
        return CriterionStrength.NOT_MET, "Gene LoF mechanism not established for splice variant"

    sp = annotation.splice
    splice_lof_predicted = False
    splice_tool_note = "no splice tool"

    if sp and sp.is_available:
        if sp.tool == "spliceai" and sp.max_delta is not None:
            splice_lof_predicted = sp.max_delta >= 0.20
            splice_tool_note = f"SpliceAI={sp.max_delta:.3f}"
        elif sp.tool == "squirls" and sp.raw_score is not None:
            splice_lof_predicted = sp.raw_score >= 0.50
            splice_tool_note = f"SQUIRLS={sp.raw_score:.3f} (approx)"

    if splice_lof_predicted:
        if not alt_rescue:
            return (
                CriterionStrength.VERY_STRONG,
                f"{pc.consequence.value}; splice LoF predicted ({splice_tool_note}); no rescue transcript",
            )
        else:
            return (
                CriterionStrength.STRONG,
                f"{pc.consequence.value}; splice LoF predicted ({splice_tool_note}); alt transcript may rescue",
            )
    else:
        # Splice impact uncertain or not predicted
        exon_skip_possible = True  # conservative assumption without RNA data
        if exon_skip_possible:
            domains = pc.domains or []
            if domains:
                return (
                    CriterionStrength.STRONG,
                    f"{pc.consequence.value}; {splice_tool_note}; exon skip affects functional domain",
                )
            return (
                CriterionStrength.MODERATE,
                f"{pc.consequence.value}; {splice_tool_note}; exon skip, domain unknown",
            )
        return (
            CriterionStrength.SUPPORTING,
            f"{pc.consequence.value}; splice impact uncertain ({splice_tool_note})",
        )
