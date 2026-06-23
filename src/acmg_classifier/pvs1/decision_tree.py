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
    lof_established: bool | None = None,
) -> tuple[CriterionStrength, str]:
    """
    Run the ClinGen 2019 PVS1 decision tree.

    ``lof_established`` overrides the LoF-mechanism question: a VCEP that
    explicitly applies PVS1 has, by definition, established LoF as the disease
    mechanism, so pass ``True`` to skip the ClinVar/LOEUF heuristic (which can
    miss under-represented genes). ``None`` (default) uses the heuristic.

    Returns (strength, evidence_string).
    strength == NOT_MET means PVS1 should not be applied.
    """
    pc = annotation.primary_consequence
    if pc is None:
        return CriterionStrength.NOT_MET, "No primary consequence"

    gd = annotation.gnomad
    loeuf = gd.loeuf if gd else None
    alt_rescue = has_alternative_transcript_rescue(annotation)

    # "Is LoF a known disease mechanism?" — when the VCEP explicitly applies PVS1
    # (lof_established=True) the mechanism is settled; otherwise the primary
    # signal is the count of P/LP null variants already reported in ClinVar for
    # the gene (Franklin-style), with LOEUF as a secondary constraint hint.
    if lof_established:
        lof_mechanism = True
    else:
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
    # Skipped when the VCEP explicitly applies PVS1 (lof_established) — the panel
    # has already established LoF as the mechanism, so the ClinVar-curation caps
    # do not apply.
    if not lof_established \
            and strength in (CriterionStrength.VERY_STRONG, CriterionStrength.STRONG) \
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
    """ClinGen 2019 PVS1 sub-tree for frameshift / stop-gained variants.

    Rationale: when NMD is predicted, the truncated mRNA is degraded so the
    allele effectively produces no protein → Very Strong. If an alternative
    transcript can rescue the LoF, we down-grade because the cell may
    still express a functional protein → Strong. When NMD is escaped, the
    truncated protein may still be expressed; severity then depends on what
    region is removed (functional-domain truncation is more damaging than
    truncation of an uncharacterised C-terminus)."""
    nmd = predicts_nmd(pc)

    # Gate: if LoF is not a known mechanism for the gene, PVS1 does not apply
    # regardless of how convincing the molecular evidence is — this is the
    # very first ClinGen 2019 decision-tree branch.
    if not lof_mechanism:
        return CriterionStrength.NOT_MET, "Gene LoF mechanism not established"

    if nmd:
        if not alt_rescue:
            return CriterionStrength.VERY_STRONG, f"{pc.consequence.value}; NMD predicted; no rescue transcript"
        else:
            return CriterionStrength.STRONG, f"{pc.consequence.value}; NMD predicted; alt transcript may rescue"
    else:
        # NMD is escaped when the premature stop is in the last exon or within
        # ~50 bp of the last exon-exon junction (penultimate exon). These two
        # cases are usually grouped because the rule of thumb cannot
        # distinguish them without splice-junction-level precision.
        last = is_last_exon(pc)
        penult = is_penultimate_exon(pc)
        note = "last exon" if last else ("penultimate exon" if penult else "NMD escape")

        if last or penult:
            # NMD escapes, so a (largely) full-length protein is still made.
            # Per the ClinGen SVI PVS1 decision tree, PVS1 then applies only when
            # the truncation removes a CRITICAL functional region; otherwise it
            # is N/A. Domain presence is our proxy for "critical region
            # truncated": with a functional domain in the truncated tail → Strong;
            # WITHOUT any domain evidence we must NOT assume criticality (the old
            # "Moderate" over-applied PVS1 to last-exon truncations the VCEPs
            # leave uncalled, e.g. APC/MYOC), so PVS1 is withheld.
            domains = pc.domains or []
            has_domain = bool(domains)
            if has_domain:
                return CriterionStrength.STRONG, f"{pc.consequence.value}; {note}; truncated region contains functional domain"
            else:
                return CriterionStrength.NOT_MET, f"{pc.consequence.value}; {note}; no critical region removed (NMD escaped) — PVS1 N/A"
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
    """ClinGen 2019 PVS1 sub-tree for canonical-splice variants.

    The variant is in a splice donor/acceptor by VEP consequence. We use
    a splice predictor score to decide whether the LoF interpretation is
    supported. When no predictor is available we still award a reduced
    strength because the canonical splice site itself is strong prior
    evidence of LoF — but we cap at Moderate to reflect the missing
    confirmation."""
    if not lof_mechanism:
        return CriterionStrength.NOT_MET, "Gene LoF mechanism not established for splice variant"

    sp = annotation.splice
    splice_lof_predicted = False
    splice_tool_note = "no splice tool"

    # Threshold differences: SpliceAI 0.20 is the Walker 2023 calibration.
    # SQUIRLS uses 0.50 (a higher bar) because its score distribution is
    # different and it is NOT Walker-calibrated; we tag the note "(approx)"
    # so reviewers can see this caveat in the evidence string.
    if sp and sp.is_available:
        if sp.tool == "spliceai" and sp.max_delta is not None:
            splice_lof_predicted = sp.max_delta >= 0.20
            splice_tool_note = f"SpliceAI={sp.max_delta:.3f}"
        elif sp.tool == "openspliceai" and sp.max_delta is not None:
            # Same 0–1 delta scale as SpliceAI; this is a LoF-prediction gate
            # (not a strength-tier calibration), so the same 0.20 cutoff applies.
            splice_lof_predicted = sp.max_delta >= 0.20
            splice_tool_note = f"OpenSpliceAI={sp.max_delta:.3f}"
        elif sp.tool == "squirls" and sp.raw_score is not None:
            splice_lof_predicted = sp.raw_score >= 0.50
            splice_tool_note = f"SQUIRLS={sp.raw_score:.3f} (approx)"
        # MMSplice DISABLED — retained, commented out, for later:
        # elif sp.tool == "mmsplice" and sp.raw_score is not None:
        #     # |delta_logit_psi| >= 2 → strong predicted splice effect (MMSplice 2019).
        #     splice_lof_predicted = abs(sp.raw_score) >= 2.0
        #     splice_tool_note = f"MMSplice delta_logit_psi={sp.raw_score:.3f}"

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
        # Splice predictor disagrees or unavailable. Without RNA-seq we cannot
        # rule out exon skipping, so we conservatively assume it is possible
        # and fall back to the domain-presence heuristic used in _nmd_branch.
        # NOTE: `exon_skip_possible` is hard-coded True — the Supporting
        # branch at the bottom is currently unreachable. See cleanup-candidates.md.
        exon_skip_possible = True
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
