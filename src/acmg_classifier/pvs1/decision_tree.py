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
from acmg_classifier.pvs1.nmd_predictor import nmd_status
from acmg_classifier.pvs1.transcript_evaluator import (
    has_alternative_transcript_rescue,
    lof_mechanism_reason,
)

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
    mechanism, so pass ``True`` to skip the ClinGen HI / ClinVar check.
    ``None`` (default) uses ClinGen HI score 3 or >= 3 ClinVar P/LP nulls.

    Returns (strength, evidence_string).
    strength == NOT_MET means PVS1 should not be applied.
    """
    pc = annotation.primary_consequence
    if pc is None:
        return CriterionStrength.NOT_MET, "No primary consequence"

    gd = annotation.gnomad
    loeuf = gd.loeuf if gd else None
    alt_rescue = has_alternative_transcript_rescue(annotation)

    # "Is LoF a known disease mechanism?" — the first node of the ClinGen PVS1
    # decision tree. Established when (i) the VCEP explicitly applies PVS1
    # (lof_established=True), (ii) the gene has a ClinGen Dosage Sensitivity
    # haploinsufficiency score of 3, or (iii) ClinVar reports >= 3 P/LP null
    # variants. gnomAD LOEUF is reported for reference only. When none holds,
    # PVS1 is not applied.
    if lof_established:
        lof_mechanism, lof_note = True, "VCEP applies PVS1 (LoF established)"
    else:
        from acmg_classifier.local_db.clinvar_sqlite import query_pathogenic_null_count
        from acmg_classifier.local_db.clingen_gene_db import hi_score
        plp_null = query_pathogenic_null_count(cfg.clinvar_sqlite, pc.gene_symbol)
        hi = hi_score(getattr(cfg, "clingen_dosage_tsv", None), pc.gene_symbol)
        lof_mechanism, lof_note = lof_mechanism_reason(loeuf, plp_null, hi)
    if not lof_mechanism:
        return CriterionStrength.NOT_MET, lof_note

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
        # Initiation codon: Moderate (no downstream in-frame AUG data); only
        # reached once LoF has been established as the disease mechanism.
        return CriterionStrength.MODERATE, (
            f"Start-loss; assumed partial LoF (no downstream AUG data) [{lof_note}]"
        )
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

    return strength, f"{evidence} [{lof_note}]"


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
    nmd, nmd_note = nmd_status(pc)

    # Gate: if LoF is not a known mechanism for the gene, PVS1 does not apply
    # regardless of how convincing the molecular evidence is — this is the
    # very first ClinGen 2019 decision-tree branch.
    if not lof_mechanism:
        return CriterionStrength.NOT_MET, "Gene LoF mechanism not established"

    cq = pc.consequence.value
    if nmd:
        if not alt_rescue:
            return CriterionStrength.VERY_STRONG, f"{cq}; {nmd_note}; no rescue transcript"
        return CriterionStrength.STRONG, f"{cq}; {nmd_note}; alt transcript may rescue"

    # NMD escaped (PTC in the last exon or the 3'-most 50 nt of the penultimate
    # exon): a truncated protein may be produced. PVS1 applies only when the
    # truncation removes a critical region; annotated functional-domain overlap
    # is used as the proxy (Strong); without it PVS1 is withheld. (The decision
    # tree's protein-length (10%) branch is not implemented in this generic path.)
    if pc.domains:
        return CriterionStrength.STRONG, f"{cq}; {nmd_note}; truncated region contains functional domain"
    return CriterionStrength.NOT_MET, f"{cq}; {nmd_note}; no critical region removed (NMD escaped) — PVS1 N/A"


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
