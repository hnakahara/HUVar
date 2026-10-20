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

def _lof_thresholds(cfg) -> dict:
    """Config.pvs1_min_plp / pvs1_min_null_fraction (defaults when unset)."""
    from acmg_classifier.pvs1.transcript_evaluator import MIN_NULL_FRACTION, MIN_PLP
    mp = getattr(cfg, "pvs1_min_plp", MIN_PLP)
    mf = getattr(cfg, "pvs1_min_null_fraction", MIN_NULL_FRACTION)
    return {"min_plp": mp if isinstance(mp, int) and not isinstance(mp, bool) else MIN_PLP,
            "min_null_fraction": float(mf) if isinstance(mf, (int, float))
            and not isinstance(mf, bool) else MIN_NULL_FRACTION}


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
    ``None`` (default) uses ClinGen HI score 3 or >= 20% null among >= 5 ClinVar P/LP variants.

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
    # haploinsufficiency score of 3, or (iii) null variants make up >= 20% of
    # >= 5 ClinVar P/LP variants (formerly >= 3 P/LP null variants). gnomAD
    # LOEUF is reported for reference only. When none holds,
    # PVS1 is not applied.
    if lof_established:
        lof_mechanism, lof_note = True, "VCEP applies PVS1 (LoF established)"
    else:
        from acmg_classifier.local_db.clinvar_sqlite import query_pathogenic_null_stats
        from acmg_classifier.local_db.clingen_gene_db import hi_score
        plp_null, plp_total = query_pathogenic_null_stats(cfg.clinvar_sqlite, pc.gene_symbol)
        hi = hi_score(getattr(cfg, "clingen_dosage_tsv", None), pc.gene_symbol)
        lof_mechanism, lof_note = lof_mechanism_reason(
            loeuf, plp_null, hi, plp_total, **_lof_thresholds(cfg))
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
    # exon): a truncated protein may be produced. ClinGen SVI tree: if the
    # truncated/altered region is critical to protein function -> Strong
    # (annotated functional-domain overlap is the proxy); otherwise ("role of
    # region unknown") the protein-length branch applies: the variant removes
    # >10% of the protein -> Strong, <=10% -> Moderate. The altered region is
    # counted from the first changed residue (frameshifts alter everything
    # downstream). Without a known protein length PVS1 is withheld.
    if pc.domains:
        return CriterionStrength.STRONG, f"{cq}; {nmd_note}; truncated region contains functional domain"
    from acmg_classifier.pvs1.nmd_predictor import transcript_structure
    st = transcript_structure(pc.transcript_id)
    if st is None or pc.protein_position is None:
        return CriterionStrength.NOT_MET, (
            f"{cq}; {nmd_note}; protein length unknown — PVS1 N/A")
    aa_len = (st.cds_end - st.cds_start + 1) // 3 - 1
    removed = max(0, aa_len - pc.protein_position + 1)
    frac = removed / aa_len if aa_len > 0 else 0.0
    strength = CriterionStrength.STRONG if frac > 0.10 else CriterionStrength.MODERATE
    return strength, (
        f"{cq}; {nmd_note}; region role unknown, {removed}/{aa_len} aa "
        f"({frac:.1%}) altered -> {'>' if frac > 0.10 else '<='}10% (SVI)")


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

    The variant is in a canonical splice donor/acceptor (+/-1,2) by VEP
    consequence. Following the ClinGen SVI PVS1 decision tree (Abou Tayoun
    2018), disruption of splicing is presumed for these positions; a splice
    predictor is not required to confirm it. The predictor score is reported
    in the evidence text for reference only. (Previously a score below 0.20
    down-graded PVS1 to Strong/Moderate, which with the short-context
    OpenSpliceAI default model affected about a third of canonical sites.)"""
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

    # SVI: canonical +/-1,2 variants are presumed to disrupt splicing; the
    # predictor score is informational only.
    splice_tool_note = f"canonical splice site, splicing disruption presumed; {splice_tool_note}"
    splice_lof_predicted = True

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
