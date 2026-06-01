"""PP3 -- computational evidence of deleterious effect (Bergquist 2024 thresholds)."""
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


def _alphamissense_pp3(score: float) -> CriterionStrength | None:
    """Map AlphaMissense (Cheng et al. 2023) to PP3 strength.

    Thresholds follow Bergquist et al. 2024 Table 2, which calibrated each
    in-silico tool to ACMG strength tiers via OddsPath. We return None below
    the Supporting cutoff so callers can distinguish "no evidence" from
    "explicitly benign-leaning" (the BP4 evaluator handles the latter side)."""
    if score >= 0.990:
        return CriterionStrength.STRONG
    if score >= 0.972:
        return CriterionStrength.THREE_POINT
    if score >= 0.906:
        return CriterionStrength.MODERATE
    if score >= 0.792:
        return CriterionStrength.SUPPORTING
    return None


def _esm1b_pp3(llr: float) -> CriterionStrength | None:
    """Bergquist 2024 Table 2 ESM1b PP3 thresholds.

    LLR convention (Brandes et al., Nat Genet 2023): more negative LLR
    indicates more pathogenic. Comparisons use <= because the threshold is
    a *minimum* magnitude on the negative side."""
    if llr <= -24.0:
        return CriterionStrength.STRONG
    if llr <= -14.0:
        return CriterionStrength.THREE_POINT
    if llr <= -12.2:
        return CriterionStrength.MODERATE
    if llr <= -10.7:
        return CriterionStrength.SUPPORTING
    return None


def _squirls_pp3(score: float) -> CriterionStrength | None:
    """SQUIRLS (Danis 2021) splice-pathogenicity threshold.

    Only ≥ 0.50 → Supporting. SQUIRLS lacks a published OddsPath calibration
    for ACMG Moderate/Strong, so only one conservative tier is applied."""
    if score >= 0.50:
        return CriterionStrength.SUPPORTING
    return None


def _mmsplice_pp3(delta_logit_psi: float) -> CriterionStrength | None:
    """MMSplice (Cheng et al., Genome Biology 2019) splice-effect threshold.

    delta_logit_psi is on a logit scale; |delta_logit_psi| > 2 is a strong
    splice effect per the MMSplice paper. We judge by absolute value because
    both strong exclusion (negative) and strong inclusion shifts (positive)
    are aberrant splicing. No published OddsPath calibration to ACMG tiers
    exists, so only the conservative Supporting tier is awarded."""
    if abs(delta_logit_psi) >= 2.0:
        return CriterionStrength.SUPPORTING
    return None


def _spliceai_pp3(max_delta: float) -> CriterionStrength | None:
    """SpliceAI (Jaganathan 2019) PP3 cutoff.

    Walker 2023 (ClinGen SVI splicing WG) recommends max_delta >= 0.20 as
    Moderate for predicted splice impact. Stronger tiers require additional
    RNA evidence and are intentionally not awarded by the predictor alone."""
    if max_delta >= 0.20:
        return CriterionStrength.MODERATE
    return None


class PP3Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.PP3, "No consequence")

        # Missense branch: the variant changes the amino acid but might also
        # affect splicing. ClinGen SVI splicing-WG says we should award PP3
        # for predicted splice impact even on missense variants because the
        # mechanism (loss of normal protein via aberrant splicing) is
        # independent of the protein-level damage signal.
        if pc.consequence == ConsequenceType.MISSENSE:
            sp = annotation.splice
            # Only SpliceAI has a Walker-calibrated cutoff for this scenario;
            # SQUIRLS is not used here to avoid uncertain dual-counting.
            if sp and sp.is_available and sp.tool == "spliceai" and sp.max_delta is not None:
                if sp.max_delta >= 0.20:
                    return CriteriaResult.met(
                        ACMGCriterion.PP3, CriterionStrength.MODERATE,
                        f"SpliceAI max_delta={sp.max_delta:.3f} (Moderate) — missense with predicted splice impact",
                    )
            # Protein-level missense predictor: the user picks exactly ONE in
            # cfg.insilico_tool to avoid combining tools that share training
            # data (which would inflate evidence). ESM1b is preferred when
            # licence-compatible because it is fully open-source.
            if self._cfg.insilico_tool == InSilicoTool.ESM1B:
                es = annotation.esm1b
                if es and es.llr is not None:
                    strength = _esm1b_pp3(es.llr)
                    if strength:
                        return CriteriaResult.met(
                            ACMGCriterion.PP3, strength,
                            f"ESM1b LLR={es.llr:.3f} ({strength.value})",
                        )
                    # We distinguish "score present but not pathogenic" from
                    # "no score at all" — the former is informative for BP4
                    # and goes into the evidence trail.
                    return CriteriaResult.not_met(
                        ACMGCriterion.PP3,
                        f"ESM1b LLR={es.llr:.3f} (indeterminate or benign)",
                    )
                return CriteriaResult.not_met(ACMGCriterion.PP3, "No in-silico score available")

            am = annotation.alphamissense
            if am and am.score is not None:
                strength = _alphamissense_pp3(am.score)
                if strength:
                    return CriteriaResult.met(
                        ACMGCriterion.PP3, strength,
                        f"AlphaMissense={am.score:.3f} ({strength.value})",
                    )
                return CriteriaResult.not_met(
                    ACMGCriterion.PP3,
                    f"AlphaMissense={am.score:.3f} (indeterminate or benign)",
                )
            return CriteriaResult.not_met(ACMGCriterion.PP3, "No in-silico score available")

        # Splice-impacting non-missense classes: synonymous, intronic, and
        # the soft "splice_region" zone are evaluated by the splice
        # predictor alone — there is no protein change to score otherwise.
        if pc.consequence in (
            ConsequenceType.SPLICE_REGION,
            ConsequenceType.INTRON,
            ConsequenceType.SYNONYMOUS,
        ):
            sp = annotation.splice
            if sp and sp.is_available:
                if sp.tool == "spliceai" and sp.max_delta is not None:
                    strength = _spliceai_pp3(sp.max_delta)
                    score_str = f"SpliceAI max_delta={sp.max_delta:.3f}"
                elif sp.tool == "squirls" and sp.raw_score is not None:
                    strength = _squirls_pp3(sp.raw_score)
                    # Suffix the score string so reviewers see the calibration
                    # caveat inline alongside the trigger evidence.
                    score_str = f"SQUIRLS={sp.raw_score:.3f}"
                elif sp.tool == "mmsplice" and sp.raw_score is not None:
                    strength = _mmsplice_pp3(sp.raw_score)
                    score_str = f"MMSplice delta_logit_psi={sp.raw_score:.3f}"
                else:
                    return CriteriaResult.not_met(ACMGCriterion.PP3, "Splice score unavailable")

                if strength:
                    return CriteriaResult.met(
                        ACMGCriterion.PP3, strength,
                        f"{score_str} ({strength.value})",
                    )
            return CriteriaResult.not_met(ACMGCriterion.PP3, "Splice score not pathogenic")

        # Everything else (UTR, intergenic, etc.) is out of scope for PP3.
        return CriteriaResult.not_met(ACMGCriterion.PP3, "Consequence not applicable for PP3")
