"""BP7 -- synonymous/deep-intronic variant with no predicted splice impact (Walker 2023 expansion)."""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# Walker 2023 splice-impact safe-distance thresholds. Variants closer to
# the splice site are inside the canonical donor/acceptor consensus and
# cannot satisfy BP7 even with a low splice-tool score. The asymmetric
# cutoffs (+7 vs -21) reflect biology: branch points and polypyrimidine
# tracts extend much further upstream of the acceptor than downstream of
# the donor.
_DEEP_INTRONIC_DONOR_MIN = 7       # >= +7  (downstream of donor)
_DEEP_INTRONIC_ACCEPTOR_MAX = -21  # <= -21 (upstream of acceptor)
# Canonical splice dinucleotides are +/-1,2; the "noncanonical" VCEP range admits
# any intronic position beyond them, i.e. |distance| >= 3.
_NONCANONICAL_MIN_ABS = 3


def _intronic_eligible(pc, mode: str = "") -> bool:
    """Whether an intronic variant is far enough from the splice site for BP7.

    Default (Walker 2023): deep-intronic only (donor >= +7, acceptor <= -21).
    ``"noncanonical"`` (RASopathy / PIK3 VCEPs): any intronic position beyond the
    canonical +/-1,2 sites (|distance| >= 3); the benign-splice gate downstream
    still guards against a predicted cryptic site. Returns False when distance is
    unknown — we refuse to assume safety without evidence."""
    dist = pc.intron_distance_from_splice
    if dist is None:
        return False
    if mode == "noncanonical":
        return abs(dist) >= _NONCANONICAL_MIN_ABS
    return dist >= _DEEP_INTRONIC_DONOR_MIN or dist <= _DEEP_INTRONIC_ACCEPTOR_MAX


def _splice_benign(annotation: AnnotationData) -> bool:
    """Splice predictor agrees the variant has no impact.

    Walker 2023 calibrates SpliceAI ≤ 0.10 as "no impact" (Supporting BP4).
    SQUIRLS lacks an equivalent calibration so we use < 0.20 as a practical
    approximation matching the lower BP4 bound — see README caveat."""
    sp = annotation.splice
    if sp is None or not sp.is_available:
        return False
    if sp.tool == "spliceai" and sp.max_delta is not None:
        return sp.max_delta <= 0.10
    # OpenSpliceAI shares SpliceAI's 0–1 delta scale, so the same <= 0.10
    # "no impact" cutoff applies.
    if sp.tool == "openspliceai" and sp.max_delta is not None:
        return sp.max_delta <= 0.10
    if sp.tool == "squirls" and sp.raw_score is not None:
        return sp.raw_score < 0.20
    # MMSplice DISABLED — retained, commented out, for later:
    # if sp.tool == "mmsplice" and sp.raw_score is not None:
    #     # |delta_logit_psi| < 0.5 → minimal predicted splice effect (matches BP4).
    #     return abs(sp.raw_score) < 0.5
    return False


class BP7Evaluator(CriterionEvaluator):
    """
    Walker 2023 expanded BP7 applies to:
    1. Synonymous variants with no predicted splice impact
    2. Intronic variants >= +7 or <= -21 bp from canonical splice site
    """

    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        from acmg_classifier.local_db.conservation import PhyloPReader
        from acmg_classifier.criteria.bp_genes import BPApplicability
        self._phylop = PhyloPReader(cfg.phylop_bigwig)
        # Per-gene phyloP "highly conserved" cutoff from the VCEP specs; falls
        # back to the global default when a gene has no VCEP cutoff.
        self._spec = BPApplicability(cfg.disease_prevalence_tsv)

    def _conservation_block(self, variant: VariantRecord, gene: str | None) -> str | None:
        """Return a not-met reason when the position is highly conserved (so BP7
        must not fire), or None when the gate passes or is unavailable.

        ACMG/Walker BP7 requires the nucleotide to be NOT highly conserved. The
        cutoff is the gene's VCEP phyloP threshold when specified (e.g. CDH1-style
        neurodev panels 0.1, VHL 0.2, RPGR 0, GP genes 1.5), else the global
        default ``bp7_phylop_max`` (2.0, phyloP100way). The gate is applied only
        when phyloP is available; otherwise it is skipped (graceful degradation).

        A VCEP that declared conservation NON-informative (TP53, GALT, the SCID
        T-/B-cell genes) skips the gate entirely — phyloP is not consulted."""
        if self._spec.bp7_conservation_na(gene):
            return None
        if not self._phylop.is_available():
            return None
        score = self._phylop.value(variant.chrom, variant.pos)
        if score is None:
            return None
        cutoff = self._spec.bp7_phylop(gene)
        if cutoff is None:
            cutoff = self._cfg.bp7_phylop_max
        if score >= cutoff:
            return f"highly conserved (phyloP={score:.2f} >= {cutoff})"
        return None

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        pc = annotation.primary_consequence
        if pc is None:
            return CriteriaResult.not_met(ACMGCriterion.BP7, "No consequence")

        # Synonymous branch: requires a splice-tool agreement because a
        # silent codon change can still create/destroy a splice site. We do
        # NOT auto-fire on "synonymous alone" — that was the pre-Walker
        # behaviour and is now known to mis-classify exonic splice variants.
        if pc.consequence == ConsequenceType.SYNONYMOUS:
            if _splice_benign(annotation):
                blocked = self._conservation_block(variant, pc.gene_symbol)
                if blocked:
                    return CriteriaResult.not_met(
                        ACMGCriterion.BP7, f"Synonymous, splice benign, but {blocked}",
                    )
                sp = annotation.splice
                return CriteriaResult.met(
                    ACMGCriterion.BP7,
                    evidence=f"Synonymous + {sp.tool} score benign + not highly conserved",
                )
            return CriteriaResult.not_met(
                ACMGCriterion.BP7,
                "Synonymous but splice impact not ruled out",
            )

        # Intronic branch: BP7 requires BOTH that the variant is outside the
        # splice consensus (distance) AND that splicing-prediction algorithms
        # predict no impact (Walker 2023 / ClinGen SVI). Distance alone is NOT
        # sufficient — a deep-intronic variant can still create a cryptic splice
        # site — so we withhold BP7 when no splice predictor confirms no impact.
        #
        # The eligible distance range is per-gene: the Walker default admits only
        # DEEP-intronic variants (+7/-21, consequence INTRON); the RASopathy /
        # PIK3 VCEPs admit any intronic position except the canonical +/-1,2
        # sites, which also pulls in the +3..+8 / -3..-8 SPLICE_REGION variants.
        mode = self._spec.bp7_intronic_mode(pc.gene_symbol)
        intronic_types = (
            (ConsequenceType.INTRON, ConsequenceType.SPLICE_REGION)
            if mode == "noncanonical"
            else (ConsequenceType.INTRON,)
        )
        if pc.consequence in intronic_types:
            if not _intronic_eligible(pc, mode):
                return CriteriaResult.not_met(
                    ACMGCriterion.BP7,
                    f"Intronic but within splice consensus (dist="
                    f"{pc.intron_distance_from_splice})",
                )
            if not _splice_benign(annotation):
                return CriteriaResult.not_met(
                    ACMGCriterion.BP7,
                    "Intronic but splice impact not ruled out "
                    "(no splice prediction of no impact)",
                )
            blocked = self._conservation_block(variant, pc.gene_symbol)
            if blocked:
                return CriteriaResult.not_met(
                    ACMGCriterion.BP7, f"Intronic, splice benign, but {blocked}",
                )
            range_desc = (
                "intronic outside canonical +/-1,2"
                if mode == "noncanonical"
                else "deep intronic (+7/-21)"
            )
            return CriteriaResult.met(
                ACMGCriterion.BP7,
                evidence=(
                    f"{range_desc} (dist={pc.intron_distance_from_splice}) "
                    "and splice score benign + not highly conserved (Walker 2023)"
                ),
            )

        return CriteriaResult.not_met(
            ACMGCriterion.BP7,
            "Not a synonymous or deep-intronic variant",
        )
