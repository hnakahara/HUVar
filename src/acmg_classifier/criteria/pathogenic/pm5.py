"""PM5 -- different missense at same codon as established pathogenic variant."""
from __future__ import annotations

import re

from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.grantham import AA3_TO_AA1, grantham_distance
from acmg_classifier.criteria.pm5_genes import GT, PM5Grantham
from acmg_classifier.models.annotation import AnnotationData, ClinVarRecord
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# A comparator protein change as ClinVar stores it (3-letter, e.g.
# "NP_000298.6:p.Arg156His"). Captures wild-type and variant residues.
_COMP_P = re.compile(r"p\.([A-Z][a-z]{2})\d+([A-Z][a-z]{2})")

# ClinVar significances that count as "pathogenic" (PM5 default Moderate). A
# "Likely pathogenic"-only comparator drops PM5 to Supporting per the VCEP texts.
_PATHOGENIC = {"Pathogenic", "Pathogenic/Likely pathogenic"}


class PM5Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Per-gene Grantham-distance gate (PIK3CD, PIK3R1, VHL, …). Genes absent
        # from the table keep the plain same-codon-different-AA behaviour.
        self._grantham = PM5Grantham(cfg.disease_prevalence_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # PM5 is the "different AA, same codon" companion to PS1's "same AA"
        # rule, so only missense is eligible — splice variants don't have a
        # comparable codon-level interpretation, and synonymous changes by
        # definition do not change the amino acid.
        pc = annotation.primary_consequence
        if pc is None or pc.consequence != ConsequenceType.MISSENSE:
            return CriteriaResult.not_met(ACMGCriterion.PM5, "Not a missense variant")

        # Look for ClinVar P/LP variants at the same protein position with a
        # DIFFERENT amino-acid substitution (pc.hgvs_p is passed so the query
        # can exclude exact-match variants — those are PS1, not PM5).
        # min_stars=1 enforces the ACMG requirement for a reviewed assertion.
        from acmg_classifier.local_db.clinvar_sqlite import query_same_codon_different_aa
        hits = query_same_codon_different_aa(
            self._cfg.clinvar_sqlite,
            pc.gene_symbol,
            pc.protein_position,
            pc.hgvs_p,
            min_stars=1,
        )
        if not hits:
            return CriteriaResult.not_met(
                ACMGCriterion.PM5, "No ClinVar >=1 star same-codon different-AA hit"
            )

        op = self._grantham.operator(pc.gene_symbol)
        if not op:
            # No VCEP Grantham gate for this gene — plain PM5 at default Moderate.
            evidence = "ClinVar same codon, diff AA: " + ", ".join(
                h.variation_id or "" for h in hits[:3]
            )
            return CriteriaResult.met(ACMGCriterion.PM5, evidence=evidence)

        return self._evaluate_grantham(pc, hits, op)

    def _evaluate_grantham(
        self, pc, hits: list[ClinVarRecord], op: str
    ) -> CriteriaResult:
        """Apply the VCEP Grantham-distance gate (PIK3CD-style genes).

        The candidate must be chemically as different (``ge``) or more different
        (``gt``) from the wild-type residue than the comparator, no benign
        variant may be known at the codon, and the strength follows the
        comparator's classification (Pathogenic → Moderate, Likely pathogenic →
        Supporting)."""
        from acmg_classifier.local_db.clinvar_sqlite import has_benign_at_codon

        # VCEP caveat: do not apply at a codon where any benign variant is known.
        if has_benign_at_codon(
            self._cfg.clinvar_sqlite, pc.gene_symbol, pc.protein_position
        ):
            return CriteriaResult.not_met(
                ACMGCriterion.PM5, "Grantham-gated gene: benign variant known at codon"
            )

        cand = _candidate_distance(pc.amino_acids)
        if cand is None:
            # Cannot run the mandated comparison → withhold PM5 (precision-first).
            return CriteriaResult.not_met(
                ACMGCriterion.PM5, "Grantham distance unavailable for candidate"
            )

        qualifying: list[tuple[ClinVarRecord, int]] = []
        for h in hits:
            comp = _comparator_distance(h.hgvs_p)
            if comp is None:
                continue
            if (cand > comp) if op == GT else (cand >= comp):
                qualifying.append((h, comp))

        if not qualifying:
            sym = ">" if op == GT else ">="
            return CriteriaResult.not_met(
                ACMGCriterion.PM5,
                f"Grantham gate failed: candidate {cand} not {sym} any comparator",
            )

        # Strength: Moderate if a Pathogenic comparator qualifies, else (only
        # Likely pathogenic comparators) Supporting, per the VCEP specifications.
        pathogenic = [
            (h, d) for h, d in qualifying if h.clinical_significance in _PATHOGENIC
        ]
        chosen = pathogenic or qualifying
        strength = (
            CriterionStrength.MODERATE if pathogenic else CriterionStrength.SUPPORTING
        )
        ids = ", ".join(h.variation_id or "" for h, _ in chosen[:3])
        evidence = (
            f"Grantham-gated PM5 (cand {cand} {'>' if op == GT else '>='} "
            f"comparator): {ids}"
        )
        return CriteriaResult.met(ACMGCriterion.PM5, strength=strength, evidence=evidence)


def _candidate_distance(amino_acids: str | None) -> int | None:
    """Grantham distance of the candidate from its "REF/ALT" 1-letter pair."""
    if not amino_acids:
        return None
    parts = amino_acids.split("/")
    if len(parts) != 2:
        return None
    return grantham_distance(parts[0].strip(), parts[1].strip())


def _comparator_distance(hgvs_p: str | None) -> int | None:
    """Grantham distance of a ClinVar comparator from its 3-letter protein change."""
    if not hgvs_p:
        return None
    m = _COMP_P.search(hgvs_p)
    if not m:
        return None
    return grantham_distance(AA3_TO_AA1.get(m.group(1)), AA3_TO_AA1.get(m.group(2)))
