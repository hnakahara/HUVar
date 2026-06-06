"""PM5 -- different missense at same codon as established pathogenic variant."""
from __future__ import annotations

import re

from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.grantham import AA3_TO_AA1, grantham_distance
from acmg_classifier.criteria.pm5_genes import GT, PM5Spec
from acmg_classifier.models.annotation import AnnotationData, ClinVarRecord
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

# A comparator protein change as ClinVar stores it (3-letter, e.g.
# "NP_000298.6:p.Arg156His"). Captures wild-type and variant residues.
_COMP_P = re.compile(r"p\.([A-Z][a-z]{2})\d+([A-Z][a-z]{2})")

# ClinVar significances that count as "pathogenic" (PM5 default Moderate). When
# every same-codon comparator is only "Likely pathogenic" the strength drops to
# Supporting, per the VCEP specifications.
_PATHOGENIC = {"Pathogenic", "Pathogenic/Likely pathogenic"}


class PM5Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Per-gene PM5 spec: Grantham gate, PM1/PS1 exclusions, strength ceiling.
        self._spec = PM5Spec(cfg.disease_prevalence_tsv)

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

        from acmg_classifier.local_db.clinvar_sqlite import (
            has_benign_at_codon,
            query_same_codon_different_aa,
        )

        # Same-codon P/LP missense comparators with a DIFFERENT amino-acid
        # change (an exact match would be PS1, not PM5). min_stars defaults to 2
        # (multiple submitters / expert) so single-submitter assertions do not
        # anchor PM5 — the dominant source of PM5 over-assignment.
        hits = query_same_codon_different_aa(
            self._cfg.clinvar_sqlite,
            pc.gene_symbol,
            pc.protein_position,
            pc.hgvs_p,
            min_stars=self._cfg.pm5_min_stars,
        )
        if not hits:
            return CriteriaResult.not_met(
                ACMGCriterion.PM5,
                f"No ClinVar >={self._cfg.pm5_min_stars} star same-codon different-AA hit",
            )

        # Do not apply at a codon where any benign variant is known (a VCEP
        # caveat for the Grantham-gated genes, applied to every gene here as a
        # conservative, false-positive-reducing constraint).
        if has_benign_at_codon(
            self._cfg.clinvar_sqlite, pc.gene_symbol, pc.protein_position,
            min_stars=self._cfg.pm5_min_stars,
        ):
            return CriteriaResult.not_met(
                ACMGCriterion.PM5, "Benign variant known at codon"
            )

        # Grantham-distance gate (PIK3CD, PIK3R1, VHL, …): keep only comparators
        # the candidate is chemically as/more different from the wild type than.
        op = self._spec.operator(pc.gene_symbol)
        if op:
            qualifying = self._grantham_filter(pc, hits, op)
            if qualifying is None:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM5, "Grantham distance unavailable for candidate"
                )
            if not qualifying:
                sym = ">" if op == GT else ">="
                return CriteriaResult.not_met(
                    ACMGCriterion.PM5, f"Grantham gate failed: candidate not {sym} comparator"
                )
            tag = f"Grantham-gated ({'>' if op == GT else '>='})"
        else:
            qualifying = hits
            tag = "same codon, diff AA"

        # Comparator-significance policy: VCEPs that offer no Supporting PM5
        # strength require the comparator to reach Pathogenic (PTEN, VHL, KCNQ1,
        # RASopathy genes, …) — a Likely-pathogenic-only comparator must not
        # trigger PM5. Genes without this restriction keep the LP->Supporting path.
        if self._spec.requires_pathogenic(pc.gene_symbol):
            qualifying = [h for h in qualifying if h.clinical_significance in _PATHOGENIC]
            if not qualifying:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM5,
                    "No Pathogenic same-codon comparator (LP not accepted for this gene)",
                )

        # Strength (ClinGen SVI):
        #   >=2 DIFFERENT pathogenic missense at the codon -> Strong
        #   1 pathogenic comparator                         -> Moderate
        #   only likely-pathogenic comparators              -> Supporting
        # then clamped to the gene's VCEP ceiling. PM5_Strong is granted only
        # when the VCEP explicitly allows it (cap == Strong); genes without a
        # VCEP keep the historical Moderate ceiling.
        pathogenic_hits = [h for h in qualifying if h.clinical_significance in _PATHOGENIC]
        n_distinct_path = len({
            _comparator_change(h.hgvs_p) for h in pathogenic_hits if _comparator_change(h.hgvs_p)
        })
        if n_distinct_path >= 2:
            base = CriterionStrength.STRONG
        elif pathogenic_hits:
            base = CriterionStrength.MODERATE
        else:
            base = CriterionStrength.SUPPORTING

        cap = self._spec.max_strength(pc.gene_symbol) or CriterionStrength.MODERATE
        strength = _min_strength(base, cap)
        ids = ", ".join(h.variation_id or "" for h in qualifying[:3])
        return CriteriaResult.met(
            ACMGCriterion.PM5, strength=strength, evidence=f"ClinVar {tag}: {ids}"
        )

    def _grantham_filter(
        self, pc, hits: list[ClinVarRecord], op: str
    ) -> list[ClinVarRecord] | None:
        """Comparators the candidate satisfies the Grantham gate against.

        Returns ``None`` when the candidate's own Grantham distance cannot be
        computed (the mandated comparison is impossible → withhold PM5)."""
        cand = _candidate_distance(pc.amino_acids)
        if cand is None:
            return None
        out: list[ClinVarRecord] = []
        for h in hits:
            comp = _comparator_distance(h.hgvs_p)
            if comp is None:
                continue
            if (cand > comp) if op == GT else (cand >= comp):
                out.append(h)
        return out


_STRENGTH_ORDER = {
    CriterionStrength.SUPPORTING: 1,
    CriterionStrength.MODERATE: 2,
    CriterionStrength.STRONG: 3,
}


def _min_strength(a: CriterionStrength, b: CriterionStrength) -> CriterionStrength:
    """The weaker of two strengths (clamp a base strength to a ceiling)."""
    return a if _STRENGTH_ORDER[a] <= _STRENGTH_ORDER[b] else b


def _comparator_change(hgvs_p: str | None) -> str | None:
    """The comparator's amino-acid change (e.g. "ArgHis") from its protein
    HGVS, used to count DISTINCT pathogenic missense at the codon for the
    PM5_Strong (>=2 different missense) tier."""
    if not hgvs_p:
        return None
    m = _COMP_P.search(hgvs_p)
    return (m.group(1) + m.group(2)) if m else None


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
