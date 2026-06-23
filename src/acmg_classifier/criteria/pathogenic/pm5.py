"""PM5 -- different missense at same codon as established pathogenic variant."""
from __future__ import annotations

import re

from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.blosum62 import blosum62_score
from acmg_classifier.criteria.grantham import AA3_TO_AA1, grantham_distance
from acmg_classifier.criteria.pm5_genes import PM5Spec
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
            query_chrom=variant.chrom,
            query_pos=variant.pos,
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

        # Chemical-severity gate (PIK3CD, PIK3R1, VHL via Grantham; PTEN via
        # BLOSUM62): keep only comparators the candidate is chemically
        # as-severe-or-more than at the wild-type residue.
        gate = self._spec.gate(pc.gene_symbol)
        if gate:
            matrix, op = gate
            qualifying = self._gate_filter(pc, hits, matrix, op)
            engine = "Grantham" if matrix == "grantham" else "BLOSUM62"
            if qualifying is None:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM5, f"{engine} score unavailable for candidate"
                )
            if not qualifying:
                sym = {"ge": ">=", "gt": ">", "le": "<=", "lt": "<"}[op]
                return CriteriaResult.not_met(
                    ACMGCriterion.PM5,
                    f"{engine} gate failed: candidate not {sym} comparator",
                )
            tag = f"{engine}-gated ({ {'ge': '>=', 'gt': '>', 'le': '<=', 'lt': '<'}[op] })"
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

        cap = self._spec.max_strength(pc.gene_symbol) or CriterionStrength.MODERATE

        # Count-threshold VCEPs (ACVRL1/ENG, HHT): PM5 applies only when >=N
        # DISTINCT same-codon LP/P missense comparators exist, then fires at the
        # VCEP ceiling (Strong). Fewer than N distinct changes → PM5 not met.
        min_count = self._spec.min_count(pc.gene_symbol)
        if min_count > 1:
            distinct = len({
                _comparator_change(h.hgvs_p) for h in qualifying
                if _comparator_change(h.hgvs_p)
            })
            if distinct < min_count:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM5,
                    f"only {distinct} distinct same-codon LP/P comparator(s) "
                    f"(< {min_count} required for {pc.gene_symbol})",
                )
            ids = ", ".join(h.variation_id or "" for h in qualifying[:3])
            return CriteriaResult.met(
                ACMGCriterion.PM5, strength=cap,
                evidence=(f"ClinVar {tag}: >={min_count} distinct same-codon LP/P "
                          f"({distinct}): {ids}"),
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

        strength = _min_strength(base, cap)
        ids = ", ".join(h.variation_id or "" for h in qualifying[:3])
        return CriteriaResult.met(
            ACMGCriterion.PM5, strength=strength, evidence=f"ClinVar {tag}: {ids}"
        )

    def _gate_filter(
        self, pc, hits: list[ClinVarRecord], matrix: str, op: str
    ) -> list[ClinVarRecord] | None:
        """Comparators the candidate satisfies the chemical-severity gate
        against, using the Grantham or BLOSUM62 engine.

        Returns ``None`` when the candidate's own score cannot be computed (the
        mandated comparison is impossible → withhold PM5)."""
        cand = _candidate_score(pc.amino_acids, matrix)
        if cand is None:
            return None
        out: list[ClinVarRecord] = []
        for h in hits:
            comp = _comparator_score(h.hgvs_p, matrix)
            if comp is None:
                continue
            if _gate_pass(cand, comp, op):
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


def _score_pair(aa_a: str | None, aa_b: str | None, matrix: str) -> int | None:
    """Grantham distance or BLOSUM62 similarity between two residues."""
    if matrix == "blosum":
        return blosum62_score(aa_a, aa_b)
    return grantham_distance(aa_a, aa_b)


def _gate_pass(cand: int, comp: int, op: str) -> bool:
    """True if the candidate score satisfies the gate against the comparator.
    Grantham distance is severity-increasing (``ge``/``gt``); BLOSUM62
    similarity is severity-decreasing (``le``/``lt``)."""
    if op == "ge":
        return cand >= comp
    if op == "gt":
        return cand > comp
    if op == "le":
        return cand <= comp
    return cand < comp  # "lt"


def _candidate_score(amino_acids: str | None, matrix: str) -> int | None:
    """Chemical-severity score of the candidate from its "REF/ALT" 1-letter pair."""
    if not amino_acids:
        return None
    parts = amino_acids.split("/")
    if len(parts) != 2:
        return None
    return _score_pair(parts[0].strip(), parts[1].strip(), matrix)


def _comparator_score(hgvs_p: str | None, matrix: str) -> int | None:
    """Chemical-severity score of a ClinVar comparator from its 3-letter change."""
    if not hgvs_p:
        return None
    m = _COMP_P.search(hgvs_p)
    if not m:
        return None
    return _score_pair(AA3_TO_AA1.get(m.group(1)), AA3_TO_AA1.get(m.group(2)), matrix)
