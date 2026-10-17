"""Gene-specific PM5 conditions from ClinGen VCEP specifications (CSpec text).

Complements :mod:`acmg_classifier.criteria.pm5_genes` (Grantham / BLOSUM gate,
exclusions, strength ceiling, Pathogenic-only comparators, distinct-change
count) with the conditions that VCEPs attach to PM5 for missense variants:

* **Splicing must be ruled out** — PM5 is a protein-level argument, so most
  VCEPs require that the candidate (and often the comparator) is not predicted
  to act through splicing. ``SPLICE`` maps a gene to
  ``(cutoff, scope, mode)``: ``mode`` ``"max"`` = the SpliceAI max delta must be
  below ``cutoff`` (``"le"`` = at or below); ``scope`` ``"both"`` = applied to the
  candidate and the comparator(s), ``"candidate"`` / ``"comparator"`` = one side.
  ``"category"`` (GUCY2D) = both scores in the same band (<0.1, 0.1-0.2, >0.2);
  ``"relative"`` (RPE65) = the two scores within 10% of each other.
* **PP3 must be met** — InSiGHT MMR genes ("Only use PM5 if PP3 is supporting"),
  RPGR / RS1 ("must meet PP3"). Enforced after all criteria are evaluated
  (``REQUIRES_PP3``; see ``CriteriaRegistry``).
* **REVEL minimum** — DYSF: "the REVEL score for the missense variant under
  curation should be >0.7".
* **Grantham downgrade** — Monogenic Diabetes VCEP (GCK / HNF1A / HNF4A):
  "Supporting ... if the previously classified variant is pathogenic but has a
  greater Grantham distance" (instead of withholding PM5).
* **Two likely-pathogenic comparators = Moderate** — F8 / F9: "1 pathogenic
  variant or 2 likely pathogenic variants at the same residue".
* **Paralogue route** — SCN1A/2A/3A/8A: PM5_Supporting for ">1 non-identical
  amino-acid change in paralogous gene(s)" classified P/LP at the analogous
  residue.

Comparator splice scores are looked up in the precomputed SpliceAI files when
the run uses SpliceAI; with other splice tools only the candidate is checked
(the evidence says so).
"""
from __future__ import annotations

# gene -> (cutoff, scope, mode)
SPLICE: dict[str, tuple[float, str, str]] = {
    # "neither change should be predicted to affect splicing (SpliceAI < 0.2)"
    "CTLA4": (0.20, "both", "max"),
    "PIK3CD": (0.20, "both", "max"),
    "PIK3R1": (0.20, "both", "max"),
    # "must not affect splicing (SpliceAI <= 0.2)" (candidate)
    "RPGR": (0.20, "candidate", "le"),
    "RS1": (0.20, "candidate", "le"),
    # "exclude likely splice effects (SpliceAI < 0.5)" for both
    "DYSF": (0.50, "both", "max"),
    # "Splice AI score must be < 0.1 for both"
    "SLC6A8": (0.10, "both", "max"),
    # "in silico predictors do not suggest a splicing defect" (PP3 cutoff 0.2)
    "F8": (0.20, "candidate", "max"),
    "F9": (0.20, "candidate", "max"),
    # "splicing abnormalities (VCEP-specified prediction algorithms) excluded"
    "ABCD1": (0.50, "candidate", "max"),      # ABCD1 PP3 splice cutoff >= 0.5
    "GALT": (0.20, "candidate", "max"),
    "VHL": (0.20, "candidate", "max"),        # "splice metapredictors"
    "DICER1": (0.20, "candidate", "max"),     # "agreement ... no splicing effects"
    "GAMT": (0.20, "candidate", "max"),       # last-3-nt splicing check
    # InSiGHT MMR: comparator "Pathogenic on the protein level and not due to
    # aberrant splicing"; candidate not spliceogenic (PP3 splice cutoff 0.2)
    "MLH1": (0.20, "both", "max"),
    "MSH2": (0.20, "both", "max"),
    "MSH6": (0.20, "both", "max"),
    "PMS2": (0.20, "both", "max"),
    # "Do not apply if the comparison variant is suspected to cause a splicing
    # defect via SpliceAI" (ABCA4 PP3 splice cutoff 0.2)
    "ABCA4": (0.20, "comparator", "max"),
    # "SpliceAI scores for both variants should be in the same category"
    "GUCY2D": (0.0, "both", "category"),
    # "SpliceAI scores for both variants should be within 10% of each other"
    "RPE65": (0.10, "both", "relative"),
}

REQUIRES_PP3: frozenset[str] = frozenset({"MLH1", "MSH2", "MSH6", "PMS2", "RPGR", "RS1"})

REVEL_MIN: dict[str, float] = {"DYSF": 0.7}

GRANTHAM_SUPPORTING: frozenset[str] = frozenset({"GCK", "HNF1A", "HNF4A"})

TWO_LP_MODERATE: frozenset[str] = frozenset({"F8", "F9"})

PARALOG_SUPPORTING: frozenset[str] = frozenset({"SCN1A", "SCN2A", "SCN3A", "SCN8A"})


def _band(x: float) -> int:
    return 0 if x < 0.1 else (1 if x <= 0.2 else 2)


def splice_ok(gene: str, value: float | None, side: str,
              other: float | None = None) -> bool | None:
    """Whether a SpliceAI score passes the gene's PM5 splice condition.

    ``side`` is "candidate" or "comparator"; ``other`` is the candidate score
    when checking a comparator (for the paired modes). Returns None when the
    gene has no condition for this side or the score is unknown."""
    spec = SPLICE.get(gene)
    if spec is None or value is None:
        return None
    cutoff, scope, mode = spec
    if scope not in ("both", side):
        return None
    if mode == "max":
        return value < cutoff
    if mode == "le":
        return value <= cutoff
    if side == "candidate":
        return None  # paired modes are judged on the comparator
    if other is None:
        return None
    if mode == "category":
        return _band(value) == _band(other)
    if mode == "relative":
        return abs(value - other) <= cutoff
    return None
