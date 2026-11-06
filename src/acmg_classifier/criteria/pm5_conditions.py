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
* **Paralogue conservation** — KCNQ1 (LQTS VCEP): the residue must be identical
  in all five human KCNQ paralogues (``PARALOG_CONSERVATION``; table built by
  ``scripts/build_kcnq1_paralog_conservation.py``).
* Count-based strengths (F8 / F9 "1 pathogenic variant or 2 likely pathogenic
  variants", MYOC, ACADVL, LGMD, VWF, RASopathy, KCNQ1 P/LP): see
  :mod:`acmg_classifier.criteria.comparator_tiers`.
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
    # LQTS VCEP (KCNQ1): "SpliceAI must be used to examine both variants for
    # similar predicted effect or lack of effect on splicing (SpliceAI delta <0.2)"
    "KCNQ1": (0.20, "both", "max"),
}

REQUIRES_PP3: frozenset[str] = frozenset({"MLH1", "MSH2", "MSH6", "PMS2", "RPGR", "RS1"})

REVEL_MIN: dict[str, float] = {"DYSF": 0.7}

GRANTHAM_SUPPORTING: frozenset[str] = frozenset({"GCK", "HNF1A", "HNF4A"})

# F8 / F9 "1 pathogenic variant or 2 likely pathogenic variants" -> PM5 Moderate:
# criteria/comparator_tiers.py (PM5_RULES), with the other count-based VCEPs.

PARALOG_SUPPORTING: frozenset[str] = frozenset({"SCN1A", "SCN2A", "SCN3A", "SCN8A"})

# HHT VCEP (ACVRL1 GN135 / ENG GN136): "Beware of changes that impact splicing
# rather than at the amino acid/protein level" -- a caveat, not an exclusion:
# PM5 is kept and the evidence carries a note when the candidate (or, with the
# precomputed SpliceAI files, a comparator) reaches the gene's PP3 splice cutoff.
SPLICE_CAVEAT: frozenset[str] = frozenset({"ACVRL1", "ENG"})


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


# Paralogue-conservation gate: gene -> table of the gene's residues and the
# residue each paralogue has at the aligned position (``conserved`` = 1 when all
# match). LQTS VCEP (KCNQ1): "Residue must be highly conserved across all 5 human
# KCNQ paralogues ... Poor conservation (ineligibility for PM5) is defined whenever
# 1 or more KCNQ's show a different amino acid at the position."
PARALOG_CONSERVATION: dict[str, str] = {"KCNQ1": "kcnq1_paralog_conservation.tsv"}


def _load_conservation(name: str) -> dict[int, tuple[str, dict[str, str], bool]]:
    import csv
    from pathlib import Path
    path = Path(__file__).resolve().parents[3] / "resources" / "shared" / name
    out: dict[int, tuple[str, dict[str, str], bool]] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8") as fh:
        rows = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        others = [c for c in (rows.fieldnames or []) if c.startswith("kcnq") and c != "kcnq1"]
        for r in rows:
            out[int(r["position"])] = (r["kcnq1"], {c.upper(): r[c] for c in others},
                                       r["conserved"] == "1")
    return out


_CONS_CACHE: dict[str, dict] = {}


def paralog_conserved(gene: str, position: int | None, ref_aa: str | None
                      ) -> tuple[bool | None, str]:
    """(conserved?, detail) for the gene's paralogue-conservation gate; None when
    the gene has no gate or the residue cannot be checked."""
    name = PARALOG_CONSERVATION.get(gene)
    if name is None or position is None:
        return None, ""
    if name not in _CONS_CACHE:
        _CONS_CACHE[name] = _load_conservation(name)
    row = _CONS_CACHE[name].get(position)
    if row is None:
        return None, f"{gene} paralogue table unavailable for residue {position}"
    aa, others, cons = row
    if ref_aa and ref_aa != aa:
        return None, f"{gene} residue {position} is {aa} in the table, variant reports {ref_aa}"
    detail = ", ".join(f"{g} {a}" for g, a in others.items())
    return cons, f"{gene} {aa}{position}: {detail}"
