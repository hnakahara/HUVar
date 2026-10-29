"""AMP/ASCO/CAP 2017 somatic variant tiers (Li et al., J Mol Diagn 2017;19:4-23).

Evidence (CIViC predictive / diagnostic / prognostic items and assertions; ClinVar
Somatic Clinical Impact) is graded against the patient's tumour type
(OncoTree code -> Disease Ontology, see tumor.py):

  CIViC level  same tumour type (incl. a more general disease such as "cancer")   other tumour type
  A            Tier I, level A                                                     Tier II, level C
  B            Tier I, level B                                                     Tier II, level C
  C            Tier II, level C                                                    Tier II, level D
  D, E         Tier II, level D                                                    Tier II, level D

Assertions carry their own AMP category, used as such for the same tumour type and
capped at Tier II level C otherwise. Evidence recorded for a *subtype* of the tumour
type counts as the same tumour type and is flagged. Only "Supports" records count.

Without Tier I/II evidence: Tier IV when the variant is benign / likely benign or
common in gnomAD (SBVS1 / SBS1), otherwise Tier III. Without a tumour type no tier
is assigned; the strongest evidence in any tumour type is listed for reference.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from acmg_classifier.somatic.tumor import DiseaseOntology, TumorType, relation

_RANK = {("I", "A"): 6, ("I", "B"): 5, ("I", ""): 5, ("II", "C"): 4, ("II", "D"): 3,
         ("II", ""): 3}
_CATEGORY = {"Predictive": "therapeutic", "Diagnostic": "diagnostic", "Prognostic": "prognostic"}


@dataclass
class TierCall:
    tier: str                         # I / II
    level: str                        # A-D ("" when the source gives none)
    category: str                     # therapeutic / diagnostic / prognostic
    text: str
    source: str                       # CIViC / ClinVar
    evidence_level: str = ""          # CIViC evidence level A-E; "assertion"; "" (ClinVar)
    rating: int = 0

    @property
    def rank(self) -> int:
        return _RANK.get((self.tier, self.level), 0)

    @property
    def sort_key(self) -> tuple:
        """Strongest first: tier/level, then assertions, CIViC level, rating, CIViC
        before ClinVar."""
        lvl = {"assertion": 0, "A": 1, "B": 2, "C": 3, "D": 4, "E": 5}.get(self.evidence_level, 6)
        return (-self.rank, lvl, -self.rating, self.source != "CIViC")

    @property
    def label(self) -> str:
        return f"Tier {self.tier}" + (f", level {self.level}" if self.level else "")


@dataclass
class AmpResult:
    tier: str = ""                    # "Tier I" ... "Tier IV"; "" = no tumour type
    level: str = ""
    tumor: str = ""
    calls: list[TierCall] = field(default_factory=list)   # supporting calls, best first
    candidates: list[str] = field(default_factory=list)   # no tumour type: reference list
    note: str = ""

    def evidence_text(self, n: int = 6) -> str:
        shown = list(self.calls[:n])
        # Always show the strongest ClinVar call when CIViC fills the list.
        cv = next((c for c in self.calls if c.source == "ClinVar"), None)
        if cv is not None and cv not in shown:
            shown.append(cv)
        return " | ".join(f"{c.label} ({c.category}, {c.source}): {c.text}" for c in shown)


def _grade(level: str, rel: str) -> tuple[str, str]:
    level = (level or "").upper()
    same = rel in ("same", "subtype")
    if level in ("A", "B"):
        return ("I", level) if same else ("II", "C")
    if level == "C":
        return ("II", "C") if same else ("II", "D")
    return ("II", "D")


def _parse_amp(cat: str) -> tuple[str, str]:
    """'Tier I - Level A' -> ('I', 'A')."""
    c = (cat or "").replace("Tier", "").strip()
    tier, _, lvl = c.partition("-")
    return tier.strip(), lvl.replace("Level", "").strip()


def civic_calls(records: list[tuple[dict, str]], tt: Optional[TumorType],
                do: DiseaseOntology) -> list[TierCall]:
    calls = []
    for r, how in records:
        if r.get("direction") != "Supports":
            continue
        rel = relation(tt, r.get("doid"), do) if tt is not None else "unknown"
        ref = ("AID" if r["record"] == "assertion" else "EID") + r["record_id"]
        disease = r.get("disease") or "?"
        drugs = f" {r['therapies']}" if r.get("therapies") else ""
        sub = " [evidence in a subtype]" if rel == "subtype" else ""
        if r["record"] == "assertion":
            tier, lvl = _parse_amp(r.get("amp_category", ""))
            if tier not in ("I", "II"):
                continue
            if rel not in ("same", "subtype") and tier == "I":
                tier, lvl = "II", "C"
            extra = "; NCCN " + r["nccn_guideline"] if r.get("nccn_guideline") else ""
            extra += "; regulatory approval" if r.get("regulatory_approval") == "true" else ""
            text = (f"CIViC {ref} {r['significance']}{drugs} in {disease}{extra} "
                    f"[{how}]{sub}")
        else:
            tier, lvl = _grade(r.get("level", ""), rel)
            text = (f"CIViC {ref} level {r.get('level') or '?'} {r['significance']}{drugs} in "
                    f"{disease} (PMID {r.get('citation_id') or '-'}, rating "
                    f"{r.get('rating') or '-'}) [{how}]{sub}")
        try:
            rating = int(float(r.get("rating") or 0))
        except ValueError:
            rating = 0
        calls.append(TierCall(tier, lvl, _CATEGORY.get(r["evidence_type"], r["evidence_type"]),
                              text, "CIViC",
                              "assertion" if r["record"] == "assertion" else
                              (r.get("level") or "").upper(), rating))
    return calls


def clinvar_calls(records: list[dict], tt: Optional[TumorType], do: DiseaseOntology) -> list[TierCall]:
    """ClinVar Somatic Clinical Impact of the variant itself (one dict per RCV:
    sci, sci_review_status, sci_stars, sci_type, conditions)."""
    calls = []
    for r in records:
        sci = r.get("sci") or ""
        if not sci.lower().startswith(("tier i ", "tier i-", "tier ii")):
            continue
        tier = "II" if sci.lower().startswith("tier ii") else "I"
        rel = "unknown"
        if tt is not None:
            rels = [relation(tt, do.lookup_name(c), do) for c in (r.get("conditions") or "").split("; ")]
            rel = "same" if "same" in rels else "subtype" if "subtype" in rels else (
                "other" if "other" in rels else "unknown")
        if tier == "I" and rel not in ("same", "subtype"):
            tier, lvl = "II", "C"
        else:
            lvl = "" if tier == "I" else ("C" if rel in ("same", "subtype") else "D")
        kind = (r.get("sci_type") or "").split(";")[0].strip().lower()
        category = ("therapeutic" if "therap" in kind else "diagnostic" if "diagn" in kind
                    else "prognostic" if "progn" in kind else "clinical impact")
        calls.append(TierCall(
            tier, lvl, category,
            f"ClinVar somatic clinical impact '{sci}' ({r.get('sci_type') or '-'}) in "
            f"{r.get('conditions') or '?'} (ClinVar record {r.get('variation_id')}, {r.get('sci_stars', 0)} "
            f"stars: {r.get('sci_review_status') or '-'})", "ClinVar"))
    return calls


def assign(calls: list[TierCall], tt: Optional[TumorType], do: DiseaseOntology,
           benign: bool, common: bool) -> AmpResult:
    calls = sorted(calls, key=lambda c: c.sort_key)
    if tt is None:
        cand = [f"({c.category}) {c.text}" for c in calls if c.rank >= _RANK[("II", "C")]][:8]
        return AmpResult(candidates=cand, note="no tumour type given (--tumor-type): tier not "
                                                "assigned; strongest evidence in any tumour type listed")
    res = AmpResult(tumor=tt.describe(do))
    best = [c for c in calls if c.rank > 0]
    if best:
        top = best[0]
        res.tier, res.level, res.calls = f"Tier {top.tier}", top.level, best
        if top.source == "ClinVar" and top.tier == "I":
            res.note = "Tier I from ClinVar somatic clinical impact (level not given)"
        return res
    if benign or common:
        res.tier = "Tier IV"
        res.note = ("benign / likely benign oncogenicity" if benign else
                    "observed at significant frequency in gnomAD (SBVS1/SBS1)")
    else:
        res.tier = "Tier III"
        res.note = "no Tier I/II evidence for this tumour type (CIViC / ClinVar)"
    return res
