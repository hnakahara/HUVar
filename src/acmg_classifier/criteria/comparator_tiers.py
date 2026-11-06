"""VCEP PS1 / PM5 strengths that depend on HOW MANY comparators are classified
Pathogenic or Likely pathogenic.

The ClinGen SVI default (any Pathogenic comparator -> the full strength, Likely
pathogenic only -> one step down; for PM5 >=2 different Pathogenic changes ->
Strong) does not fit VCEPs that let several LP comparators add up, or that count
LP like P. Each rule here is written as points per comparator (P / LP) and the
minimum total for each strength, which expresses every published combination
exactly; the CSpec wording is quoted next to each rule.

Comparators are counted as distinct variants: for PM5 one per different
amino-acid change at the codon (its best classification), for PS1 one per
different nucleotide change giving the same amino-acid change. ClinVar records
(>= the configured star level) stand in for "classified by the VCEP", as for the
other PS1 / PM5 rules.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Optional

from acmg_classifier.models.enums import CriterionStrength

S = CriterionStrength


@dataclass(frozen=True)
class TierRule:
    p: float                                    # points per Pathogenic comparator
    lp: float                                   # points per Likely pathogenic comparator
    tiers: tuple[tuple[CriterionStrength, float], ...]   # strongest first: (strength, min points)
    source: str                                 # VCEP / CSpec, for the evidence text
    cap: Optional[CriterionStrength] = None     # ceiling below the top tier (RASopathy)
    note: str = ""                              # added to the evidence when the top tier is capped


_LGMD = "LGMD VCEP v2.0.0"
_RAS_NOTE = ("RASopathy VCEP: PM5_Strong needs >=2 different (likely) pathogenic changes "
             "observed in >=5 probands; the proband count is not verifiable here, so PM5 is "
             "capped at Moderate - review the probands")

PM5_RULES: dict[str, TierRule] = {
    # GN019 MYOC 2.1.0: Strong "2 P, or 1 P and 2 LP, or 4 LP"; Moderate "1 P or
    # 2 LP"; Supporting "1 LP" (all assessed independently of PM5).
    "MYOC": TierRule(1.0, 0.5, ((S.STRONG, 2.0), (S.MODERATE, 1.0), (S.SUPPORTING, 0.5)),
                     "Glaucoma VCEP v2.1.0"),
    # GN021 ACADVL 2.1.0, PM5 table (05/01/25): "each additional variant at the same
    # codon: LP = 0.5, P = 1.0"; Strong >= 2.0, Moderate >= 1.0, Supporting 0.5.
    "ACADVL": TierRule(1.0, 0.5, ((S.STRONG, 2.0), (S.MODERATE, 1.0), (S.SUPPORTING, 0.5)),
                       "ACADVL VCEP v2.1.0 PM5 table"),
    # GN081 / GN090 VWF 1.0.0: Moderate "a previously reported variant reaches P ...
    # or two previously reported variants reach LP"; Supporting "LP". No Strong.
    "VWF": TierRule(2.0, 1.0, ((S.MODERATE, 2.0), (S.SUPPORTING, 1.0)), "VWD VCEP v1.0.0"),
    # GN071 F8 2.0.0 / GN080 F9 2.1.0: Moderate "1 pathogenic variant or 2 likely
    # pathogenic variants at the same residue"; Supporting "1 likely pathogenic".
    "F8": TierRule(2.0, 1.0, ((S.MODERATE, 2.0), (S.SUPPORTING, 1.0)), "CFD VCEP v2.0.0"),
    "F9": TierRule(2.0, 1.0, ((S.MODERATE, 2.0), (S.SUPPORTING, 1.0)), "CFD VCEP v2.1.0"),
    # GN112 KCNQ1 (LQTS VCEP): Moderate only; "Comparison variants must have a
    # Pathogenic or Likely Pathogenic classification"; "Not used at the PM5_Strong
    # level" -> one P or LP comparator gives Moderate.
    "KCNQ1": TierRule(1.0, 1.0, ((S.MODERATE, 1.0),), "LQTS VCEP"),
}
# GN180/184-189 LGMD 2.0.0: Strong "2 pathogenic or 3 likely pathogenic variants
# resulting in different amino acid changes at the same residue"; Moderate "1 P or
# 2 LP"; Supporting "1 LP".
for _g in ("DYSF", "SGCA", "SGCB", "SGCG", "SGCD", "CAPN3", "ANO5"):
    PM5_RULES[_g] = TierRule(3.0, 2.0, ((S.STRONG, 6.0), (S.MODERATE, 3.0),
                                        (S.SUPPORTING, 2.0)), _LGMD)
# RASopathy VCEP v2.3.0 (GN038-049, GN087, GN094 LZTR1, GN127/128): Strong ">=2
# different [likely] pathogenic residue changes at the same codon observed in >=5
# probands"; Moderate "1 [likely] pathogenic residue change at the same codon".
# LP counts like P; Strong is capped (probands not available) with a note.
for _g in ("PTPN11", "SOS1", "SOS2", "RAF1", "KRAS", "NRAS", "HRAS", "MAP2K1", "MAP2K2",
           "BRAF", "RIT1", "SHOC2", "MRAS", "RRAS2", "PPP1CB", "LZTR1"):
    PM5_RULES[_g] = TierRule(1.0, 1.0, ((S.STRONG, 2.0), (S.MODERATE, 1.0)),
                             "RASopathy VCEP", cap=S.MODERATE, note=_RAS_NOTE)

# PS1: "1 pathogenic or 2 likely pathogenic variants resulting in the same amino acid
# change" -> Strong, "1 likely pathogenic" -> Moderate (LGMD GN180/184-189, F8 GN071,
# F9 GN080; ACADVL GN021: "a P variant ... OR >=2 LP variants").
_PS1_1P_2LP = ((S.STRONG, 2.0), (S.MODERATE, 1.0))
PS1_RULES: dict[str, TierRule] = {
    "ACADVL": TierRule(2.0, 1.0, _PS1_1P_2LP, "ACADVL VCEP v2.1.0"),
    "F8": TierRule(2.0, 1.0, _PS1_1P_2LP, "CFD VCEP v2.0.0"),
    "F9": TierRule(2.0, 1.0, _PS1_1P_2LP, "CFD VCEP v2.1.0"),
}
for _g in ("DYSF", "SGCA", "SGCB", "SGCG", "SGCD", "CAPN3", "ANO5"):
    PS1_RULES[_g] = TierRule(2.0, 1.0, _PS1_1P_2LP, _LGMD)

_P = {"pathogenic", "pathogenic/likely pathogenic"}
_LP = {"likely pathogenic"}
_RANK = {S.SUPPORTING: 1, S.MODERATE: 2, S.STRONG: 3, S.VERY_STRONG: 4}


def classify_comparator(sig: Optional[str]) -> Optional[str]:
    """'P', 'LP' or None for a ClinVar clinical significance."""
    s = (sig or "").strip().lower()
    if s in _P:
        return "P"
    if s in _LP:
        return "LP"
    return None


def count_distinct(hits: Iterable, key: Callable) -> tuple[int, int]:
    """(#P, #LP) over distinct comparators; a comparator seen as both P and LP
    (several ClinVar records of the same change) counts once as P."""
    best: dict = {}
    for h in hits:
        cls = classify_comparator(getattr(h, "clinical_significance", None))
        if cls is None:
            continue
        k = key(h)
        if k is None:
            k = id(h)
        if best.get(k) != "P":
            best[k] = cls
    n_p = sum(1 for v in best.values() if v == "P")
    return n_p, len(best) - n_p


def tier_strength(rule: TierRule, n_p: int, n_lp: int
                  ) -> tuple[Optional[CriterionStrength], str]:
    """Strength for ``n_p`` / ``n_lp`` distinct comparators under ``rule`` (None =
    below the lowest tier) and the evidence text."""
    pts = n_p * rule.p + n_lp * rule.lp
    strength = next((st for st, need in rule.tiers if pts >= need), None)
    text = f"{n_p} P + {n_lp} LP distinct comparator(s) = {pts:g} points ({rule.source})"
    if strength is not None and rule.cap is not None and _RANK[strength] > _RANK[rule.cap]:
        strength = rule.cap
        if rule.note:
            text += f" [{rule.note}]"
    return strength, text
