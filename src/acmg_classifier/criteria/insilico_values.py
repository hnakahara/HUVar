"""Score / cutoff notes for PP3 and BP4 (reported in the TSV's inference columns).

The tier tables mirror the strength functions in pathogenic/pp3.py and
benign/bp4.py (Bergquist 2024 Table 2; FuncVEP-CTI: Kayaalp 2026); a unit test
checks every cutoff against those functions so the two cannot drift apart.
"""
from __future__ import annotations

from acmg_classifier.models.enums import CriterionStrength, InSilicoTool

S = CriterionStrength
_ORDER = [S.SUPPORTING, S.MODERATE, S.THREE_POINT, S.STRONG, S.VERY_STRONG]

AM_PP3 = {S.STRONG: 0.990, S.THREE_POINT: 0.972, S.MODERATE: 0.906, S.SUPPORTING: 0.792}
AM_BP4 = {S.THREE_POINT: 0.070, S.MODERATE: 0.099, S.SUPPORTING: 0.169}
ESM1B_PP3 = {S.STRONG: -24.0, S.THREE_POINT: -14.0, S.MODERATE: -12.2, S.SUPPORTING: -10.7}
ESM1B_BP4 = {S.THREE_POINT: 8.8, S.MODERATE: -3.2, S.SUPPORTING: -6.3}


def tiers_text(table: dict, op: str) -> str:
    """"Supporting>=0.6903; Moderate>=0.8023; …" (weakest tier first)."""
    return "; ".join(f"{s.value}{op}{table[s]:g}" for s in _ORDER if s in table)


def missense_tool_values(cfg, annotation, gene: str | None, revel_spec, side: str) -> dict:
    """Tool, score and tier cutoffs of the configured missense predictor.

    ``side`` is "pp3" or "bp4". Returns {} when the tool has no score."""
    from acmg_classifier.criteria.benign.bp4 import _FUNCVEP_BP4, _REVEL_BP4_DEFAULT
    from acmg_classifier.criteria.pathogenic.pp3 import _FUNCVEP_PP3, _REVEL_PP3_DEFAULT

    pp3 = side == "pp3"
    tool = cfg.insilico_tool
    key = f"{side}_insilico_cutoffs"
    if tool == InSilicoTool.REVEL:
        rv = annotation.revel
        if not (rv and rv.score is not None):
            return {"insilico_tool": "REVEL"}
        rule = revel_spec.get(gene) if gene else None
        tiers = (rule.pp3 if pp3 else rule.bp4) if rule else None
        label = f"REVEL [{gene} VCEP cutoffs]" if tiers else "REVEL"
        table = tiers or (_REVEL_PP3_DEFAULT if pp3 else _REVEL_BP4_DEFAULT)
        return {"insilico_tool": label, "insilico_score": f"{rv.score:.3f}",
                key: tiers_text(table, ">=" if pp3 else "<=")}
    if tool == InSilicoTool.ESM1B:
        es = annotation.esm1b
        if not (es and es.llr is not None):
            return {"insilico_tool": "ESM1b"}
        return {"insilico_tool": "ESM1b", "insilico_score": f"{es.llr:.3f}",
                key: tiers_text(ESM1B_PP3 if pp3 else ESM1B_BP4, "<=" if pp3 else ">=")}
    if tool == InSilicoTool.FUNCVEP:
        fv = annotation.funcvep
        if not (fv and fv.score is not None):
            return {"insilico_tool": "FuncVEP-CTI"}
        return {"insilico_tool": "FuncVEP-CTI", "insilico_score": f"{fv.score:.4f}",
                key: tiers_text(_FUNCVEP_PP3 if pp3 else _FUNCVEP_BP4, ">=" if pp3 else "<=")}
    am = annotation.alphamissense
    if not (am and am.score is not None):
        return {"insilico_tool": "AlphaMissense"}
    return {"insilico_tool": "AlphaMissense", "insilico_score": f"{am.score:.3f}",
            key: tiers_text(AM_PP3 if pp3 else AM_BP4, ">=" if pp3 else "<=")}
