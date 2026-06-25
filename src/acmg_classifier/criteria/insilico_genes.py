"""Per-gene auxiliary in-silico PP3/BP4 rules (BayesDel / CADD / 2-of-3 combos).

Some VCEPs define PP3/BP4 on predictors beyond the single tool chosen via
``--insilico-tool``: BayesDel (ENIGMA BRCA1/2, TP53) or an *agreement* of REVEL
with CADD / AlphaMissense (Antibody-Deficiency CTLA4/PIK3CD/PIK3R1, Pulmonary-
Hypertension BMPR2, ABCA4). Those predictors are licence-encumbered, so they are
loaded only under ``--insilico-tool revel/alphamissense`` together with
``--with-bayesdel`` / ``--with-cadd`` (see ``AnnotationOrchestrator``).

This module holds the cutoffs (mined from the ClinGen cspec PP3/BP4 texts) and,
when active, is **authoritative** for the gene: it REPLACES the genome-wide
single-tool PP3/BP4 path so that, e.g., CTLA4 cannot meet PP3 on REVEL alone when
its VCEP requires REVEL∧CADD agreement (the per-gene ``revel_*`` cutoffs that the
default path would otherwise apply are deliberately superseded here).

Each rule returns:
  * ``(strength, note)`` — criterion met at that strength (authoritative);
  * ``(None, note)``     — criterion authoritatively NOT met (suppress default);
  * ``None``             — rule does not govern this consequence; the caller
                           falls through to the default single-tool path.

The caller (PP3/BP4 evaluators) only consults a gene when its required auxiliary
predictor is enabled, so a return value is always safe to treat as final.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from acmg_classifier.models.enums import ConsequenceType, CriterionStrength

# Outcome: None = fall through; (strength|None, note) = authoritative met/not-met.
Outcome = Optional[tuple[Optional[CriterionStrength], str]]

_SUP = CriterionStrength.SUPPORTING
_MOD = CriterionStrength.MODERATE

# Which auxiliary predictor a gene's rule needs loaded to activate. The caller
# checks this against cfg.use_bayesdel / cfg.use_cadd before consulting the rule.
_REQUIRES: dict[str, str] = {
    "BRCA1": "bayesdel",
    "BRCA2": "bayesdel",
    "TP53": "bayesdel",
    "CTLA4": "cadd",
    "PIK3CD": "cadd",
    "PIK3R1": "cadd",
    "BMPR2": "cadd",
    "ABCA4": "cadd",
}

# ENIGMA "(potentially) clinically important functional domains" (1-based aa
# ranges) — PP3/BP4 BayesDel only apply to variants inside these.
_BRCA1_DOMAINS = ((2, 101), (1391, 1424), (1650, 1857))
_BRCA2_DOMAINS = ((10, 40), (2481, 3186))

_INFRAME = (ConsequenceType.INFRAME_INSERTION, ConsequenceType.INFRAME_DELETION)


@dataclass(frozen=True)
class InSilicoScores:
    """The predictor scores + variant context a per-gene rule needs."""

    consequence: ConsequenceType
    protein_position: Optional[int]
    revel: Optional[float]
    alphamissense: Optional[float]
    cadd: Optional[float]
    bayesdel: Optional[float]
    splice_max: Optional[float]  # max splice delta (SpliceAI/OpenSpliceAI) or None
    hgvs_c: Optional[str] = None  # transcript change (TP53 precomputed-code key)
    hgvs_p: Optional[str] = None  # protein change (TP53 fallback key)


def _f(v: Optional[float]) -> str:
    return f"{v:.3f}" if v is not None else "NA"


def _in_domain(pp: Optional[int], domains: tuple[tuple[int, int], ...]) -> bool:
    return pp is not None and any(lo <= pp <= hi for lo, hi in domains)


def _splice_impact(sc: InSilicoScores, thr: float) -> bool:
    """True when a splice predictor reports impact at/above ``thr`` (benign-side
    gate). Unknown splice (no predictor) is treated as no impact — graceful
    degradation consistent with the rest of the pipeline."""
    return sc.splice_max is not None and sc.splice_max >= thr


# ---------------------------------------------------------------------------
# BayesDel genes (ENIGMA BRCA1/2, TP53)
# ---------------------------------------------------------------------------


def _brca_pp3(sc: InSilicoScores, domains, cutoff: float) -> Outcome:
    # ENIGMA BayesDel PP3 is scoped to missense in a clinically-important domain.
    if sc.consequence != ConsequenceType.MISSENSE and sc.consequence not in _INFRAME:
        return None
    if not _in_domain(sc.protein_position, domains):
        return (None, "BayesDel PP3 N/A (outside clinically-important domain)")
    if sc.bayesdel is None:
        return (None, "BayesDel PP3 not met (no BayesDel score)")
    if sc.bayesdel >= cutoff:
        return (_SUP, f"BayesDel no-AF={sc.bayesdel:.3f} >= {cutoff} in functional domain (Supporting)")
    return (None, f"BayesDel no-AF={sc.bayesdel:.3f} < {cutoff} (PP3 not met)")


def _brca_bp4(sc: InSilicoScores, domains, cutoff: float) -> Outcome:
    if sc.consequence != ConsequenceType.MISSENSE and sc.consequence not in _INFRAME:
        return None
    if not _in_domain(sc.protein_position, domains):
        return (None, "BayesDel BP4 N/A (outside clinically-important domain)")
    if sc.bayesdel is None:
        return (None, "BayesDel BP4 not met (no BayesDel score)")
    if _splice_impact(sc, 0.2):
        return (None, f"BP4 blocked by predicted splice impact (delta={_f(sc.splice_max)})")
    if sc.bayesdel <= cutoff:
        return (_SUP, f"BayesDel no-AF={sc.bayesdel:.3f} <= {cutoff} in functional domain, no splice impact (Supporting)")
    return (None, f"BayesDel no-AF={sc.bayesdel:.3f} > {cutoff} (BP4 not met)")


# The TP53 VCEP gates PP3/BP4 on the Align-GVGD class (which this pipeline does
# not compute) combined with BayesDel. Rather than reimplement aGVGD, we consult
# the VCEP's published per-missense code table (TP53Codes), which already bakes in
# the aGVGD class. The table column is the protein-level ("missense-only") code;
# predicted splice impact is handled separately by the PP3/BP4 splice branches.
_TP53_PP3_STRENGTH = {"PP3": _SUP, "PP3_moderate": _MOD}
_TP53_BP4_STRENGTH = {"BP4": _SUP, "BP4_moderate": _MOD}


def _norm_hgvs(h: Optional[str]) -> Optional[str]:
    """Strip a transcript/protein accession prefix: 'NM_000546.6:c.4G>C' -> 'c.4G>C'."""
    if not h:
        return None
    return h.split(":", 1)[1].strip() if ":" in h else h.strip()


@dataclass(frozen=True)
class TP53Entry:
    """One VCEP precomputed-code row: the call plus the evidence behind it."""

    code: str                 # PP3 / PP3_moderate / BP4 / BP4_moderate / No evidence
    agvgd: str                # Align-GVGD class (e.g. "Class C65")
    bayesdel: str             # BayesDel score (string, as published)


class TP53Codes:
    """ClinGen TP53 VCEP precomputed PP3/BP4 missense codes, loaded from the TSV
    built by scripts/build_tp53_codes.py.

    Keyed primarily on the transcript change (hgvs_c, unique per variant); the
    protein change (hgvs_p) is a fallback, populated only for unambiguous protein
    changes (a few aa changes map to >1 code via nucleotide-level BayesDel, so they
    are resolvable by hgvs_c only). Each entry carries the Align-GVGD class and
    BayesDel score so a MET call can explain *why* in its evidence. A missing file
    degrades to "no codes"."""

    def __init__(self, tsv_path: Path) -> None:
        self._by_c: dict[str, TP53Entry] = {}
        self._by_p: dict[str, Optional[TP53Entry]] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        try:
            if not tsv_path.exists():
                return
            with tsv_path.open(encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    c = (row.get("hgvs_c") or "").strip()
                    p = (row.get("hgvs_p") or "").strip()
                    code = (row.get("code") or "").strip()
                    if not c or not code:
                        continue
                    entry = TP53Entry(
                        code=code,
                        agvgd=(row.get("agvgd") or "").strip(),
                        bayesdel=(row.get("bayesdel") or "").strip(),
                    )
                    self._by_c[c] = entry
                    if p:
                        # Mark a protein change ambiguous (None) once it maps to a
                        # second, different code.
                        if p in self._by_p and (self._by_p[p] is None
                                                or self._by_p[p].code != code):
                            self._by_p[p] = None
                        elif p not in self._by_p:
                            self._by_p[p] = entry
        except Exception:
            # A malformed/locked file must not break evaluator construction.
            self._by_c.clear()
            self._by_p.clear()

    def __bool__(self) -> bool:
        return bool(self._by_c)

    def lookup(self, hgvs_c: Optional[str], hgvs_p: Optional[str]) -> Optional[TP53Entry]:
        c = _norm_hgvs(hgvs_c)
        if c and c in self._by_c:
            return self._by_c[c]
        p = _norm_hgvs(hgvs_p)
        if p:
            return self._by_p.get(p)  # None when absent or ambiguous
        return None


def _tp53_evidence(entry: TP53Entry, strength: CriterionStrength) -> str:
    return (f"TP53 VCEP precomputed code '{entry.code}' "
            f"(Align-GVGD={entry.agvgd or 'NA'}, BayesDel={entry.bayesdel or 'NA'}) "
            f"({strength.value})")


def _tp53_pp3(sc: InSilicoScores, codes: Optional[TP53Codes]) -> Outcome:
    if not codes:
        return (None, "TP53 PP3 unavailable (VCEP code table not loaded)")
    if sc.consequence != ConsequenceType.MISSENSE:
        return None
    entry = codes.lookup(sc.hgvs_c, sc.hgvs_p)
    if entry is None:
        return (None, "TP53 PP3 not met (variant not in VCEP code table)")
    strength = _TP53_PP3_STRENGTH.get(entry.code)
    if strength:
        return (strength, _tp53_evidence(entry, strength))
    return (None, f"TP53 VCEP code '{entry.code}' "
                  f"(Align-GVGD={entry.agvgd or 'NA'}, BayesDel={entry.bayesdel or 'NA'}) — PP3 not met")


def _tp53_bp4(sc: InSilicoScores, codes: Optional[TP53Codes]) -> Outcome:
    if not codes:
        return (None, "TP53 BP4 unavailable (VCEP code table not loaded)")
    if sc.consequence != ConsequenceType.MISSENSE:
        return None
    entry = codes.lookup(sc.hgvs_c, sc.hgvs_p)
    if entry is None:
        return (None, "TP53 BP4 not met (variant not in VCEP code table)")
    strength = _TP53_BP4_STRENGTH.get(entry.code)
    if strength:
        return (strength, _tp53_evidence(entry, strength))
    return (None, f"TP53 VCEP code '{entry.code}' "
                  f"(Align-GVGD={entry.agvgd or 'NA'}, BayesDel={entry.bayesdel or 'NA'}) — BP4 not met")


# ---------------------------------------------------------------------------
# CADD two-tool-agreement genes (CTLA4 / PIK3CD / PIK3R1)
# ---------------------------------------------------------------------------


def _two_tool_pp3(sc: InSilicoScores, revel_min: float, cadd_min: float,
                  gene: str) -> Outcome:
    if sc.consequence != ConsequenceType.MISSENSE:
        return None
    r_ok = sc.revel is not None and sc.revel >= revel_min
    c_ok = sc.cadd is not None and sc.cadd >= cadd_min
    if r_ok and c_ok:
        return (_SUP, f"{gene}: REVEL={sc.revel:.3f}>={revel_min} ∧ CADD={sc.cadd:.1f}>={cadd_min} "
                      f"(2-tool agreement, Supporting)")
    return (None, f"{gene} PP3 not met (needs REVEL>={revel_min} ∧ CADD>={cadd_min}; "
                  f"REVEL={_f(sc.revel)} CADD={_f(sc.cadd)})")


def _two_tool_bp4(sc: InSilicoScores, revel_max: float, cadd_max: float,
                  gene: str, strict: bool) -> Outcome:
    if sc.consequence != ConsequenceType.MISSENSE:
        return None
    if _splice_impact(sc, 0.1):
        return (None, f"{gene} BP4 blocked by predicted splice impact (delta={_f(sc.splice_max)})")
    if strict:
        r_ok = sc.revel is not None and sc.revel < revel_max
        c_ok = sc.cadd is not None and sc.cadd < cadd_max
        op = "<"
    else:
        r_ok = sc.revel is not None and sc.revel <= revel_max
        c_ok = sc.cadd is not None and sc.cadd <= cadd_max
        op = "<="
    if r_ok and c_ok:
        return (_SUP, f"{gene}: REVEL={sc.revel:.3f}{op}{revel_max} ∧ CADD={sc.cadd:.1f}{op}{cadd_max} "
                      f"(2-tool agreement, Supporting)")
    return (None, f"{gene} BP4 not met (needs REVEL{op}{revel_max} ∧ CADD{op}{cadd_max}; "
                  f"REVEL={_f(sc.revel)} CADD={_f(sc.cadd)})")


# ---------------------------------------------------------------------------
# BMPR2 — 2-of-3 (CADD / AlphaMissense / REVEL)
# ---------------------------------------------------------------------------


def _bmpr2_pp3(sc: InSilicoScores) -> Outcome:
    if sc.consequence != ConsequenceType.MISSENSE:
        return None
    passes = sum((
        sc.cadd is not None and sc.cadd >= 25.3,
        sc.alphamissense is not None and sc.alphamissense >= 0.792,
        sc.revel is not None and sc.revel >= 0.644,
    ))
    if passes >= 2:
        return (_SUP, f"BMPR2: {passes}/3 of CADD>=25.3/AM>=0.792/REVEL>=0.644 "
                      f"(CADD={_f(sc.cadd)} AM={_f(sc.alphamissense)} REVEL={_f(sc.revel)}, Supporting)")
    return (None, f"BMPR2 PP3 not met (<2 of 3 pass; CADD={_f(sc.cadd)} "
                  f"AM={_f(sc.alphamissense)} REVEL={_f(sc.revel)})")


def _bmpr2_bp4(sc: InSilicoScores) -> Outcome:
    # Synonymous: CADD <= 22.7 alone suffices (REVEL/AM not applicable).
    if sc.consequence == ConsequenceType.SYNONYMOUS:
        if _splice_impact(sc, 0.1):
            return (None, f"BMPR2 synonymous BP4 blocked by splice impact (delta={_f(sc.splice_max)})")
        if sc.cadd is not None and sc.cadd <= 22.7:
            return (_SUP, f"BMPR2 synonymous CADD={sc.cadd:.1f}<=22.7 (Supporting)")
        return (None, f"BMPR2 synonymous BP4 not met (CADD={_f(sc.cadd)})")
    if sc.consequence != ConsequenceType.MISSENSE:
        return None
    if _splice_impact(sc, 0.1):
        return (None, f"BMPR2 BP4 blocked by predicted splice impact (delta={_f(sc.splice_max)})")
    passes = sum((
        sc.cadd is not None and sc.cadd <= 22.7,
        sc.alphamissense is not None and sc.alphamissense <= 0.169,
        sc.revel is not None and sc.revel <= 0.29,
    ))
    if passes >= 2:
        return (_SUP, f"BMPR2: {passes}/3 of CADD<=22.7/AM<=0.169/REVEL<=0.29 "
                      f"(CADD={_f(sc.cadd)} AM={_f(sc.alphamissense)} REVEL={_f(sc.revel)}, Supporting)")
    return (None, f"BMPR2 BP4 not met (<2 of 3 pass; CADD={_f(sc.cadd)} "
                  f"AM={_f(sc.alphamissense)} REVEL={_f(sc.revel)})")


# ---------------------------------------------------------------------------
# ABCA4 — REVEL for missense (default path), CADD for synonymous / in-frame indel
# ---------------------------------------------------------------------------


def _abca4_pp3(sc: InSilicoScores) -> Outcome:
    # Missense uses REVEL only — already covered by the default per-gene REVEL
    # path (revel_genes.py), so fall through.
    if sc.consequence == ConsequenceType.MISSENSE:
        return None
    if sc.consequence != ConsequenceType.SYNONYMOUS and sc.consequence not in _INFRAME:
        return None
    if sc.cadd is None:
        return (None, "ABCA4 PP3 not met (no CADD score)")
    if sc.cadd >= 28.1:
        return (_MOD, f"ABCA4 synonymous/indel CADD={sc.cadd:.1f}>=28.1 (Moderate)")
    if sc.cadd >= 25.3:
        return (_SUP, f"ABCA4 synonymous/indel CADD={sc.cadd:.1f} in [25.3, 28.0] (Supporting)")
    return (None, f"ABCA4 PP3 not met (CADD={sc.cadd:.1f} < 25.3)")


def _abca4_bp4(sc: InSilicoScores) -> Outcome:
    if sc.consequence == ConsequenceType.MISSENSE:
        return None
    if sc.consequence != ConsequenceType.SYNONYMOUS and sc.consequence not in _INFRAME:
        return None
    if sc.cadd is None:
        return (None, "ABCA4 BP4 not met (no CADD score)")
    if sc.cadd <= 17.3:
        return (_MOD, f"ABCA4 synonymous/indel CADD={sc.cadd:.1f}<=17.3 (Moderate)")
    if sc.cadd <= 20.0:
        return (_SUP, f"ABCA4 synonymous/indel CADD={sc.cadd:.1f} in [17.4, 20.0] (Supporting)")
    return (None, f"ABCA4 BP4 not met (CADD={sc.cadd:.1f} > 20.0)")


class InSilicoGeneSpec:
    """Registry of per-gene auxiliary PP3/BP4 rules.

    ``tp53_codes`` is the optional ClinGen TP53 VCEP precomputed-code table; when
    absent, TP53 PP3/BP4 degrade to not-applied."""

    def __init__(self, tp53_codes: Optional[TP53Codes] = None) -> None:
        self._tp53 = tp53_codes

    def requires(self, gene: Optional[str]) -> Optional[str]:
        """The auxiliary predictor ('bayesdel' / 'cadd') the gene's rule needs,
        or None when no per-gene auxiliary rule exists."""
        return _REQUIRES.get(gene) if gene else None

    def covers(self, gene: Optional[str]) -> bool:
        return bool(gene) and gene in _REQUIRES

    def pp3(self, gene: str, sc: InSilicoScores) -> Outcome:
        if gene == "BRCA1":
            return _brca_pp3(sc, _BRCA1_DOMAINS, 0.28)
        if gene == "BRCA2":
            return _brca_pp3(sc, _BRCA2_DOMAINS, 0.30)
        if gene == "TP53":
            return _tp53_pp3(sc, self._tp53)
        if gene == "CTLA4":
            return _two_tool_pp3(sc, 0.75, 20.0, gene)
        if gene == "PIK3CD":
            return _two_tool_pp3(sc, 0.644, 25.3, gene)
        if gene == "PIK3R1":
            return _two_tool_pp3(sc, 0.644, 26.0, gene)
        if gene == "BMPR2":
            return _bmpr2_pp3(sc)
        if gene == "ABCA4":
            return _abca4_pp3(sc)
        return None

    def bp4(self, gene: str, sc: InSilicoScores) -> Outcome:
        if gene == "BRCA1":
            return _brca_bp4(sc, _BRCA1_DOMAINS, 0.15)
        if gene == "BRCA2":
            return _brca_bp4(sc, _BRCA2_DOMAINS, 0.18)
        if gene == "TP53":
            return _tp53_bp4(sc, self._tp53)
        if gene == "CTLA4":
            return _two_tool_bp4(sc, 0.25, 20.0, gene, strict=True)
        if gene == "PIK3CD":
            return _two_tool_bp4(sc, 0.290, 22.7, gene, strict=False)
        if gene == "PIK3R1":
            return _two_tool_bp4(sc, 0.290, 21.5, gene, strict=False)
        if gene == "BMPR2":
            return _bmpr2_bp4(sc)
        if gene == "ABCA4":
            return _abca4_bp4(sc)
        return None


def build_scores(annotation, pc) -> InSilicoScores:
    """Collect the predictor scores a per-gene rule needs from an AnnotationData
    and its primary ConsequenceInfo."""
    sp = annotation.splice
    splice_max = sp.max_delta if (sp and sp.is_available and sp.max_delta is not None) else None
    return InSilicoScores(
        consequence=pc.consequence,
        protein_position=pc.protein_position,
        revel=annotation.revel.score if annotation.revel else None,
        alphamissense=annotation.alphamissense.score if annotation.alphamissense else None,
        cadd=annotation.cadd.phred if annotation.cadd else None,
        bayesdel=annotation.bayesdel.score if annotation.bayesdel else None,
        splice_max=splice_max,
        hgvs_c=pc.hgvs_c,
        hgvs_p=pc.hgvs_p,
    )


def combo_active(spec: InSilicoGeneSpec, gene: Optional[str], cfg) -> bool:
    """True when the gene's auxiliary rule should govern PP3/BP4 for this run.

    Gated to the licence-encumbered missense path (insilico_tool REVEL/
    AlphaMissense) AND the required auxiliary predictor being enabled — so it is
    NEVER active under ESM1B (where BayesDel/CADD are not loaded at all)."""
    from acmg_classifier.models.enums import InSilicoTool
    if cfg.insilico_tool not in (InSilicoTool.REVEL, InSilicoTool.ALPHAMISSENSE):
        return False
    req = spec.requires(gene)
    return (req == "cadd" and cfg.use_cadd) or (req == "bayesdel" and cfg.use_bayesdel)
