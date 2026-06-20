#!/usr/bin/env python3
"""Build disease_prevalence.tsv from ClinGen cspec GN*.json exports.

Each GN*.json is a single ClinGen VCEP criteria specification (JSON-LD from
cspec.genome.network). This script reads every `GN*.json` in a directory,
extracts the gene(s), mode(s) of inheritance, and the BA1 / BS1 filtering
allele-frequency thresholds from the *applicable* evidence-strength
descriptions, and writes one row per gene to a tab-separated
disease_prevalence.tsv consumable by `acmg_classifier.criteria.allele_frequency`.

Usage:
    python scripts/build_disease_thresholds.py \
        --json-dir resources/clingen \
        --out resources/shared/disease_prevalence.tsv

Notes / safety:
  * Only thresholds that are explicitly *Applicable* in the spec are emitted.
  * If a gene appears in several specs, the most conservative (highest) BA1 is
    kept and the conflict is recorded in `notes` for manual review.
  * Descriptions that carry more than one numeric cutoff (e.g. separate AR/AD
    values in one string) are flagged in `notes` — verify those rows by hand.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
from typing import Optional

COLUMNS = [
    "gene_symbol", "inheritance", "prevalence", "allelic_het", "genetic_het",
    "penetrance", "bs1_threshold", "bs1_strength", "bs1_exclude", "ba1_threshold", "af_basis",
    "pm2_threshold", "pm2_strength", "pm2_basis", "pm2_subpop", "pm2_zygosity",
    "pm4", "pm4_supporting_max_aa", "pp2",
    "pp2_requires", "pm5_grantham", "pm5_excludes", "pm5_max", "pm5_lp", "bs2",
    "bs2_count", "bs2_strength", "bs2_female_only", "bs2_hom_only", "pvs1",
    "ps1", "ps1_splice", "ps1_max", "ps1_paralog_group", "ps1_paralog_strength",
    "bp1", "bp1_target", "bp1_exclude", "bp1_strength",
    "bp1_no_splice", "bp3", "bp3_regions", "bp7_phylop", "bp7_intronic",
    "bp4_splice_cutoff", "bp7_splice_cutoff",
    "revel_pp3_supporting", "revel_pp3_moderate", "revel_pp3_strong",
    "revel_bp4_supporting", "revel_bp4_moderate", "revel_bp4_strong",
    "source_vcep", "cspec_url", "notes",
]

# Spec-preference overrides: gene -> the ONLY spec (GN id) allowed to supply its
# values. Normally the most gene-specific spec wins, but a few genes have empty
# single-gene specs (In-Prep/Pilot, 0 criteriaCodes) that would shadow a fully
# populated grouped spec. Pinning the authoritative spec makes it win every
# criterion (and the base row) for that gene.
#   ITGA2B / ITGB3: GN011 (Platelet Disorders VCEP, Released, 28 codes) states
#   e.g. BP1 not_applicable ("missense variants are a known cause of disease");
#   the empty GN059/GN060/GN221-223 single-gene specs must NOT override it.
_FORCE_SPEC: dict[str, str] = {
    "ITGA2B": "GN011",
    "ITGB3": "GN011",
}

# PS1 paralogue / analogous-residue groups, transcribed from the specs (the group
# text is prose, so curated rather than mined). gene -> (sibling genes, fixed
# paralog strength). RASopathy GN004 grants PS1 at the analogous residue across
# each "highly analogous grouping" with the full (comparator-derived) strength;
# HBA2 GN173 grants only PS1_Moderate from its paralogue HBA1.
_PS1_PARALOG: dict[str, tuple[str, str]] = {
    "HRAS": ("NRAS,KRAS", ""), "NRAS": ("HRAS,KRAS", ""), "KRAS": ("HRAS,NRAS", ""),
    "MAP2K1": ("MAP2K2", ""), "MAP2K2": ("MAP2K1", ""),
    "SOS1": ("SOS2", ""), "SOS2": ("SOS1", ""),
    "HBA2": ("HBA1", "Moderate"),
}

# Genes whose BA1/BS1 spec defines the cutoff on the *male* allele frequency
# rather than the overall population FAF — detected from "in males" wording that
# qualifies the FREQUENCY (RPGR, RS1). Emitted in af_basis so the BA1/BS1
# evaluators compare against gnomAD AF_XY for these genes.
# NB: the bare token "hemizyg" must NOT be used — several X-linked specs that
# actually use the overall "Total Grpmax FAF" merely mention "hemizygotes" in a
# count rule (OTC, SLC6A8) or a prevalence note (ABCD1), and were wrongly
# flagged males, making BA1/BS1 compare AF_XY and under-fire.
_MALES_BASIS = re.compile(r"\bin males\b", re.IGNORECASE)

# A frequency cutoff is the number directly after a ≥ / > operator, optionally
# followed by '%'. Anchoring on the operator avoids the many *non-threshold*
# numbers in these free-text descriptions — allele counts ("≥2,000 alleles",
# ">2000"), penetrance/heterogeneity ("85%", "100%"), and CI ("99.99% CI") —
# which either lack the operator or resolve to a value outside (0, 1).
# A frequency number: a decimal (leading zero optional, e.g. ".5"), an integer,
# or scientific "N x 10^-k" (optionally HTML <sup>). Each may carry a '%'.
_SCI = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s*[x×]\s*10\s*(?:<sup>)?\s*(-?[0-9]+)")
# Operator-anchored: a ≥ / > symbol OR a spelled-out form, then a number(+%).
_OP_NUM = re.compile(
    r"(?:≥|≧|>=|&ge;|>|≳|greater than or equal to|greater than|above|at least)"
    r"\s*([0-9]*\.?[0-9]+)\s*(%?)",
    re.IGNORECASE,
)
# Parenthetical proportion, e.g. "(or 0.000333)", "(0.0000111)", "(0.0043%)".
_PAREN = re.compile(r"\(\s*(?:or\s*)?([0-9]*\.[0-9]+)\s*(%?)\s*\)")
# Any number(+%) — used only as a last resort on short descriptions.
_ANY = re.compile(r"([0-9]*\.?[0-9]+)\s*(%?)")
# Legacy ACMG-2015 stand-alone-benign boilerplate: "above 5% in [Exome
# Sequencing Project | 1000 Genomes | Exome Aggregation Consortium]". This 5%
# is a pre-gnomAD catch-all, NOT the VCEP's operative gnomAD cutoff. When a spec
# also lists a gnomAD-specific number (e.g. KCNQ1 BA1 "≥0.004"), the 5% must be
# dropped or it shadows the real, lower cutoff (causing BA1 to never fire in the
# 0.4–5% band). Removed only when something else remains to parse (RPGR-style
# specs whose *only* number is this males "5%" keep it via the fallback below).
_LEGACY_5PCT = re.compile(
    r"(?:is\s+)?(?:above|over|greater than(?:\s+or\s+equal\s+to)?|>=?|≥)\s*5\s*%"
    r"[^.\n]*?(?:Exome Sequencing Project|1000\s*Genomes|"
    r"Exome Aggregation Consortium|\bESP\b|\bExAC\b)[^.\n]*",
    re.IGNORECASE,
)
# Sub-population application rule (Rett/Angelman-like panels, e.g. GN032-037):
# "...present at ≥0.000083 (0.0083%) in any sub-population." This is the VCEP's
# operative gnomAD cutoff and takes precedence over a generic "Allele frequency
# above 0.05%" headline in the same description (the headline is ACMG
# boilerplate, exactly like KCNQ1's legacy "5%"). Captures the number (and an
# optional '%') immediately preceding "in any sub[-]population".
_SUBPOP = re.compile(
    r"(?:≥|≧|>=|&ge;|>|greater than or equal to)\s*"
    r"([0-9]*\.?[0-9]+)\s*(%?)\s*(?:\([^)]*\))?\s*in any sub-?population",
    re.IGNORECASE,
)


def _to_prop(num: str, pct: str) -> Optional[float]:
    """A numeric token (+ optional '%') as a proportion in (0, 1), else None."""
    try:
        v = float(num)
    except ValueError:
        return None
    if pct:
        v /= 100.0
    return v if 0.0 < v < 1.0 else None


def _collect_cands(desc: str) -> list[float]:
    """Frequency candidates in (0, 1) from a description, by parser tier."""
    _pct = _to_prop

    # Tier 1: scientific notation (very specific; e.g. "≥8 x 10^-3", "1.11 x 10^-5").
    sci = [a * 10 ** int(b) for a, b in
           ((float(m.group(1)), m.group(2)) for m in _SCI.finditer(desc))]
    cands = [v for v in sci if 0.0 < v < 1.0]

    # Tier 2: operator-anchored (symbol or spelled-out). Safe on long narrative
    # descriptions because penetrance/CI/heterogeneity percentages there are not
    # introduced by a ≥/> or "greater than" operator.
    if not cands:
        cands = [v for v in (_pct(n, p) for n, p in _OP_NUM.findall(desc)) if v is not None]

    # Tier 3: parenthetical proportion, e.g. "(or 0.000333)", "(0.0043%)".
    if not cands:
        cands = [v for v in (_pct(n, p) for n, p in _PAREN.findall(desc)) if v is not None]

    # Tier 4 (last resort): first frequency-valued token anywhere. Only reached
    # by short cutoff strings ("MAF cutoff of 0.2%."); narrative descriptions are
    # already resolved above, so stray penetrance %s are not a concern here.
    if not cands:
        cands = [v for v in (_pct(n, p) for n, p in _ANY.findall(desc)) if v is not None]
    return cands


def _applicable_strengths(code: dict) -> list[dict]:
    return [
        es for es in code.get("evidenceStrengths", [])
        if str(es.get("applicability", "")).lower().startswith("applic")
        and es.get("description")
    ]


def _threshold_from_desc(desc: str) -> tuple[Optional[float], bool]:
    """Return (threshold, multi_value_flag) parsed from a description.

    Collects every operator-anchored frequency in (0, 1) — converting a
    trailing '%' to a proportion — and returns the FIRST one (specs state the
    operative cutoff first; for multi-MOI descriptions the autosomal-recessive,
    i.e. higher/conservative, value is listed first). A second distinct value
    sets the flag so the row can be spot-checked.

    Exception — *range* descriptions ("Allele frequency between X and Y"):
    these state a BS1 band whose operative cutoff is the LOWER edge (BS1 fires
    when AF ≥ that edge, up to the BA1 ceiling). The textual order of the two
    bounds is inconsistent across specs (RUNX1 lists lower-first, RPE65 lists
    upper-first), so the lower bound is taken as min(...) rather than by
    position to avoid silently adopting the upper bound (an effective BS1 = BA1
    that never fires).

    Exception — *legacy 5%* boilerplate ("above 5% in ESP/1000 Genomes/ExAC"):
    dropped when a gnomAD-specific cutoff is also present, so the operative
    gnomAD number wins (e.g. KCNQ1 BA1 → 0.004, not 0.05). If that legacy clause
    is the *only* number (RPGR-style males "5%"), it is retained as a fallback.

    Exception — *sub-population* application rule ("present at ≥X in any
    sub-population"): X is the VCEP's operative gnomAD cutoff and wins outright
    over any generic "above 0.05%" headline (e.g. Rett/Angelman-like panels →
    BA1 0.000083, not 0.0005).
    """
    # Highest precedence: an explicit "≥X in any sub-population" gnomAD rule.
    sub = _SUBPOP.search(desc)
    if sub:
        v = _to_prop(sub.group(1), sub.group(2))
        if v is not None:
            return v, False

    # Prefer the description with the legacy 5%/ESP clause removed; fall back to
    # the original only if stripping leaves nothing parseable.
    cands = _collect_cands(_LEGACY_5PCT.sub(" ", desc))
    if not cands:
        cands = _collect_cands(desc)

    if not cands:
        return None, False
    # Range band ("between X and Y"): the BS1 cutoff is the lower edge,
    # independent of which bound is written first.
    if re.search(r"\bbetween\b", desc, re.IGNORECASE) and len(set(cands)) > 1:
        return min(cands), True
    return cands[0], len(set(cands)) > 1


def _criterion_threshold(rule_set: dict, label: str) -> tuple[Optional[float], str]:
    """Extract the BA1/BS1 threshold (and a note) for one criterion label."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != label:
            continue
        strengths = _applicable_strengths(code)
        if not strengths:
            return None, f"{label}: not applicable"
        # Prefer the canonical strength (BA1→Stand Alone, BS1→Strong); else first.
        preferred = "Stand Alone" if label == "BA1" else "Strong"
        chosen = next(
            (s for s in strengths if s.get("label") == preferred), strengths[0]
        )
        thr, multi = _threshold_from_desc(chosen.get("description", ""))
        note = ""
        if thr is None:
            note = f"{label}: could not parse threshold from desc"
        elif multi:
            note = f"{label}: multiple cutoffs in desc — verify"
        return thr, note
    return None, f"{label}: code absent"


# cspec strength label → CriterionStrength value emitted in the TSV.
_CSPEC_TO_STRENGTH = {
    "Very Strong": "VeryStrong", "Strong": "Strong",
    "Moderate": "Moderate", "Supporting": "Supporting",
}


def _bs1_strength(rule_set: dict) -> str:
    """The strength of the BS1 tier whose threshold was chosen by
    _criterion_threshold (same selection: prefer Strong, else the first
    applicable). Lets the evaluator emit the spec's tier strength instead of
    always defaulting to BS1_Strong — e.g. MYO15A/OTOF (GN023) whose only BS1
    tiers are Very Strong (>=0.3%) and Supporting. "" when BS1 is not applicable
    (the evaluator then uses the Strong default)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BS1":
            continue
        strengths = _applicable_strengths(code)
        if not strengths:
            return ""
        chosen = next(
            (s for s in strengths if s.get("label") == "Strong"), strengths[0]
        )
        return _CSPEC_TO_STRENGTH.get(chosen.get("label", ""), "")
    return ""


# --- PP2 gene-level applicability (from the per-VCEP PP2 criteria code) ---
# A VCEP marks PP2 Applicable / Not Applicable per gene. The applicability flag
# lives on the evidence strengths, but gene-specific decisions are *also* in the
# free-text description: a whole-description "not applicable" (KCNQ1: "Not
# applicable due to … z-score 1.83") or a single-gene exclusion (GN018:
# "applicable to MTOR, PIK3CA and AKT3 but not PIK3R2"). Both are parsed below.
# These heuristics are validated against the current Released specs.
_PP2_NOT_APPLICABLE = re.compile(r"not applicable", re.IGNORECASE)
# A *positive* "applicable" not preceded by "not" — its presence means a
# "not applicable" elsewhere is a scoped exception, not a blanket negation.
_PP2_POSITIVE_APPLICABLE = re.compile(r"(?<!not )applicable", re.IGNORECASE)


def _pp2_gene_excluded(desc: str, gene: str) -> bool:
    """True if *desc* singles out *gene* as not-applicable ("(but) not GENE",
    "not applicable to GENE"). The gene may be wrapped in markdown underscores."""
    pat = r"\bnot\s+(?:applicable\s+(?:to|for|in)\s+)?_?" + re.escape(gene) + r"_?\b"
    return re.search(pat, desc, re.IGNORECASE) is not None


def _pp2_blanket_negation(desc: str) -> bool:
    """True if the description negates PP2 for the whole rule set — a "not
    applicable" with no positive "applicable to …" clause to scope it."""
    return bool(_PP2_NOT_APPLICABLE.search(desc)) and not _PP2_POSITIVE_APPLICABLE.search(desc)


def _pp2_applicability(rule_set: dict) -> dict[str, str]:
    """{gene: "applicable"|"not_applicable"} from the rule set's PP2 code.

    A gene is `applicable` when the PP2 code has an Applicable strength and the
    description neither blanket-negates nor excludes that gene; `not_applicable`
    when the VCEP carries a PP2 code but declined it (no Applicable strength) or
    negated/excluded the gene. Genes whose VCEP has no PP2 code are omitted."""
    genes = [g.get("label") for g in rule_set.get("genes", []) if g.get("label")]
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PP2":
            continue
        applic = _applicable_strengths(code)
        if not applic:
            return {g: "not_applicable" for g in genes}
        desc = applic[0].get("description", "") or ""
        out: dict[str, str] = {}
        for g in genes:
            if _pp2_gene_excluded(desc, g) or _pp2_blanket_negation(desc):
                out[g] = "not_applicable"
            else:
                out[g] = "applicable"
        return out
    return {}


# A PP2 description may make the criterion *conditional* on other criteria
# (BMPR2/GN125: "PM2_supporting and PP3 must be met."). We capture the required
# ACMG codes so the registry can suppress PP2 unless they are also triggered.
_PP2_REQ_TRIGGER = re.compile(r"must be met|must also be|required|requires", re.IGNORECASE)
_ACMG_CODE = re.compile(
    r"\b(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])(?:_[a-z]+)?\b",
    re.IGNORECASE,
)


def _pp2_more_specific(cand: tuple[int, str, str], cur: tuple[int, str, str]) -> bool:
    """True if PP2 decision *cand* should replace *cur*. Fewer-genes (more
    gene-specific) spec wins; on a specificity tie the conservative
    "not_applicable" wins (minimises false-positive PP2)."""
    if cand[0] != cur[0]:
        return cand[0] < cur[0]
    return cand[1] == "not_applicable" and cur[1] != "not_applicable"


def _pp2_requires(rule_set: dict) -> str:
    """Comma-joined ACMG codes PP2 is conditional on (e.g. "PM2,PP3"), or ""."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PP2":
            continue
        applic = _applicable_strengths(code)
        if not applic:
            return ""
        desc = applic[0].get("description", "") or ""
        if not _PP2_REQ_TRIGGER.search(desc):
            return ""
        out: list[str] = []
        for m in _ACMG_CODE.finditer(desc):
            c = m.group(1).upper()
            if c != "PP2" and c not in out:
                out.append(c)
        return ",".join(out)
    return ""


# --- PM5 Grantham-distance gate (from the per-VCEP PM5 criteria code) ---
# A subset of VCEPs require PM5 to clear a Grantham-distance comparison against
# the same-codon pathogenic/likely-pathogenic comparator. The comparison
# operator is encoded per gene: "ge" (candidate >= comparator; "equal or
# greater/worse" wording, the common case) or "gt" (strictly greater — PIK3R1
# "higher Grantham score", RYR1 comparator "must be less than" the candidate).
_PM5_GT = re.compile(r"higher grantham|less than", re.IGNORECASE)


def _pm5_grantham_op(rule_set: dict) -> str:
    """The PM5 chemical-severity gate token for this rule set's gene(s):
    ``ge``/``gt`` (Grantham distance), ``blosum_le``/``blosum_lt`` (BLOSUM62
    similarity, e.g. PTEN), or "" when no applicable PM5 strength mentions a
    matrix gate.

    Grantham: an "equal" clause is inclusive (``ge``); a strict "higher/less
    than" phrasing with no "equal" yields ``gt``; the inclusive ``ge`` wins on a
    cross-strength conflict. BLOSUM62 runs the opposite direction (lower score =
    more severe): "equal to or less" → ``blosum_le``, strict "less" → ``blosum_lt``."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PM5":
            continue
        op = ""
        for es in _applicable_strengths(code):
            desc = es.get("description", "") or ""
            low = desc.lower()
            if "blosum" in low:
                # BLOSUM62 gate: "equal to or less than" → inclusive ``blosum_le``.
                return "blosum_le" if "equal" in low else "blosum_lt"
            if "grantham" not in low:
                continue
            if "equal" in low:
                return "ge"  # inclusive clause is the operative gate
            if _PM5_GT.search(desc):
                op = "gt"
            elif not op:
                op = "ge"  # Grantham-conditioned but no explicit operator → inclusive
        return op
    return ""


# A PM5 description may forbid combining PM5 with PM1 (RUNX1: "PM5 cannot be
# used if PM1 was applied"; RASopathy/Cardiomyopathy: "PM5 should not be combined
# with PM1") and, for DICER1, also PS1. "Not mutually exclusive with PM1" is the
# opposite (PM5 *may* combine) and must NOT be read as an exclusion.
_PM5_EXCL_PM1 = re.compile(
    r"(?:should|can|must|may)\s*not\s+be\s+(?:combined|applied|used)"
    r"[^.]*?\b(?:with|if|in combination with|in conjunction with)\b[^.]*?\bPM1\b"
    r"|cannot be used if PM1",
    re.IGNORECASE,
)
_PM5_EXCL_PS1 = re.compile(r"PM1 or PS1|PS1 or PM1", re.IGNORECASE)
_PM5_NOT_MUTEX = re.compile(r"not mutually exclusive with pm1", re.IGNORECASE)
_STRENGTH_RANK = {"Supporting": 1, "Moderate": 2, "Strong": 3,
                  "Very Strong": 4, "Stand Alone": 5}


def _pm5_excludes(rule_set: dict) -> str:
    """Comma-joined criteria PM5 may not be combined with (e.g. "PM1",
    "PM1,PS1"), or "". Enforced as a registry post-hoc PM5 suppression."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PM5":
            continue
        descs = " ".join(es.get("description", "") or "" for es in _applicable_strengths(code))
        if _PM5_NOT_MUTEX.search(descs) or not _PM5_EXCL_PM1.search(descs):
            return ""
        out = ["PM1"]
        if _PM5_EXCL_PS1.search(descs):
            out.append("PS1")
        return ",".join(out)
    return ""


def _pm5_max(rule_set: dict) -> str:
    """Per-spec PM5 strength ceiling: the highest applicable PM5 strength —
    "Supporting" (ATM/CDH1/PALB2), "Moderate", or "Strong" (VCEPs that grant
    PM5_Strong for >=2 different pathogenic missense at the codon — RUNX1,
    HNF1A, GCK, HNF4A, …); "" when PM5 is not applicable. Aggregated across
    specs in main() by MAX rank, so the most permissive ceiling wins."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PM5":
            continue
        ranks = [_STRENGTH_RANK[es["label"]] for es in _applicable_strengths(code)
                 if es.get("label") in _STRENGTH_RANK]
        if not ranks:
            return ""
        top = max(ranks)
        if top <= _STRENGTH_RANK["Supporting"]:
            return "Supporting"
        if top == _STRENGTH_RANK["Moderate"]:
            return "Moderate"
        return "Strong"  # Strong or higher applicable
    return ""


def _pm5_lp_comparator(rule_set: dict) -> str:
    """PM5 same-codon comparator-significance policy for this rule set:

    * ``"no"``  — PM5 offers no applicable Supporting strength (only Moderate /
      Strong), so per these specs the comparator must reach **Pathogenic**; a
      merely Likely-pathogenic comparator must NOT trigger PM5 (e.g. PTEN, VHL,
      KCNQ1, the RASopathy single-gene specs).
    * ``"yes"`` — a Supporting PM5 strength is applicable, so a Likely-pathogenic
      comparator may trigger PM5 (at Supporting; e.g. ABCA4, RPE65).
    * ``""``    — no applicable PM5 strength (skipped in cross-spec resolution).
    """
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PM5":
            continue
        applic = _applicable_strengths(code)
        if not applic:
            return ""
        labels = {es.get("label") for es in applic}
        return "yes" if "Supporting" in labels else "no"
    return ""


# BS2 (observed in a healthy individual) gene-level applicability. Our BS2
# evaluator is gnomAD-based, so a VCEP must be treated as not_applicable whenever
# its BS2 cannot be derived from gnomAD population counts:
#   (a) it explicitly declines / bars general-population data, e.g.
#         GN004 RASopathy : "general population data should not be used ..."
#         GN120 RPE65     : "Presence in databases such as gnomAD are not considered"
#         GN112 KCNQ1     : "Not applicable due to incomplete penetrance"
#   (b) it scores BS2 purely on clinical phenotyping / point-based family
#       evaluation with NO gnomAD-countable rule, e.g. the single-gene RASopathy
#       specs (GN039 NRAS etc., "-4 Points." inherited from GN004) and the
#       Fanconi-Anemia "points per proband" cancer specs (BRCA1/2, PALB2).
# A spec that DOES offer a homozygote/hemizygote/gnomAD count path stays
# applicable even if it also mentions points (e.g. APC ">=10 points OR >=2 in
# homozygous state"; PDHA1 ">=16 hemizygotes in gnomAD").
_BS2_DECLINE = re.compile(
    r"should not be used|not considered|not applicable due to",
    re.IGNORECASE,
)
_BS2_POINTS = re.compile(r"\bpoints?\b", re.IGNORECASE)
_BS2_GNOMAD_COUNTABLE = re.compile(r"homozyg|hemizyg|gnomad", re.IGNORECASE)

# Genes whose VCEP BS2 cannot be derived from gnomAD population counts and is
# therefore forced not_applicable (a gnomAD-based evaluator would FALSELY fire
# BS2 on a pathogenic variant). Two reasons, neither auto-detectable by regex
# (the disqualifying phrases co-occur with "homozygous"/"healthy adult", so the
# generic countable heuristic would wrongly keep them applicable):
#
# (1) The rule needs phase, a lab/functional assay, or a specific clinical
#     phenotype gnomAD does not carry:
#   HNF4A/HNF1A GN085/GN017 — normoglycemic + age 70+ (fasting glucose, age)
#   GCK      GN086 — fasting glucose < 100 mg/dL (lab)
#   LDLR     GN013 — well-phenotyped, untreated, normolipidemic (lipids, Tx history)
#   GAA      GN010 — normal GAA enzyme activity (assay)
#   IDUA     GN091 — IDUA enzyme activity in unaffected range (assay + reference)
#   ITGA2B/ITGB3/GP1BA/GP1BB/GP9 GN011/079/082/083 — unaffected proven by
#            aggregometry / flow cytometry + platelet count & size (function)
#   CDH1     GN007 — tumor-free (GC/DGC/SRC/LBC) + family history not HDGC
#   DICER1   GN024 — tumor-free through age 50 + PS4-ratio caveat (internal cohort)
#   TP53     GN009 — cancer-free females ≥60y from a single source (internal cohort)
#   VHL      GN078 — ≥65y, full phenotyping & screening for VHL cancers
#   PTEN     GN003 — homozygous in healthy / PHTS-unaffected (clinical PHTS exclusion)
#   SERPINC1 GN084 — normal antithrombin level > 0.8 IU/mL (lab value)
#   MLH1/MSH2/MSH6/PMS2 GN115/137/138/139 — co-occurrence IN TRANS with a known
#            pathogenic variant + Lynch-cancer phenotype (gnomAD has no phase)
#   Hearing Loss GN005/GN023 (CDH23, GJB2, MYO6, MYO7A, SLC26A4, TECTA, USH2A,
#            MYO15A, OTOF) — BIALLELIC with a known pathogenic variant (phase)
#   UBE3A    GN016/GN037 — unaffected het, parent-of-origin / internal data
#
# (2) The rule scores BS2 on an individual confirmed to be a "healthy/unaffected
#     adult" — a per-individual clinical status gnomAD's aggregate counts do not
#     establish (gnomAD is *presumed*-healthy, not phenotyped; it also includes
#     non-adult and disease-cohort samples). Same wording that disqualified RYR1.
#     Genes that INSTEAD explicitly sanction gnomAD homozygote counting (BMPR2,
#     PIK3R2: "≥3 homozygotes in gnomAD") are NOT listed — they stay applicable.
#   Congenital Myopathies: ACTA1, DNM2, NEB (GN147/148/146)
#   Epilepsy Na-channel:  SCN1A, SCN2A, SCN3A, SCN8A (GN067-070)
#   SCID:                 DCLRE1C (GN116)
#   Mito/CCDS:            GATM (GN025)
#   Hemoglobinopathy:     HBB, HBA2 (GN170/173)
#   Rett/Angelman-like:   FOXG1, TCF4 (GN035/032)
#   Other:                APC (GN089)
#
# NOT listed (reclassified to *applicable*): IL7R, RAG1, RAG2 (SCID GN119/123/
# 124), ADA (GN114), PAH (GN006), POLG, ETHE1 (Mito GN014), GAMT (CCDS GN026).
# Their VCEP BS2 has an explicit gnomAD homozygote-count path ("≥N homozygotes" /
# "observed in the homozygous state in a healthy adult"), so a gnomAD-derived BS2
# is legitimate and was previously a FALSE NEGATIVE. They are recessive (or
# mode-agnostic) genes, so the evaluator counts gnomAD homozygotes (nhomalt).
_BS2_CLINICAL_CONFIRMATION = frozenset({
    # batch 1
    "HNF4A", "RYR1", "LDLR", "GAA", "ITGA2B", "ITGB3",
    "CDH1", "UBE3A", "DICER1", "SERPINC1", "TP53",
    # batch 2 — reason (1): phase / lab / specific phenotype
    "HNF1A", "GCK", "IDUA", "GP1BA", "GP1BB", "GP9",
    "VHL", "PTEN", "MLH1", "MSH2", "MSH6", "PMS2",
    "CDH23", "GJB2", "MYO6", "MYO7A", "SLC26A4", "TECTA",
    "USH2A", "MYO15A", "OTOF",
    # batch 2 — reason (2): "healthy/unaffected adult" not gnomAD-confirmable
    "ACTA1", "DNM2", "NEB", "SCN1A", "SCN2A", "SCN3A", "SCN8A",
    "GATM", "HBB", "HBA2",
    "FOXG1", "TCF4", "APC",
    # NB: DCLRE1C (GN116) is NOT here — its SCID VCEP takes the homozygote count
    # "in gnomAD" (BS2_Strong >=3 / Supporting >=1 hom), so it is population-count
    # eligible; routed applicable with bs2_strength tiers.
})


# X-linked genes whose VCEP BS2 requires an INTERNAL, clinically-phenotyped
# cohort that gnomAD cannot supply — a documented unaffected adult male with
# functional/lab confirmation, or phenotyped unaffected het/hemi relatives:
#   RPGR  GN106 / RS1 GN126 — males >30y with eye exam + functional studies
#                             (normal ERG/FAF; "without retinoschisis")
#   F8    GN071 / F9 GN080  — male with normal factor VIII/IX activity (>40% IU)
#   SLC6A8 GN027            — male with documented normal creatine transport study
#   PDHA1 GN014             — well-characterized phenotype (explicitly "not just
#                             seen in database") / Pyruvate enzyme assay
#   CDKL5/MECP2/SLC9A6 GN016 — "N unaffected (related/unrelated) Het/Hemi" Rett-
#                             like individuals (phenotyped, internal/family data)
# gnomAD is presumed-healthy and unphenotyped (and, per the hemizygote-count
# limitation, an unreliable source for X-linked male counts), so a gnomAD-count
# BS2 would FALSELY fire on a pathogenic X-linked variant. Forced not_applicable.
_BS2_XLINKED_INTERNAL = frozenset({
    "CDKL5", "RS1", "PDHA1", "RPGR",
    "SLC9A6", "F9", "SLC6A8", "MECP2", "F8",
    # NB: IL2RG (GN129) is NOT here — its SCID VCEP takes the hemizygote count
    # ">=3 in gnomAD" (Strong) / ">=2 in gnomAD" (Supporting), so it is
    # gnomAD-eligible; routed applicable with bs2_strength tiers.
})


# BP7 genes whose VCEP declared evolutionary conservation NON-informative (so no
# phyloP gate), but whose cspec JSON-LD export does not carry that statement in
# any parseable criteriaCode description. _bp7_phylop() auto-detects the in-text
# phrasings ("conservation is not required/considered/informative", ...); this
# curated set covers genes whose policy lives only in the published spec prose.
# Forced to the "na" sentinel so BP7 skips the conservation gate, like the
# auto-detected genes.
#   Leber Congenital Amaurosis / early-onset Retinal Dystrophy VCEP — RPE65
#   (GN120), GUCY2D (GN167), AIPL1 (GN208). Their BP7/BP4 criteriaCode
#   descriptions specify only the SpliceAI gate and excluded splice positions;
#   the "Evolutionary conservation is not considered informative for application
#   of this code" statement is absent from the JSON-LD (it is in the spec text).
#   The only phyloP in these files is an unrelated PM4 rule (PhyloP>2.0 for a
#   conserved residue), so nothing is auto-extracted for BP7.
#   BRCA1 (GN092) / BRCA2 (GN097) — the ENIGMA VCEP BP7 is BP4-driven (REVEL /
#   SpliceAI) and states no phyloP nucleotide-conservation requirement, so
#   conservation is not applied; nothing is auto-extractable from the text.
#   RUNX1 (GN008) — "Conservation is no longer a requirement for BP7 ... (Cheung
#   2019; Walker 2023)" appears in the spec notes, not the JSON-LD criteriaCode.
#   MYOC (GN019) — BP7 is purely BP4-driven ("Apply to intronic/noncoding ...
#   variants if BP4 is met"), with no phyloP conservation requirement.
_BP7_CONSERVATION_NA = frozenset({
    "RPE65", "GUCY2D", "AIPL1",
    "BRCA1", "BRCA2",
    "RUNX1", "MYOC",
})


# BP7 genes whose intronic applicability is BROAD — any intronic position gated
# only by SpliceAI (no +7/-21 deep-distance restriction), equivalent to the
# "noncanonical" mode but phrased without the "except canonical splice sites"
# wording the auto-detector (_bp7_intronic) keys on. Manually verified from the
# BP7 criteriaCode text:
#   RUNX1 (GN008) — "Intronic variants with SpliceAI ∆ scores <= 0.20."
#   MYOC  (GN019) — "Apply to intronic/noncoding ... variants if BP4 is met."
#   VHL   (GN078) — "BP7 can be applied to ... intronic variants where the
#                    PhyloP score is <=0.2." (still keeps its phyloP 0.2 gate)
# Genes that phrase a deep restriction differently (PIK3R1 "+1 to +6 and -1 to
# -20"; the LGMD panel "outside the splice donor/acceptor regions designated in
# Walker") are NOT here — they remain the +7/-21 default.
_BP7_INTRONIC_NONCANONICAL = frozenset({
    "RUNX1", "MYOC", "VHL",
})


# Per-gene, VARIANT-LEVEL BS1 exclusions: specific protein changes the VCEP bars
# from BS1 despite their population frequency, because they are the recurrent
# disease allele (a founder/common pathogenic variant whose frequency must NOT
# be read as benign evidence). Stored as the bare protein change; the BS1
# evaluator withholds BS1 (and the gnomAD-frequency path) for a matching variant.
#   MYOC (GN019) — "Does not apply to p.Gln368Ter" (the common pathogenic POAG
#   allele, ~2.6% allelic contribution; its frequency is disease-driven, not
#   benign). Not present in the JSON-LD BS1 description, so curated here.
_BS1_EXCLUDE_VARIANTS: dict[str, str] = {
    "MYOC": "p.Gln368Ter",
}


# Curated threshold corrections for genes whose cspec JSON-LD carries a typo'd or
# outdated AF threshold that disagrees with the VCEP's published classifications.
# Applied like a manual --override but persisted here so a rebuild keeps them.
# Keys are TSV column names; values overwrite the resolved row verbatim.
#   RPGR (GN106, X-linked IRD) — the cspec lists BS1 ">=8.3x10^-5" and a legacy
#   ACMG "5%" BA1 boilerplate, but the VCEP's published RPGR classifications use
#   BS1 > 5x10^-6 and BA1 > 5x10^-5 (BA1 = 10 x BS1), on the male hemizygous AF.
#   Both cspec numbers are typos; gnomAD male (XY) AF basis is already correct.
_CURATED_OVERRIDES: dict[str, dict[str, str]] = {
    "RPGR": {"bs1_threshold": "0.000005", "ba1_threshold": "0.00005"},
}


def _bs2_applicability(rule_set: dict) -> str:
    """"applicable" / "not_applicable" / "" for the rule set's BS2 code.

    not_applicable when the VCEP declines BS2, bars population data, or scores it
    purely on clinical phenotype/points with no gnomAD-countable rule — none of
    which a gnomAD-based evaluator can assess. "" when no BS2 code exists."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BS2":
            continue
        applic = _applicable_strengths(code)
        if not applic:
            return "not_applicable"
        desc = " ".join(es.get("description", "") or "" for es in applic)
        if _BS2_DECLINE.search(desc):
            return "not_applicable"
        # Pure point/phenotype scoring with no gnomAD-countable observation rule.
        if _BS2_POINTS.search(desc) and not _BS2_GNOMAD_COUNTABLE.search(desc):
            return "not_applicable"
        return "applicable"
    return ""


# BS2 minimum observation count — the LOWEST count at which the VCEP fires BS2 at
# any strength. Our evaluator is binary (met / not-met), so the operative bar is
# the gene's Supporting threshold, not its Strong one: a tiered spec (GUCY2D /
# AIPL1: Supporting ">=3 homozygotes in gnomAD", Strong ">=6") must fire BS2 at
# 3, and taking the MAX (6) caused a FALSE NEGATIVE at nhomalt=3. Strict cancer
# panels that demand many observations (CDH1, TP53) are handled separately —
# they resolve to not_applicable via _BS2_CLINICAL_CONFIRMATION, so the lower
# count here cannot leak a false-benign there. Only operator-anchored integers
# tied to an observation noun are taken (never "allele" counts or "20x
# coverage"), the MIN is kept, and the value is sanity-bounded. ">N" means N+1.
_BS2_COUNT = re.compile(
    r"(≥|≧|>=|&ge;|>|at least)\s*(\d{1,3})\b"
    r"(?!\s*[:/])"                                  # exclude ratios ("≥ 40:1")
    r"(?!\s*(?:year|yr|y\b|yo\b|month|week|day|%|x\b|×))"  # exclude ages/units (">18 years")
    r"(?=[^.]{0,30}(?:homozyg|hemizyg|heterozyg|individual|carrier|observ|"
    r"unrelated|proband|male|female|adult|case))",
    re.IGNORECASE,
)
_BS2_COUNT_MAX = 200  # reject mis-parsed large numbers (allele counts, ages)


def _pm4_applicability(rule_set: dict) -> str:
    """"applicable" / "not_applicable" / "" for the rule set's PM4 code.

    PM4 (protein-length change from in-frame indel / stop-loss) is declined by
    several VCEPs — notably cancer panels (BRCA1/2, the MMR genes, TP53, APC,
    PALB2) and the PI3K-pathway specs — where firing PM4 on an in-frame indel
    would wrongly add pathogenic weight. "" when the spec carries no PM4 code."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PM4":
            continue
        return "applicable" if _applicable_strengths(code) else "not_applicable"
    return ""


# Size-based PM4 downgrade: many VCEPs apply PM4 at Supporting (not the default
# Moderate) for a small in-frame indel — "single/one/1 amino acid" (→ <=1 aa) or
# "< 3 amino acid residues" (→ <=2 aa). Read from the PM4 Supporting strength.
_PM4_LT3_AA = re.compile(
    r"(?:less than|fewer than|under|<)\s*(?:3|three)\s+amino\s*acid", re.IGNORECASE
)
_PM4_SINGLE_AA = re.compile(r"(?:single|one|1)\s+amino\s*acid", re.IGNORECASE)


def _pm4_supporting_max_aa(rule_set: dict) -> str:
    """"1" / "2" / "" — the in-frame-indel size (in amino acids) at or below
    which PM4 downgrades to Supporting for the rule set's gene(s). Read from the
    PM4 Supporting strength description; "" when PM4 has no size-based Supporting
    tier (PM4 then fires at its default Moderate for any in-frame indel)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PM4":
            continue
        for es in _applicable_strengths(code):
            if es.get("label") != "Supporting":
                continue
            desc = es.get("description", "") or ""
            if _PM4_LT3_AA.search(desc):
                return "2"
            if _PM4_SINGLE_AA.search(desc):
                return "1"
        return ""
    return ""


def _bs2_count(rule_set: dict) -> str:
    """The VCEP's minimum BS2 observation count for the gene, or "" (use the
    global default). The LOWEST (min) operator-anchored count across all
    applicable strengths wins — that is the Supporting-level bar at which a
    binary BS2 first fires (see the comment above)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BS2":
            continue
        strs = _applicable_strengths(code)
        if not strs:
            return ""
        joined = " ".join(es.get("description", "") for es in strs)
        vals: list[int] = []
        for op, num in _BS2_COUNT.findall(joined):
            n = int(num) + (1 if op == ">" else 0)
            if 1 <= n <= _BS2_COUNT_MAX:
                vals.append(n)
        return str(min(vals)) if vals else ""
    return ""


# Count→strength tiers (e.g. GUCY2D "Strong >=6 homozygotes / Supporting >=3";
# BMPR2 "Strong >=3 / Moderate >=2 / Supporting >=1"). Each Applicable BS2
# strength carries its own operator-anchored count in its own description; the
# MIN count per strength is taken (the bar at which that strength first fires).
# Emitted only when >=2 distinct strengths carry a count (a real tiering) — a
# single tier is already captured by bs2_count (which fires Strong by default).
_BS2_STRENGTH_ORDER = ("Very Strong", "Strong", "Moderate", "Supporting")


def _bs2_tier_count(desc: str) -> Optional[int]:
    """The gnomAD-anchored BS2 count for one strength's description.

    A single strength may cite two thresholds — a phenotyped-literature count and
    a (usually higher) gnomAD count (GUCY2D Strong: ">=3 homozytes [literature] …
    or >=6 homozygotes in gnomAD"). The app scores BS2 on gnomAD, so a count
    immediately followed by "gnomAD" wins; only when none is gnomAD-anchored does
    the MIN of all counts apply (the bar at which the strength first fires)."""
    all_vals: list[int] = []
    gnomad_vals: list[int] = []
    for m in _BS2_COUNT.finditer(desc):
        n = int(m.group(2)) + (1 if m.group(1) == ">" else 0)
        if not (1 <= n <= _BS2_COUNT_MAX):
            continue
        all_vals.append(n)
        if "gnomad" in desc[m.start():m.end() + 40].lower():
            gnomad_vals.append(n)
    pool = gnomad_vals or all_vals
    return min(pool) if pool else None


def _bs2_strength(rule_set: dict) -> str:
    """"Strong:N,Moderate:M,Supporting:K" — the BS2 count→strength tiers for the
    rule set's gene(s), or "" when fewer than two strengths carry a count."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BS2":
            continue
        tiers: dict[str, int] = {}
        for es in _applicable_strengths(code):
            label = es.get("label", "")
            if label not in _BS2_STRENGTH_ORDER:
                continue
            n = _bs2_tier_count(es.get("description", "") or "")
            if n is not None:
                tiers[label] = n
        if len(tiers) >= 2:
            return ",".join(
                f"{s.replace(' ', '')}:{tiers[s]}"
                for s in _BS2_STRENGTH_ORDER if s in tiers
            )
        return ""
    return ""


# A VCEP whose BS2 counts only females (TP53: ">=8 unrelated females ... without
# cancer"; DICER1: "40+ unrelated females ... tumor-free"). Detected when the
# applicable BS2 descriptions mention "female"/"women" but NOT "male"/"men" —
# OTC ("(female) homozygotes or (male) hemizygotes") names both, so it stays a
# standard mode-based count, not female-only.
_BS2_FEMALE = re.compile(r"\bfemales?\b|\bwom[ae]n\b", re.IGNORECASE)
_BS2_MALE = re.compile(r"\bmales?\b|\bmen\b", re.IGNORECASE)


def _bs2_female_only(rule_set: dict) -> str:
    """"1" if the rule set's BS2 counts only females, else "" (count all sexes)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BS2":
            continue
        joined = " ".join(
            es.get("description", "") for es in _applicable_strengths(code)
        )
        if _BS2_FEMALE.search(joined) and not _BS2_MALE.search(joined):
            return "1"
        return ""
    return ""


# A dominant gene with incomplete penetrance whose BS2 is scored on HOMOZYGOTES,
# not heterozygous carriers (BMPR2/GN125, PIK3R2: ">=3 homozygotes in gnomAD
# controls"). Detected when the applicable BS2 descriptions mention homozygotes
# but NOT heterozygotes/carriers: a healthy het of an incompletely-penetrant
# dominant gene is not benign evidence, so the evaluator must count homozygotes
# (the default AD path counts carriers and would FALSELY fire BS2). Recessive /
# X-linked genes already count homozygotes/hemizygotes, so the flag is a no-op
# there — it only redirects the dominant (AD) path.
_BS2_HOM = re.compile(r"homozyg", re.IGNORECASE)
_BS2_HET = re.compile(r"heterozyg|\bcarrier", re.IGNORECASE)


def _bs2_hom_only(rule_set: dict) -> str:
    """"1" if the rule set's BS2 counts homozygotes only, else ""."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BS2":
            continue
        joined = " ".join(
            es.get("description", "") for es in _applicable_strengths(code)
        )
        if _BS2_HOM.search(joined) and not _BS2_HET.search(joined):
            return "1"
        return ""
    return ""


# PS1 splice handling, classified per gene into three states:
#   ""            -- no splice extension; PS1 is missense-only (original ACMG;
#                    e.g. GAA, the HCM genes). A splice variant must NOT get PS1.
#   "canonical"   -- splice extension covering canonical sites too. This is the
#                    DEFAULT for any extension: the ClinGen SVI splicing
#                    framework (Walker 2023, PMID 37352859) applies PS1 at
#                    canonical ±1/±2 sites in conjunction with PVS1 (RPE65, DYSF,
#                    IDUA, ACADVL, ADA, RPGR, HNF1A, RS1, …).
#   "noncanonical"-- splice extension EXPLICITLY limited to non-canonical
#                    positions (InSiGHT MMR / DICER1 "non-canonical splice"; BMPR2
#                    "outside the splice donor/acceptor +/-1,2") with no canonical
#                    handling.
# Caveat / exclusion sentences mentioning splicing (the standard "beware of
# changes that impact splicing", "should be excluded", comparison-variant notes)
# are stripped first so they are not mistaken for an extension.
_PS1_CAVEAT_SENT = re.compile(
    r"[^.]*\b(?:beware|should be excluded|should not be used|should not be provided"
    r"|not a predicted or confirmed splice defect"
    r"|rather than (?:at )?the amino acid|truly is a splice defect"
    r"|investigated[^.]*splic|predictions? \(by spliceai\)|spliceai scores? for both)"
    r"[^.]*\.?",
    re.IGNORECASE,
)
_PS1_SPLICE_MENTION = re.compile(r"splic|intronic", re.IGNORECASE)
# Explicit restriction to non-canonical positions only.
_PS1_NONCANON_ONLY = re.compile(r"non[- ]?canonical (?:splice|intronic)", re.IGNORECASE)
_PS1_OUTSIDE_ONLY = re.compile(
    r"outside[^.]{0,30}splice (?:donor|acceptor)[^.]{0,25}1,?\s*2", re.IGNORECASE
)
# Evidence that canonical ±1/±2 sites ARE handled (so the rule is not
# non-canonical-only): explicit "canonical splice", PS1 used with PVS1, or PS1
# applied AT/WITHIN the splice donor/acceptor ±1,2 positions.
_PS1_CANON_HANDLING = re.compile(
    r"(?<!non-)(?<!non )canonical[^.]{0,40}splic"
    r"|in conjunction with pvs1"
    r"|(?:located )?(?:at|within)[^.]{0,60}splice (?:donor|acceptor)[^.]{0,40}1,?\s*2",
    re.IGNORECASE,
)


def _ps1_applicability(rule_set: dict) -> str:
    """"applicable" / "not_applicable" / "" for the rule set's PS1 code. A VCEP
    that carries a PS1 code with no applicable strength declined it (e.g. CDH1)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PS1":
            continue
        return "applicable" if _applicable_strengths(code) else "not_applicable"
    return ""


def _pvs1_applicability(rule_set: dict) -> str:
    """"applicable" / "not_applicable" / "" for the rule set's PVS1 code. A VCEP
    that carries a PVS1 code with NO applicable strength has declined it for the
    gene — typically because loss-of-function is not the disease mechanism
    (gain-of-function / dominant-negative genes: MYOC, the RASopathy and
    cardiomyopathy panels, the activating PIK3 genes, RYR1, VWF, …). PVS1 must
    then never fire on a null variant in that gene."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PVS1":
            continue
        return "applicable" if _applicable_strengths(code) else "not_applicable"
    return ""


def _ps1_splice(rule_set: dict) -> str:
    """PS1 splice-extension state for the gene: "" (missense-only) / "canonical"
    / "noncanonical" — see the module-level note above."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PS1":
            continue
        descs = " ".join(es.get("description", "") or "" for es in _applicable_strengths(code))
        stripped = _PS1_CAVEAT_SENT.sub(" ", descs)
        if not _PS1_SPLICE_MENTION.search(stripped):
            return ""
        restricted = _PS1_NONCANON_ONLY.search(stripped) or _PS1_OUTSIDE_ONLY.search(stripped)
        if restricted and not _PS1_CANON_HANDLING.search(stripped):
            return "noncanonical"
        return "canonical"
    return ""


def _ps1_max(rule_set: dict) -> str:
    """Per-gene PS1 strength ceiling, emitted only when the VCEP caps PS1 BELOW
    its Strong default: "Supporting" (RMRP — "Downgraded to PS1_Supporting") or
    "Moderate". "" when the highest applicable PS1 strength is Strong or higher
    (no cap) or PS1 is not applicable. The PS1 counterpart of ``_pm5_max``."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PS1":
            continue
        ranks = [_STRENGTH_RANK[es["label"]] for es in _applicable_strengths(code)
                 if es.get("label") in _STRENGTH_RANK]
        if not ranks:
            return ""
        top = max(ranks)
        if top <= _STRENGTH_RANK["Supporting"]:
            return "Supporting"
        if top == _STRENGTH_RANK["Moderate"]:
            return "Moderate"
        return ""  # Strong (or higher) applicable → no cap (PS1 default is Strong)
    return ""


# BP1 ("missense in a gene where another mechanism dominates") gene-level
# applicability and TARGET consequence. Most VCEPs decline BP1; those that apply
# it target either missense (PALB2, APC, BRCA1/2) or — for gain-of-function
# RASopathy genes where loss-of-function is benign — TRUNCATING variants.
_BP1_TRUNCATING = re.compile(r"truncating", re.IGNORECASE)
_BP1_GOF = re.compile(r"gain[- ]of[- ]function", re.IGNORECASE)
# Some VCEPs leave a BP1/BP3 strength flagged "Applicable" but state in the text
# that the rule does not apply for the gene (KCNQ1, ITGA2B, ITGB3: missense also
# causes disease). The free-text decision overrides the flag.
_BP_TEXT_NA = re.compile(r"not applicable|does not apply", re.IGNORECASE)
# A positive BP3 clause ("BP3 can be applied to … / applies to …"): when present
# a trailing "not applicable" is a scoped exception (the rest of the gene), not a
# gene-level negation (VHL).
_BP3_POSITIVE = re.compile(r"can be applied|is applicable|applies to|may be applied", re.IGNORECASE)


def _bp1_applicability(rule_set: dict) -> tuple[str, str]:
    """(status, target) for the rule set's BP1 code: status is
    "applicable"/"not_applicable"/""; target is "truncating" (RASopathy GoF) or
    "missense"/"broad" (default) when applicable, else ""."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BP1":
            continue
        applic = _applicable_strengths(code)
        if not applic:
            return "not_applicable", ""
        desc = " ".join(es.get("description", "") or "" for es in applic)
        if _BP_TEXT_NA.search(desc):
            return "not_applicable", ""
        if "silent" in desc.lower():
            target = "broad"  # BRCA1/2: silent + missense + in-frame
        elif _BP1_TRUNCATING.search(desc) and _BP1_GOF.search(desc):
            target = "truncating"  # gain-of-function RASopathy genes
        else:
            target = "missense"
        return "applicable", target
    return "", ""


def _bp3_applicability(rule_set: dict) -> str:
    """"applicable"/"not_applicable"/"" for the rule set's BP3 code (in-frame
    indel in a repetitive region). A VCEP that declined BP3 → not_applicable."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BP3":
            continue
        applic = _applicable_strengths(code)
        if not applic:
            return "not_applicable"
        desc = " ".join(es.get("description", "") or "" for es in applic)
        # A "not applicable" sub-clause does NOT negate the gene when a positive
        # "BP3 can be applied to … <region>" clause is also present (VHL: applies
        # to the GXEEX repeat AA14-48, "not applicable" only for the rest).
        if _BP_TEXT_NA.search(desc) and not _BP3_POSITIVE.search(desc):
            return "not_applicable"
        return "applicable"
    return ""


# Residue ranges from a criterion description (BP1 excluded domains, BP3 allowed
# regions). Strips thousands separators, SpliceAI decimals, figure/PMID refs, and
# "N-M bp" indel-size ranges so only residue ranges remain.
_BP_NUM_RANGE = re.compile(r"(?<![\d.])(\d{1,4})\s*[-–]\s*(\d{1,4})(?!\d)")
_BP_AA = "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val"
_BP_AA_RANGE = re.compile(
    rf"(?:p\.)?(?:{_BP_AA})(\d{{1,4}})\s*[-–]\s*(?:p\.)?(?:{_BP_AA})(\d{{1,4}})", re.IGNORECASE
)
# "AA<n>-(AA)<m>" notation (VHL "AA14-AA48"). _BP_NUM_RANGE misses these because
# the dash is followed by "AA48", not a bare digit.
_BP_AA_PREFIX_RANGE = re.compile(r"AA\s?(\d{1,4})\s*[-–]\s*(?:AA\s?)?(\d{1,4})", re.IGNORECASE)


def _bp_residue_ranges(text: str) -> str:
    """";"-joined residue ranges parsed from *text* (e.g. "1021-1035",
    "2-101;1391-1424"). Empty string when none are found."""
    t = re.sub(r"(\d),(\d{3})\b", r"\1\2", text)
    t = re.sub(r"\d+\s*[-–]\s*\d+\s*bp|\d+\s*bp", " ", t, flags=re.IGNORECASE)  # indel sizes
    t = re.sub(r"\d+\.\d+|figure\s*\d+\w*|pmid:?\s*\d+|[≤<>]=?\s*[\d.]+", " ", t, flags=re.IGNORECASE)
    ranges: set[tuple[int, int]] = set()
    for m in (list(_BP_AA_PREFIX_RANGE.finditer(t))
              + list(_BP_AA_RANGE.finditer(t))
              + list(_BP_NUM_RANGE.finditer(t))):
        a, b = int(m.group(1)), int(m.group(2))
        if 0 < a < b <= 9999 and b > 4 and b - a < 4000:
            ranges.add((a, b))
    return ";".join(f"{a}-{b}" for a, b in sorted(ranges))


def _bp1_exclude(rule_set: dict) -> str:
    """Residue ranges where BP1 is NOT applied — APC's β-catenin repeat
    (1021-1035) and the BRCA1/2 clinically-important functional domains."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") == "BP1":
            return _bp_residue_ranges(
                " ".join(es.get("description", "") or "" for es in _applicable_strengths(code))
            )
    return ""


def _bp1_strength(rule_set: dict) -> str:
    """"Strong" when the VCEP applies BP1 at Strong (BRCA1/2); else "" (default
    Supporting)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") == "BP1":
            labels = {es.get("label") for es in _applicable_strengths(code)}
            return "Strong" if "Strong" in labels else ""
    return ""


_BP1_NOSPLICE = re.compile(r"no splic\w*\s+predicted|spliceai\s*[≤<]", re.IGNORECASE)


def _bp1_no_splice(rule_set: dict) -> str:
    """"yes" when BP1 requires no predicted splice impact (BRCA1/2: SpliceAI <=0.1)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") == "BP1":
            desc = " ".join(es.get("description", "") or "" for es in _applicable_strengths(code))
            return "yes" if _BP1_NOSPLICE.search(desc) else ""
    return ""


def _bp3_regions(rule_set: dict) -> str:
    """Residue ranges BP3 is RESTRICTED to (RPGR ORF15 585-1078; FOXG1 poly-AA
    tracts). Empty when BP3 is the generic repeat rule (any repeat region)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") == "BP3":
            return _bp_residue_ranges(
                " ".join(es.get("description", "") or "" for es in _applicable_strengths(code))
            )
    return ""


# --- BP7 per-gene phyloP "not highly conserved" cutoff -----------------------
# Walker 2023 BP7 requires the nucleotide to be NOT highly conserved. Most VCEPs
# state a phyloP cutoff X for this, phrased either as the BP7-eligible side
# ("not highly conserved = phyloP < X") or the conserved side ("conservation =
# phyloP > X"); either way the boundary value X is the cutoff the evaluator uses
# (block BP7 when phyloP >= X). All HUHVar phyloP is phyloP100way, matching the
# specs that say "phyloP100way < 2.0" (== the global default, a no-op). Real
# per-gene values seen: 0.1 (neurodev / coagulation panels), 0.2 (VHL), 1.5 (the
# platelet GP genes), 0 (RPGR, accelerated-only). The SpliceAI score and the
# alignment "OR" clause (1 primate / 3 mammals) are intentionally NOT modelled
# here — see BP7Evaluator.
_BP7_PHYLOP_RE = re.compile(
    r"phylop[^.]{0,45}?(?:<=|>=|<|>|less than|greater than|=)\s*([0-9]+(?:\.[0-9]+)?)",
    re.IGNORECASE,
)

# Some VCEPs declare conservation NON-informative for BP7, so phyloP must NOT be
# required at all (TP53 / ABCD1 / GALT / the SCID T-and-B-cell-development genes,
# whose poor cross-vertebrate conservation makes phyloP uninformative). These
# resolve to the sentinel "na" — the evaluator then skips the conservation gate
# entirely rather than applying any cutoff. Catches the several phrasings seen:
# "conservation is not required/considered/informative", "no conservation
# requirement", "Conservation does not have to be considered", "No requirement
# to assess for nucleotide conservation".
_BP7_NO_CONS_RE = re.compile(
    r"conservation[^.]{0,60}?\bnot\b[^.]{0,45}?"
        r"(considered|informative|required|necessary|relevant|assess|evaluat|applied|used)|"
    r"conservation[^.]{0,40}?no longer[^.]{0,30}?"
        r"(requirement|required|necessary|considered|relevant)|"
    r"\bno\b[^.]{0,25}?conservation[^.]{0,25}?(requirement|required|necessary)|"
    r"no requirement[^.]{0,45}?conservation|"
    r"conservation[^.]{0,40}?does not (?:have to|need)|"
    r"(?:do not|don't|not)[^.]{0,30}?(?:use|consider|assess|require|apply)[^.]{0,30}?conservation",
    re.IGNORECASE,
)


def _bp7_phylop(rule_set: dict) -> str:
    """The VCEP's BP7 phyloP policy for a gene:

    * a numeric cutoff string ("0.1", "0", "2.0") — block BP7 when phyloP >= it;
    * ``"na"`` — conservation is declared non-informative, so no phyloP gate;
    * ``""`` — the spec states nothing, so the evaluator keeps its global default.
    """
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BP7":
            continue
        for es in _applicable_strengths(code):
            desc = (es.get("description") or "")
            desc = desc.replace("≤", "<=").replace("≥", ">=").replace("\\", "")
            m = _BP7_PHYLOP_RE.search(desc)
            if m:
                return m.group(1)
            if _BP7_NO_CONS_RE.search(desc):
                return "na"
    return ""


# BP7 intronic applicability range. The Walker 2023 default ("") admits only
# DEEP-intronic variants (donor >= +7, acceptor <= -21). Some VCEPs (the
# RASopathy and PIK3-pathway panels) instead state BP7 applies to "intronic
# positions (except canonical splice sites)" — i.e. anywhere beyond the canonical
# +/-1,2 dinucleotides, a much broader range. Those resolve to "noncanonical" and
# the evaluator relaxes its distance gate to |distance| >= 3 (still gated on a
# benign SpliceAI prediction).
_BP7_NONCANONICAL_RE = re.compile(r"except[^.]{0,30}?canonical splice", re.IGNORECASE)
# Parametric BP7 intronic range. Two phrasings seen: Cardiomyopathy "(-4 and +7
# outward)" and SLC6A8 "(beyond -4bp or +7 bp)". Captures the two signed offsets
# (optional "bp" suffix); the negative is the acceptor cutoff, the positive the
# donor.
_BP7_PARAMETRIC_RE = re.compile(
    r"\(\s*(?:beyond\s*)?([+-]?\d+)(?:\s*bp)?\s+(?:and|or)\s+"
    r"([+-]?\d+)(?:\s*bp)?\s*(?:outward)?\s*\)",
    re.IGNORECASE,
)


def _bp7_intronic(rule_set: dict) -> str:
    """BP7 intronic range mode for a gene: ``"noncanonical"`` when the VCEP
    admits any intronic position except the canonical +/-1,2 sites;
    ``"donor:N,acceptor:-M"`` when it states explicit outward offsets (e.g. the
    Cardiomyopathy panel's "-4 and +7 outward"); ``""`` for the Walker
    deep-intronic (+7/-21) default."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BP7":
            continue
        for es in _applicable_strengths(code):
            desc = (es.get("description") or "").replace("\\", "")
            if _BP7_NONCANONICAL_RE.search(desc):
                return "noncanonical"
            m = _BP7_PARAMETRIC_RE.search(desc)
            if m:
                a, b = int(m.group(1)), int(m.group(2))
                donor = max(a, b)       # the positive (downstream-of-donor) offset
                acceptor = min(a, b)    # the negative (upstream-of-acceptor) offset
                return f"donor:{donor},acceptor:{acceptor}"
    return ""


# Per-gene SpliceAI "no-impact" cutoff for the BP4/BP7 splice gate. Only the
# benign-side operator (≤ / <) is captured — a "predicts impact ≥0.2" clause is
# the pathogenic side and must NOT be read as the no-impact cutoff. Emitted only
# when the gene's cutoff differs from the Walker 0.10 default.
_SPLICEAI_CUTOFF_RE = re.compile(
    r"spliceai[^.]{0,40}?(?:≤|<=|&le;|&lt;|<)\s*(0?\.\d+|0)\b",
    re.IGNORECASE,
)


def _splice_cutoff(rule_set: dict, label: str) -> str:
    """The non-default SpliceAI no-impact cutoff stated in the gene's *label*
    (``BP4`` / ``BP7``) code, or "" when it uses the 0.10 default / states none."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != label:
            continue
        for es in _applicable_strengths(code):
            desc = (es.get("description") or "").replace("\\", "")
            m = _SPLICEAI_CUTOFF_RE.search(desc)
            if m:
                val = float(m.group(1))
                if 0 <= val <= 1 and abs(val - 0.10) > 1e-9:
                    return m.group(1)
        return ""
    return ""


# --- REVEL per-gene PP3/BP4 thresholds (from the PP3 / BP4 criteria codes) ---
# Many VCEPs state a gene-specific REVEL cutoff for PP3/BP4 in place of the
# genome-wide Pejaver 2022 defaults. Each *applicable* evidence strength
# (Supporting/Moderate/Strong/Very Strong) of the PP3 (pathogenic, REVEL>=X) or
# BP4 (benign, REVEL<=X) code carries the cutoff for that tier. All numbers are
# anchored on the word "REVEL" so co-mentioned SpliceAI/CADD/AlphaMissense
# cutoffs in the same description are not mistaken for the REVEL value.
_REVEL_STRENGTH_COL = {
    "Supporting": "supporting", "Moderate": "moderate",
    "Strong": "strong", "Very Strong": "strong",  # no separate PP3 VS tier
}
# A REVEL score is a decimal in (0,1) — the leading-dot form excludes citation
# years ("2016"), superscript refs ("14"), and PMIDs, which are integers.
_REVEL_NUM = r"([01]?\.[0-9]+)"
# Other predictor names co-mentioned in the same PP3/BP4 description. The REVEL
# window is truncated at the first of these so e.g. "SpliceAI ≥0.2" or "CADD
# ≥20" is never mistaken for the REVEL cutoff.
_OTHER_TOOLS = re.compile(
    r"SpliceAI|CADD|AlphaMissense|BayesDel|phyloP|phastCons|GERP|MaxEntScan|"
    r"PolyPhen|\bSIFT\b|MutPred|\bVEST|PrimateAI|MetaLR|MetaSVM|PROVEAN|"
    r"MutationTaster|\bM-CAP\b",
    re.IGNORECASE)
_REVEL_TOKEN = re.compile(r"REVEL", re.IGNORECASE)
_PP3_OP = re.compile(
    r"(?:≥|≧|&ge;|&gt;=|\\?&gt;|>=|>|greater than or equal to|greater than|"
    r"at or above|above)\s*" + _REVEL_NUM, re.IGNORECASE)
_BP4_OP = re.compile(
    r"(?:≤|≦|&le;|&lt;=|\\?&lt;|<=|<|less than or equal to|less than|"
    r"at or below|below)\s*" + _REVEL_NUM, re.IGNORECASE)
_REVEL_RANGE = re.compile(_REVEL_NUM + r"\s*[-–]\s*" + _REVEL_NUM)
_REVEL_BETWEEN = re.compile(r"between\s+" + _REVEL_NUM + r"\s+and\s+" + _REVEL_NUM, re.IGNORECASE)
# Operator-free cutoff phrasing ("Use 0.75 as a discriminatory cut-off value",
# "a threshold of 0.7"). Used only as a fallback after operator/range parsing.
_REVEL_CUTOFF = re.compile(
    r"(?:cut[- ]?off|cutoff|threshold|discriminat\w*|\buse\b|\bat\b|score of)"
    r"\D{0,25}?" + _REVEL_NUM, re.IGNORECASE)


def _revel_windows(desc: str):
    """Yield text windows that begin at each "REVEL" mention and stop at the
    next other-tool name, so cutoff numbers stay attributed to REVEL alone."""
    for m in _REVEL_TOKEN.finditer(desc):
        w = desc[m.start(): m.start() + 200]
        other = _OTHER_TOOLS.search(w, m.end() - m.start())
        if other:
            w = w[: other.start()]
        yield w


def _revel_value(desc: str, pathogenic: bool) -> Optional[float]:
    """The REVEL firing cutoff in a description for the given direction.

    For a tier stated as a range/band the *firing edge* is taken: the lower
    bound for PP3 (REVEL >= edge) and the upper bound for BP4 (REVEL <= edge).
    Otherwise the directional operator value (">=" for PP3, "<=" for BP4) is
    used, and as a last resort an operator-free "cut-off" number. Numbers are
    only read from REVEL-attributed windows. Returns None when no REVEL value
    in (0, 1) is present."""
    op_re = _PP3_OP if pathogenic else _BP4_OP
    windows = list(_revel_windows(desc))
    for w in windows:
        for rng in (_REVEL_BETWEEN.search(w), _REVEL_RANGE.search(w)):
            if rng:
                vals = [v for v in (float(rng.group(1)), float(rng.group(2))) if 0.0 < v < 1.0]
                if vals:
                    return min(vals) if pathogenic else max(vals)
        op = op_re.search(w)
        if op:
            v = float(op.group(1))
            if 0.0 < v < 1.0:
                return v
    # Fallback: operator-free cutoff phrasing (e.g. FBN1 "Use 0.75 as a cut-off").
    for w in windows:
        cut = _REVEL_CUTOFF.search(w)
        if cut:
            v = float(cut.group(1))
            if 0.0 < v < 1.0:
                return v
    return None


def _revel_tiers(rule_set: dict, label: str) -> dict[str, float]:
    """{"supporting"/"moderate"/"strong": cutoff} parsed from the gene's PP3 or
    BP4 code, per applicable strength tier. Empty when the code is absent, not
    applicable, or carries no numeric REVEL cutoff (e.g. SCN1A's "follow ClinGen
    recommendations" prose → the evaluator then keeps the Pejaver defaults).

    A monotonicity guard drops any tier that contradicts its neighbours
    (PP3 cutoffs must rise Supporting<=Moderate<=Strong; BP4 cutoffs must fall),
    which rejects data-entry errors such as GN208's "0.774 - 0.092" Moderate."""
    pathogenic = label == "PP3"
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != label:
            continue
        out: dict[str, float] = {}
        for es in _applicable_strengths(code):
            col = _REVEL_STRENGTH_COL.get(es.get("label", ""))
            desc = es.get("description", "") or ""
            if not col or "REVEL" not in desc.upper():
                continue
            v = _revel_value(desc, pathogenic)
            if v is not None:
                # A stronger tier overrides a weaker one if both map to "strong"
                # (Strong + Very Strong); keep the more extreme firing edge.
                if col in out:
                    out[col] = min(out[col], v) if pathogenic else max(out[col], v)
                else:
                    out[col] = v
        return _revel_monotonic(out, pathogenic)
    return {}


def _revel_monotonic(tiers: dict[str, float], pathogenic: bool) -> dict[str, float]:
    """Drop tiers that break the expected ordering (guards parse/typo errors)."""
    order = ["supporting", "moderate", "strong"]
    present = [(t, tiers[t]) for t in order if t in tiers]
    kept: dict[str, float] = {}
    prev: Optional[float] = None
    for name, val in present:
        if prev is not None:
            # PP3: each stronger tier's cutoff must be >= the previous; BP4: <=.
            if (pathogenic and val < prev) or (not pathogenic and val > prev):
                continue  # contradictory tier — skip it
        kept[name] = val
        prev = val
    return kept


def _af_basis(rule_set: dict) -> str:
    """"males" if the applicable BA1/BS1 descriptions define the cutoff on the
    male (XY/hemizygous) allele frequency; otherwise "" (overall population)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") not in ("BA1", "BS1"):
            continue
        for es in _applicable_strengths(code):
            if _MALES_BASIS.search(es.get("description", "")):
                return "males"
    return ""


# PM2 per-gene extraction. PM2 strength is Supporting for almost every VCEP;
# only a handful set Moderate (GAA, LDLR, ETHE1, PDHA1, POLG, SLC19A3, ITGA2B,
# ITGB3). The cutoff is gene-specific and several VCEPs compare against the
# GrpMax Filtering Allele Frequency (FAF) rather than the raw popmax AF.
_PM2_FAF = re.compile(r"filtering allele frequency|grpmax\s+faf|\bFAF\b", re.IGNORECASE)
_PM2_ABSENT = re.compile(r"\babsent\b", re.IGNORECASE)
_PM2_AR = re.compile(r"recessive", re.IGNORECASE)
_PM2_AD = re.compile(r"dominant", re.IGNORECASE)
# PM2 cutoffs are stated with a LESS-THAN operator ("<", "≤", "less than",
# "below"). The shared _OP_NUM only matches GREATER-than (it was built for
# BA1/BS1), so PM2 needs its own. Using ONLY operator-anchored numbers — never
# the last-resort "any number" tier — is essential: PM2 descriptions are often
# long narratives ("prevalence 1/500", "accounts for ≤2% of variants") whose
# stray numbers would otherwise be mistaken for the cutoff.
_PM2_LT = re.compile(
    r"(?:≤|≦|<=|&le;|<|≲|less than or equal to|less than|below|at or below)\s*"
    r"([0-9]*\.?[0-9]+)\s*(%?)",
    re.IGNORECASE,
)
# Operator-led fraction cutoffs only ("≤ 1/300,000", "< 1:333,000",
# "≤ One out of 100,000"). A leading operator is required so a bare prevalence
# fraction in narrative ("(1/500 or lower)") is NOT picked up.
_PM2_FRACTION = re.compile(
    r"(?:≤|≦|<=|<|less than or equal to|less than)\s*"
    r"(?:1\s*[/:]\s*|one out of\s*)([0-9][0-9,]{2,})",
    re.IGNORECASE,
)
# PM2 is a rarity criterion; a credible cutoff is well below this. The bound
# rejects mis-parsed narrative numbers (e.g. "2%" → 0.02, "1/500" → 0.002 when
# operator-led). The few VCEPs that set a higher PM2 cutoff (cardiomyopathy
# 0.02, IDUA 0.025) state it only in an external table, so they fall back to
# the evaluator's global default — safe (stricter, never over-fires PM2).
_PM2_MAX_CREDIBLE = 0.01


def _pm2_cands(desc: str) -> list[float]:
    """PM2 frequency candidates (0, _PM2_MAX_CREDIBLE), operator-anchored only."""
    out: list[float] = []
    for m in _SCI.finditer(desc):
        v = float(m.group(1)) * 10 ** int(m.group(2))
        out.append(v)
    for n, p in _PM2_LT.findall(desc):
        v = _to_prop(n, p)
        if v is not None:
            out.append(v)
    for m in _PM2_FRACTION.finditer(desc):
        try:
            out.append(1.0 / int(m.group(1).replace(",", "")))
        except (ValueError, ZeroDivisionError):
            pass
    return sorted({v for v in out if 0.0 < v < _PM2_MAX_CREDIBLE})


def _pm2_code(rule_set: dict) -> Optional[dict]:
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") == "PM2":
            return code
    return None


def _pm2_strength(rule_set: dict) -> str:
    """VCEP PM2 strength: "Moderate" when the spec sets PM2 at Moderate (or
    higher), else "Supporting" when PM2 is applicable at Supporting, else ""
    (PM2 not applicable / no code). Aggregated across specs in main()."""
    code = _pm2_code(rule_set)
    if code is None:
        return ""
    ranks = [
        _STRENGTH_RANK[es["label"]]
        for es in _applicable_strengths(code)
        if es.get("label") in _STRENGTH_RANK
    ]
    if not ranks:
        return ""
    return "Moderate" if max(ranks) >= _STRENGTH_RANK["Moderate"] else "Supporting"


def _pm2_basis(rule_set: dict) -> str:
    """"faf" when the VCEP states the PM2 cutoff on the GrpMax Filtering Allele
    Frequency (HNF1A/HNF4A/GCK/SLC6A8/FBN1/LDLR…); else "" (raw popmax AF)."""
    code = _pm2_code(rule_set)
    if code is None:
        return ""
    desc = " ".join(es.get("description", "") or "" for es in _applicable_strengths(code))
    return "faf" if _PM2_FAF.search(desc) else ""


# PM2 subpopulation metric mode — the FAF95 (95% CI LOWER bound) is deflated for
# low-AC variants, so a variant elevated in a single subpopulation can slip under
# it. Two VCEP families correct this with a stricter highest-subpopulation metric:
#  * "point"  — RUNX1: "use GrpMax FAF when available, else require ALL
#               subpopulations meet the threshold." Evaluator additionally
#               requires the GrpMax POINT AF (popmax_af) <= threshold.
#  * "ci95"   — Cardiomyopathy/HCM: "<= threshold in the highest subpopulation
#               using the UPPER bound of the 95% CI." gnomAD shows only the FAF
#               (lower bound), so the evaluator reconstructs the upper bound from
#               the GrpMax AC/AN.
_PM2_SUBPOP_CAVEAT = re.compile(
    r"all sub-?populations?[^.]{0,30}?(?:meet|threshold)|"
    r"if a grpmax faf[^.]{0,40}?not\s+available",
    re.IGNORECASE,
)
_PM2_CI_UPPER = re.compile(r"upper bound of\s+(?:the\s+)?95\s*%?\s*(?:ci|confidence)", re.IGNORECASE)


def _pm2_subpop(rule_set: dict) -> str:
    """PM2 subpopulation mode: "ci95" (HCM upper-95%-CI), "point" (RUNX1
    all-subpopulations point AF), or "" (no special subpopulation rule)."""
    code = _pm2_code(rule_set)
    if code is None:
        return ""
    desc = " ".join(
        es.get("description", "") or "" for es in _applicable_strengths(code)
    ).replace("\\", "")
    if _PM2_CI_UPPER.search(desc):
        return "ci95"
    if _PM2_SUBPOP_CAVEAT.search(desc):
        return "point"
    return ""


# PM2 homozygote/hemizygote ceiling — several VCEPs additionally require few/no
# homozygotes or hemizygotes in gnomAD before PM2 applies, because a recurrent
# homozygous/hemizygous observation argues against rarity-for-pathogenicity even
# when the allele frequency is low. Encoded "<scope>:<max>" where scope is hom /
# hemi / homhemi (combined) and max is the highest tolerated count:
#   SLC6A8 "0 homo- or hemizygotes"  -> homhemi:0
#   OTC    "<=1 homo- or hemizygote" -> homhemi:1
#   ADA/DCLRE1C/IL7R/JAK3/RAG1/RAG2 "no homozygotes" -> hom:0
#   GATM/GAMT "Homozygotes should not be seen"        -> hom:0
#   ABCD1  "absent in hemizygotes"   -> hemi:0
_PM2_HOMHEMI = re.compile(
    r"(no|zero|0|≤\s*1|<=\s*1|one|≤\s*0|1)\s*homo-?\s*(?:or|/)\s*hemizygot",
    re.IGNORECASE,
)
_PM2_HEMI_ONLY = re.compile(
    r"absent[^.]{0,25}hemizygot|hemizygot[^.]{0,15}(?:absent|0\b|zero|not)",
    re.IGNORECASE,
)
_PM2_HOM_ONLY = re.compile(
    r"(?:no|zero|0|not)[^.]{0,30}homozygot|"
    r"homozygot[^.]{0,30}(?:not be seen|not been observed|should not)",
    re.IGNORECASE,
)


def _pm2_zygosity(rule_set: dict) -> str:
    """PM2 homozygote/hemizygote ceiling as "<scope>:<max>" (hom/hemi/homhemi),
    or "" when the VCEP states no such requirement."""
    code = _pm2_code(rule_set)
    if code is None:
        return ""
    desc = " ".join(
        es.get("description", "") or "" for es in _applicable_strengths(code)
    ).replace("\\", "")
    m = _PM2_HOMHEMI.search(desc)
    if m:
        tok = m.group(1).lower().replace(" ", "")
        return "homhemi:1" if ("1" in tok or "one" in tok) else "homhemi:0"
    if _PM2_HEMI_ONLY.search(desc):
        return "hemi:0"
    if _PM2_HOM_ONLY.search(desc):
        return "hom:0"
    return ""


def _pm2_threshold(rule_set: dict, moi: str) -> str:
    """PM2 allele-frequency cutoff for a gene with mode-of-inheritance *moi*.

    Returns "" when no PM2 cutoff can be resolved (the evaluator then keeps its
    global default). "0" means the VCEP requires the variant to be ABSENT (only
    a truly absent/AC=0 variant qualifies). For MOI-split descriptions
    ("≤X for recessive, ≤Y for dominant") the gene's own MOI selects the value.
    """
    code = _pm2_code(rule_set)
    if code is None:
        return ""
    strs = _applicable_strengths(code)
    if not strs:
        return ""
    desc = strs[0].get("description", "")
    cands = _pm2_cands(desc)
    if not cands:
        # The first tier has no numeric cutoff. Before falling back to
        # "absent"=0, scan the OTHER applicable tiers for a gnomAD numeric
        # cutoff: some VCEPs state a legacy "Absent from ESP/1000G/ExAC"
        # boilerplate at one strength and the operative gnomAD number at another
        # (ITGA2B/ITGB3 GN011: Moderate "absent" + Supporting "<0.0001 in
        # gnomAD"). The numeric cutoff governs whether PM2 fires at all, so it
        # must win over the boilerplate "absent". (strs[0]-with-a-number is left
        # untouched, so no other gene's resolved threshold changes.)
        rest = " ".join(es.get("description", "") or "" for es in strs[1:])
        rest_cands = _pm2_cands(rest)
        if rest_cands:
            desc, cands = rest, rest_cands
        else:
            return "0" if _PM2_ABSENT.search(desc) else ""
    # MOI-split: the recessive cutoff is the higher value, dominant the lower.
    if len(cands) > 1 and _PM2_AR.search(desc) and _PM2_AD.search(desc):
        if "AR" in moi and "AD" not in moi:
            return _fmt(max(cands))
        if "AD" in moi or "XL" in moi:
            return _fmt(min(cands))
        return _fmt(min(cands))  # unknown MOI → stricter (FP-minimising) value
    # Single operative cutoff (specs state the operative number first).
    return _fmt(cands[0])


def _moi(gene: dict) -> str:
    out: set[str] = set()
    for dis in gene.get("diseases", []):
        for m in dis.get("modeOfInheritance", []):
            lbl = str(m.get("@label", "")).lower()
            if "recessive" in lbl:
                out.add("AR")
            elif "dominant" in lbl:  # includes semidominant
                out.add("AD")
            elif "x-linked" in lbl:
                out.add("XL")
    order = {"AD": 0, "AR": 1, "XL": 2}
    return ",".join(sorted(out, key=lambda x: order.get(x, 9)))


def _ui_url(doc_id: str, gn: str) -> str:
    return f"https://cspec.genome.network/cspec/ui/svi/svi/{gn}"


def parse_spec(path: str) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    gn = os.path.splitext(os.path.basename(path))[0]
    vcep = (d.get("affiliation") or {}).get("label") or d.get("label", gn)
    url = _ui_url(d.get("@id", ""), gn)
    status = d.get("cspecStatus", "")
    # Spec specificity: how many distinct genes this whole spec covers. A
    # single-gene VCEP (e.g. GN035 FOXG1) is more gene-specific — and therefore
    # authoritative for that gene — than a grouped panel (e.g. GN016, 6 genes).
    # Used to resolve a gene appearing in multiple specs (see main()).
    n_spec_genes = len({
        g.get("label") for rs in d.get("ruleSets", [])
        for g in rs.get("genes", []) if g.get("label")
    })
    rows: list[dict] = []
    for rs in d.get("ruleSets", []):
        # A rule set with no criteriaCodes is a Pilot/In-Prep placeholder: it
        # carries no thresholds or criteria, so it must never shadow a populated
        # spec's values even when it is more gene-specific (see the base-row
        # resolution in main()).
        rs_empty = len(rs.get("criteriaCodes", [])) == 0
        bs1, bs1_note = _criterion_threshold(rs, "BS1")
        ba1, ba1_note = _criterion_threshold(rs, "BA1")
        bs1_strength = _bs1_strength(rs)
        af_basis = _af_basis(rs)
        pm2_strength = _pm2_strength(rs)
        pm2_basis = _pm2_basis(rs)
        pm2_subpop = _pm2_subpop(rs)
        pm2_zygosity = _pm2_zygosity(rs)
        pm4 = _pm4_applicability(rs)
        pm4_smax = _pm4_supporting_max_aa(rs)
        pp2_map = _pp2_applicability(rs)
        pp2_req = _pp2_requires(rs)
        pm5_op = _pm5_grantham_op(rs)
        pm5_excl = _pm5_excludes(rs)
        pm5_max = _pm5_max(rs)
        pm5_lp = _pm5_lp_comparator(rs)
        bs2 = _bs2_applicability(rs)
        bs2_count = _bs2_count(rs)
        bs2_strength = _bs2_strength(rs)
        bs2_female_only = _bs2_female_only(rs)
        bs2_hom_only = _bs2_hom_only(rs)
        ps1 = _ps1_applicability(rs)
        ps1_splice = _ps1_splice(rs)
        ps1_max = _ps1_max(rs)
        pvs1 = _pvs1_applicability(rs)
        bp1, bp1_target = _bp1_applicability(rs)
        bp1_exclude = _bp1_exclude(rs)
        bp1_strength = _bp1_strength(rs)
        bp1_no_splice = _bp1_no_splice(rs)
        bp3 = _bp3_applicability(rs)
        bp3_regions = _bp3_regions(rs)
        bp7_phylop = _bp7_phylop(rs)
        bp7_intronic = _bp7_intronic(rs)
        bp4_splice_cutoff = _splice_cutoff(rs, "BP4")
        bp7_splice_cutoff = _splice_cutoff(rs, "BP7")
        revel_pp3 = _revel_tiers(rs, "PP3")
        revel_bp4 = _revel_tiers(rs, "BP4")
        notes = "; ".join(n for n in (ba1_note, bs1_note) if n and "not applicable" not in n and "absent" not in n)
        for gene in rs.get("genes", []):
            sym = gene.get("label")
            if not sym:
                continue
            moi = _moi(gene)
            pm2_threshold = _pm2_threshold(rs, moi)
            rows.append({
                "gene_symbol": sym,
                "inheritance": moi,
                "prevalence": "", "allelic_het": "", "genetic_het": "", "penetrance": "",
                "bs1_threshold": "" if bs1 is None else _fmt(bs1),
                "bs1_strength": bs1_strength,
                "bs1_exclude": "",
                "ba1_threshold": "" if ba1 is None else _fmt(ba1),
                "af_basis": af_basis,
                "pm2_threshold": "",
                "pm2_strength": "",
                "pm2_basis": "",
                "pm2_subpop": "",
                "pm2_zygosity": "",
                "pm4": "",
                "pm4_supporting_max_aa": "",
                "pp2": "",
                "pp2_requires": "",
                "pm5_grantham": "",
                "pm5_excludes": "",
                "pm5_max": "",
                "pm5_lp": "",
                "bs2": "",
                "bs2_count": "",
                "bs2_strength": "",
                "bs2_female_only": "",
                "bs2_hom_only": "",
                "pvs1": "",
                "ps1": "",
                "ps1_splice": "",
                "ps1_max": "",
                "ps1_paralog_group": "",
                "ps1_paralog_strength": "",
                "bp1": "",
                "bp1_target": "",
                "bp1_exclude": "",
                "bp1_strength": "",
                "bp1_no_splice": "",
                "bp3": "",
                "bp3_regions": "",
                "bp7_phylop": "",
                "bp7_intronic": "",
                "bp4_splice_cutoff": "",
                "bp7_splice_cutoff": "",
                "revel_pp3_supporting": "",
                "revel_pp3_moderate": "",
                "revel_pp3_strong": "",
                "revel_bp4_supporting": "",
                "revel_bp4_moderate": "",
                "revel_bp4_strong": "",
                "source_vcep": vcep,
                "cspec_url": url,
                "notes": "; ".join(x for x in (f"{gn} {status}", notes) if x),
                # Transient (not TSV columns; dropped at write time via
                # extrasaction="ignore"): drive multi-spec resolution and the
                # cross-spec PP2/PM5 aggregation in main().
                "_gn": gn,
                "_specificity": n_spec_genes,
                "_empty": rs_empty,
                "_pm2_threshold": pm2_threshold,
                "_pm2_strength": pm2_strength,
                "_pm2_basis": pm2_basis,
                "_pm2_subpop": pm2_subpop,
                "_pm2_zygosity": pm2_zygosity,
                "_pm4": pm4,
                "_pm4_smax": pm4_smax,
                "_pp2": pp2_map.get(sym, ""),
                "_pp2_requires": pp2_req,
                "_pm5_grantham": pm5_op,
                "_pm5_excludes": pm5_excl,
                "_pm5_max": pm5_max,
                "_pm5_lp": pm5_lp,
                "_bs2": bs2,
                "_bs2_count": bs2_count,
                "_bs2_strength": bs2_strength,
                "_bs2_female_only": bs2_female_only,
                "_bs2_hom_only": bs2_hom_only,
                "_ps1": ps1,
                "_ps1_splice": ps1_splice,
                "_ps1_max": ps1_max,
                "_pvs1": pvs1,
                "_bp1": bp1,
                "_bp1_target": bp1_target,
                "_bp1_exclude": bp1_exclude,
                "_bp1_strength": bp1_strength,
                "_bp1_no_splice": bp1_no_splice,
                "_bp3": bp3,
                "_bp3_regions": bp3_regions,
                "_bp7_phylop": bp7_phylop,
                "_bp7_intronic": bp7_intronic,
                "_bp4_splice_cutoff": bp4_splice_cutoff,
                "_bp7_splice_cutoff": bp7_splice_cutoff,
                "_revel_pp3": revel_pp3,
                "_revel_bp4": revel_bp4,
            })
    return rows


def _fmt(x: float) -> str:
    # Compact but exact-ish representation (avoids 0.001 -> 1e-3 surprises).
    return ("%.10f" % x).rstrip("0").rstrip(".") if x < 1 else str(x)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-dir", default="resources/clingen")
    ap.add_argument("--out", default="resources/shared/disease_prevalence.tsv")
    ap.add_argument("--released-only", action="store_true",
                    help="Only emit specs whose cspecStatus is 'Released'.")
    ap.add_argument(
        "--override", action="append", default=[],
        metavar="GENE:field=val[,field=val]",
        help="Manually override a gene's resolved values, applied AFTER "
             "multi-spec resolution. Repeatable. Fields: ba1, bs1, af_basis, "
             "inheritance. Example: --override RYR1:ba1=0.0038,bs1=0.0007 . "
             "Use this to pin a disease-appropriate cutoff for genes whose "
             "multi-spec ambiguity is kept conservative (highest threshold) by "
             "default.",
    )
    args = ap.parse_args()
    overrides = _parse_overrides(args.override)

    files = sorted(glob.glob(os.path.join(args.json_dir, "GN*.json")))
    by_gene: dict[str, dict] = {}
    # PP2 applicability resolves to the MOST GENE-SPECIFIC spec's decision, the
    # same specificity rule used for BA1/BS1: a single-gene VCEP supersedes a
    # grouped panel. This is essential — the grouped RASopathy spec (GN004, 12
    # genes) declares PP2 "applicable to all", but the single-gene RASopathy
    # specs (NRAS/GN039, SOS1/GN041, …) explicitly decline it; the gene-specific
    # "not_applicable" must win. On a specificity tie the conservative
    # (FP-minimising) "not_applicable" wins. pp2_choice[g] = (specificity,
    # decision, requires); only specs that carry an actual PP2 decision count.
    pp2_choice: dict[str, tuple[int, str, str]] = {}
    # PM5 Grantham gate, aggregated across specs: inclusive "ge" outranks strict
    # "gt" (the spec's primary, less-strict rule wins on a cross-spec conflict).
    # Inclusive gates outrank strict ones within an engine; the only BLOSUM gene
    # (PTEN) is a single-gene VCEP so cross-engine ties do not arise in practice.
    pm5_rank = {"ge": 2, "gt": 1, "blosum_le": 2, "blosum_lt": 1, "": 0}
    pm5_by_gene: dict[str, str] = {}
    # PM5 exclusions (union of criteria across specs) and strength ceiling
    # ("Moderate" outranks "Supporting": a gene any VCEP allows at Moderate is
    # not capped to Supporting).
    pm5_excludes_by_gene: dict[str, set[str]] = {}
    pm5_max_rank = {"Strong": 3, "Moderate": 2, "Supporting": 1, "": 0}
    pm5_max_by_gene: dict[str, str] = {}
    # PM5 comparator-significance policy, resolved to the most gene-specific spec
    # (a single-gene VCEP supersedes a grouped panel), like PP2. Stores
    # (specificity, "yes"/"no"); a "no" gene requires a Pathogenic comparator.
    pm5_lp_choice: dict[str, tuple[int, str]] = {}
    # BS2 applicability, resolved to the most gene-specific spec (single-gene
    # VCEP over a grouped panel), like PP2; on a tie the conservative
    # not_applicable wins. bs2_choice[g] = (specificity, decision).
    bs2_choice: dict[str, tuple[int, str, str]] = {}
    # PS1 splice "non-canonical only" restriction (union across specs).
    # PS1 splice-extension state, resolved to the most gene-specific spec that
    # carries a PS1 code; on a tie an explicit extension beats "" (missense-only).
    ps1_splice_choice: dict[str, tuple[int, str]] = {}
    # PS1 applicability, resolved to the most gene-specific spec (like PP2/BS2).
    ps1_choice: dict[str, tuple[int, str, str]] = {}
    # PVS1 applicability (gain-of-function / dominant-negative genes decline it),
    # most gene-specific spec; on a tie the conservative not_applicable wins.
    pvs1_choice: dict[str, tuple[int, str, str]] = {}
    # BP1 / BP3 applicability, most-specific spec (like PP2); BP1 also carries the
    # target consequence (missense / truncating).
    bp1_choice: dict[str, tuple[int, str, str]] = {}
    bp1_fields_by_gene: dict[str, dict] = {}
    bp3_choice: dict[str, tuple[int, str, str]] = {}
    bp3_regions_by_gene: dict[str, str] = {}
    # BP7 per-gene phyloP cutoff, resolved to the most gene-specific spec that
    # states a phyloP number (like PM2/REVEL). bp7_phylop_choice[g] =
    # (specificity, cutoff_string).
    bp7_phylop_choice: dict[str, tuple[int, str]] = {}
    # BP7 intronic range mode ("noncanonical"), most gene-specific spec wins.
    bp7_intronic_choice: dict[str, tuple[int, str]] = {}
    # Per-gene SpliceAI no-impact cutoff for BP4 / BP7, most gene-specific wins.
    bp4_splice_choice: dict[str, tuple[int, str]] = {}
    bp7_splice_choice: dict[str, tuple[int, str]] = {}
    # PM2 threshold/strength/basis, resolved to the most gene-specific spec that
    # carries an applicable PM2 code (like PP2/BS2/PS1). pm2_choice[g] =
    # (specificity, threshold, strength, basis).
    pm2_choice: dict[str, tuple[int, str, str, str]] = {}
    # PM4 applicability, most-specific spec (like PP2/BS2); on a tie the
    # conservative not_applicable wins.
    pm4_choice: dict[str, tuple[int, str, str]] = {}
    # REVEL PP3/BP4 per-gene cutoffs, resolved to the most gene-specific spec
    # that states a numeric REVEL threshold (like PM2). revel_choice[g] =
    # (specificity, pp3_tiers, bp4_tiers). On a specificity tie the first spec
    # (file order) is kept; genuinely conflicting same-specificity specs (e.g.
    # RYR1 across GN012/GN150) need a manual --override.
    revel_choice: dict[str, tuple[int, dict, dict]] = {}
    n_specs = 0
    for path in files:
        try:
            with open(path, encoding="utf-8") as fh:
                status = json.load(fh).get("cspecStatus", "")
        except Exception as e:
            print(f"  [skip] {path}: {e}")
            continue
        if args.released_only and status != "Released":
            continue
        n_specs += 1
        for row in parse_spec(path):
            g = row["gene_symbol"]
            # Spec-preference override: a handful of genes have empty single-gene
            # specs (In-Prep/Pilot, 0 criteriaCodes) that would otherwise win the
            # "most gene-specific spec" rule and shadow a fully-populated grouped
            # spec. For those genes we PIN the authoritative spec so it wins every
            # criterion (and the base row). E.g. ITGA2B/ITGB3 → GN011 (Platelet
            # Disorders, Released, 28 codes) over the empty GN059/GN060/GN221-223.
            if g in _FORCE_SPEC and row["_gn"] != _FORCE_SPEC[g]:
                continue
            if row["_pp2"] in ("applicable", "not_applicable"):
                cand = (row["_specificity"], row["_pp2"], row["_pp2_requires"])
                if g not in pp2_choice or _pp2_more_specific(cand, pp2_choice[g]):
                    pp2_choice[g] = cand
            if pm5_rank[row["_pm5_grantham"]] > pm5_rank[pm5_by_gene.get(g, "")]:
                pm5_by_gene[g] = row["_pm5_grantham"]
            if row["_pm5_excludes"]:
                pm5_excludes_by_gene.setdefault(g, set()).update(
                    row["_pm5_excludes"].split(",")
                )
            if pm5_max_rank[row["_pm5_max"]] > pm5_max_rank[pm5_max_by_gene.get(g, "")]:
                pm5_max_by_gene[g] = row["_pm5_max"]
            if row["_pm5_lp"] in ("yes", "no"):
                spec = row["_specificity"]
                cur = pm5_lp_choice.get(g)
                # Most gene-specific spec wins; on a tie the conservative
                # "no" (Pathogenic comparator required) wins.
                if (cur is None or spec < cur[0]
                        or (spec == cur[0] and row["_pm5_lp"] == "no" and cur[1] != "no")):
                    pm5_lp_choice[g] = (spec, row["_pm5_lp"])
            if row["_pm4"] in ("applicable", "not_applicable"):
                cand = (row["_specificity"], row["_pm4"], row["_pm4_smax"])
                if g not in pm4_choice or _pp2_more_specific(cand, pm4_choice[g]):
                    pm4_choice[g] = cand
            if row["_bs2"] in ("applicable", "not_applicable"):
                # Carry bs2_count (3rd slot), female-only (4th), hom-only (5th)
                # and the count→strength tiers (6th) so all come from the same
                # (most-specific) spec that decided applicability.
                cand = (row["_specificity"], row["_bs2"], row["_bs2_count"],
                        row["_bs2_female_only"], row["_bs2_hom_only"],
                        row["_bs2_strength"])
                if g not in bs2_choice or _pp2_more_specific(cand, bs2_choice[g]):
                    bs2_choice[g] = cand
            # PM2: only specs with an applicable PM2 code contribute. Most
            # gene-specific spec wins; on a specificity tie keep the first
            # (deterministic by file order).
            if row["_pm2_strength"]:
                spec = row["_specificity"]
                cur = pm2_choice.get(g)
                if cur is None or spec < cur[0]:
                    pm2_choice[g] = (
                        spec, row["_pm2_threshold"], row["_pm2_strength"],
                        row["_pm2_basis"], row["_pm2_subpop"], row["_pm2_zygosity"],
                    )
            # REVEL: only specs that state a numeric REVEL cutoff contribute.
            # Most gene-specific spec wins; on a tie keep the first (file order).
            if row["_revel_pp3"] or row["_revel_bp4"]:
                spec = row["_specificity"]
                cur = revel_choice.get(g)
                if cur is None or spec < cur[0]:
                    revel_choice[g] = (spec, row["_revel_pp3"], row["_revel_bp4"])
            if row["_bp1"] in ("applicable", "not_applicable"):
                cand = (row["_specificity"], row["_bp1"], "")
                if g not in bp1_choice or _pp2_more_specific(cand, bp1_choice[g]):
                    bp1_choice[g] = cand
                    bp1_fields_by_gene[g] = (
                        {
                            "bp1_target": row["_bp1_target"],
                            "bp1_exclude": row["_bp1_exclude"],
                            "bp1_strength": row["_bp1_strength"],
                            "bp1_no_splice": row["_bp1_no_splice"],
                        }
                        if row["_bp1"] == "applicable"
                        else {}
                    )
            if row["_bp3"] in ("applicable", "not_applicable"):
                cand = (row["_specificity"], row["_bp3"], "")
                if g not in bp3_choice or _pp2_more_specific(cand, bp3_choice[g]):
                    bp3_choice[g] = cand
                    bp3_regions_by_gene[g] = (
                        row["_bp3_regions"] if row["_bp3"] == "applicable" else ""
                    )
            # BP7 phyloP cutoff: only specs that state a phyloP number contribute;
            # most gene-specific spec wins, ties keep the first (file order).
            if row["_bp7_phylop"]:
                spec = row["_specificity"]
                cur = bp7_phylop_choice.get(g)
                if cur is None or spec < cur[0]:
                    bp7_phylop_choice[g] = (spec, row["_bp7_phylop"])
            # BP7 intronic range mode (noncanonical), most gene-specific spec wins.
            if row["_bp7_intronic"]:
                spec = row["_specificity"]
                cur = bp7_intronic_choice.get(g)
                if cur is None or spec < cur[0]:
                    bp7_intronic_choice[g] = (spec, row["_bp7_intronic"])
            # Per-gene BP4/BP7 SpliceAI cutoff, most gene-specific spec wins.
            for col, choice in (
                ("_bp4_splice_cutoff", bp4_splice_choice),
                ("_bp7_splice_cutoff", bp7_splice_choice),
            ):
                if row[col]:
                    spec = row["_specificity"]
                    cur = choice.get(g)
                    if cur is None or spec < cur[0]:
                        choice[g] = (spec, row[col])
            if row["_ps1"] in ("applicable", "not_applicable"):
                cand = (row["_specificity"], row["_ps1"], row["_ps1_max"])
                if g not in ps1_choice or _pp2_more_specific(cand, ps1_choice[g]):
                    ps1_choice[g] = cand
                # Splice mode comes from a spec that actually carries a PS1 code.
                spec, mode = row["_specificity"], row["_ps1_splice"]
                cur = ps1_splice_choice.get(g)
                if (cur is None or spec < cur[0]
                        or (spec == cur[0] and mode and not cur[1])):
                    ps1_splice_choice[g] = (spec, mode)
            if row["_pvs1"] in ("applicable", "not_applicable"):
                cand = (row["_specificity"], row["_pvs1"], "")
                if g not in pvs1_choice or _pp2_more_specific(cand, pvs1_choice[g]):
                    pvs1_choice[g] = cand
            if g not in by_gene:
                by_gene[g] = row
                continue
            # Duplicate gene across specs. Prefer the more gene-specific spec
            # (fewer genes in scope): a single-gene VCEP supersedes a grouped
            # panel for that gene. On a specificity tie (e.g. RYR1/ACTA1 — the
            # same gene in distinct single-gene VCEPs for different diseases)
            # the ambiguity cannot be auto-resolved, so we default to the most
            # CONSERVATIVE cutoff — the highest BA1 (then highest BS1) — which
            # minimises false-positive benign calls. Use --override to pin a
            # disease-appropriate value instead.
            prev = by_gene[g]
            # An empty Pilot/In-Prep placeholder spec (0 criteriaCodes) must
            # never shadow a populated spec's thresholds, even when it is more
            # gene-specific. This generalises the per-gene _FORCE_SPEC pins: e.g.
            # the empty single-gene GN061 AKT3 / GN050 CDH23 specs no longer
            # blank out the Released grouped GN018 / GN005 BS1/BA1 values.
            if row["_empty"] != prev["_empty"]:
                if prev["_empty"]:  # incoming row is populated → it wins
                    row["notes"] = (row["notes"] + "; populated over empty placeholder").strip("; ")
                    by_gene[g] = row
                continue
            new_spec = row["_specificity"]
            old_spec = prev["_specificity"]
            if new_spec < old_spec:
                row["notes"] = (row["notes"] + "; multiple specs (gene-specific kept)").strip("; ")
                by_gene[g] = row
            elif new_spec > old_spec:
                prev["notes"] = (prev["notes"] + "; multiple specs").strip("; ")
            else:
                new_ba1 = _to_float(row["ba1_threshold"])
                old_ba1 = _to_float(prev["ba1_threshold"])
                new_bs1 = _to_float(row["bs1_threshold"])
                old_bs1 = _to_float(prev["bs1_threshold"])
                # Conservative (FP-minimising) tie-break: higher BA1 wins; if
                # BA1 ties, higher BS1 wins.
                new_key = (new_ba1 if new_ba1 is not None else -1.0,
                           new_bs1 if new_bs1 is not None else -1.0)
                old_key = (old_ba1 if old_ba1 is not None else -1.0,
                           old_bs1 if old_bs1 is not None else -1.0)
                if new_key > old_key:
                    row["notes"] = (row["notes"] + "; multiple specs (kept conservative)").strip("; ")
                    by_gene[g] = row
                else:
                    prev["notes"] = (prev["notes"] + "; multiple specs").strip("; ")

    # Stamp the cross-spec PP2 applicability onto each gene's resolved row
    # (before overrides, so a manual --override pp2=… still wins).
    for g, row in by_gene.items():
        choice = pp2_choice.get(g)
        row["pp2"] = choice[1] if choice else ""
        row["pp2_requires"] = choice[2] if choice and choice[1] == "applicable" else ""
        row["pm5_grantham"] = pm5_by_gene.get(g, "")
        # PM1[,PS1] union; canonical order. "Supporting" ceiling only when no
        # VCEP allowed Moderate (pm5_max_by_gene resolved to "Supporting").
        excl = pm5_excludes_by_gene.get(g, set())
        row["pm5_excludes"] = ",".join(c for c in ("PM1", "PS1") if c in excl)
        # Surface the actual ceiling (Supporting caps PM5; Strong unlocks the
        # PM5_Strong tier; Moderate/"" keep the default Moderate ceiling).
        row["pm5_max"] = pm5_max_by_gene.get(g, "")
        # Only "no" (Pathogenic comparator required) is recorded; "yes"/none
        # leave the column blank (LP comparator accepted — the default).
        lp = pm5_lp_choice.get(g)
        row["pm5_lp"] = "no" if lp and lp[1] == "no" else ""
        bchoice = bs2_choice.get(g)
        row["bs2"] = bchoice[1] if bchoice else ""
        # Genes whose VCEP BS2 needs phenotype/lab/functional confirmation or
        # internal-only clinical data that gnomAD lacks: force not_applicable so a
        # gnomAD-count BS2 is never (falsely) fired on a pathogenic variant. The
        # X-linked set additionally requires a phenotyped male/hemizygote cohort
        # gnomAD cannot supply (see _BS2_XLINKED_INTERNAL).
        if (g in _BS2_CLINICAL_CONFIRMATION or g in _BS2_XLINKED_INTERNAL) and row["bs2"]:
            row["bs2"] = "not_applicable"
        # Emit the per-gene BS2 count and female-only flag only when BS2 is
        # applicable for the gene.
        _bs2_applic = row["bs2"] == "applicable"
        row["bs2_count"] = bchoice[2] if (_bs2_applic and bchoice) else ""
        row["bs2_female_only"] = bchoice[3] if (_bs2_applic and bchoice) else ""
        row["bs2_hom_only"] = bchoice[4] if (_bs2_applic and bchoice) else ""
        row["bs2_strength"] = bchoice[5] if (_bs2_applic and bchoice) else ""
        pm4c = pm4_choice.get(g)
        row["pm4"] = pm4c[1] if pm4c else ""
        row["pm4_supporting_max_aa"] = (
            pm4c[2] if (pm4c and pm4c[1] == "applicable") else ""
        )
        pm2c = pm2_choice.get(g)
        # Emit pm2_strength only when Moderate (Supporting is the global default
        # the evaluator already applies); always emit the resolved threshold/basis.
        row["pm2_threshold"] = pm2c[1] if pm2c else ""
        row["pm2_strength"] = pm2c[2] if (pm2c and pm2c[2] == "Moderate") else ""
        row["pm2_basis"] = pm2c[3] if pm2c else ""
        row["pm2_subpop"] = pm2c[4] if pm2c else ""
        row["pm2_zygosity"] = pm2c[5] if pm2c else ""
        pchoice = ps1_choice.get(g)
        row["ps1"] = pchoice[1] if pchoice else ""
        row["ps1_splice"] = ps1_splice_choice.get(g, (0, ""))[1]
        # PS1 strength cap: only meaningful when PS1 is applicable for the gene.
        row["ps1_max"] = pchoice[2] if (pchoice and pchoice[1] == "applicable") else ""
        # PS1 paralogue / analogous-residue group (curated).
        row["ps1_paralog_group"], row["ps1_paralog_strength"] = _PS1_PARALOG.get(g, ("", ""))
        vchoice = pvs1_choice.get(g)
        row["pvs1"] = vchoice[1] if vchoice else ""
        bp1c = bp1_choice.get(g)
        row["bp1"] = bp1c[1] if bp1c else ""
        bp1f = bp1_fields_by_gene.get(g, {}) if bp1c and bp1c[1] == "applicable" else {}
        row["bp1_target"] = bp1f.get("bp1_target", "")
        row["bp1_exclude"] = bp1f.get("bp1_exclude", "")
        row["bp1_strength"] = bp1f.get("bp1_strength", "")
        row["bp1_no_splice"] = bp1f.get("bp1_no_splice", "")
        bp3c = bp3_choice.get(g)
        row["bp3"] = bp3c[1] if bp3c else ""
        row["bp3_regions"] = bp3_regions_by_gene.get(g, "") if bp3c and bp3c[1] == "applicable" else ""
        bp7c = bp7_phylop_choice.get(g)
        row["bp7_phylop"] = bp7c[1] if bp7c else ""
        # Curated correction for genes whose "conservation non-informative" BP7
        # policy is in the spec prose but absent from the JSON-LD (e.g. the LCA
        # genes, BRCA1/2).
        if g in _BP7_CONSERVATION_NA:
            row["bp7_phylop"] = "na"
        bp7ic = bp7_intronic_choice.get(g)
        row["bp7_intronic"] = bp7ic[1] if bp7ic else ""
        b4s = bp4_splice_choice.get(g)
        row["bp4_splice_cutoff"] = b4s[1] if b4s else ""
        b7s = bp7_splice_choice.get(g)
        row["bp7_splice_cutoff"] = b7s[1] if b7s else ""
        # Curated broad-intronic genes (RUNX1/MYOC/VHL) whose "any intronic
        # position" policy is not phrased as "except canonical splice sites".
        if g in _BP7_INTRONIC_NONCANONICAL:
            row["bp7_intronic"] = "noncanonical"
        # Curated variant-level BS1 exclusion (e.g. MYOC p.Gln368Ter).
        row["bs1_exclude"] = _BS1_EXCLUDE_VARIANTS.get(g, "")
        rvc = revel_choice.get(g)
        rv_pp3 = rvc[1] if rvc else {}
        rv_bp4 = rvc[2] if rvc else {}
        for tier in ("supporting", "moderate", "strong"):
            row[f"revel_pp3_{tier}"] = _fmt(rv_pp3[tier]) if tier in rv_pp3 else ""
            row[f"revel_bp4_{tier}"] = _fmt(rv_bp4[tier]) if tier in rv_bp4 else ""

    # Curated typo/outdated-threshold corrections first, then CLI --override (so
    # an explicit CLI override still wins over a curated default).
    _apply_overrides(by_gene, _CURATED_OVERRIDES)
    _apply_overrides(by_gene, overrides)
    rows = [by_gene[g] for g in sorted(by_gene)]
    # newline="\n" forces LF on every platform; csv otherwise emits CRLF on
    # Windows, which would rewrite every line of the committed LF-normalised TSV.
    with open(args.out, "w", newline="\n", encoding="utf-8") as fh:
        # extrasaction="ignore": rows carry a transient "_specificity" key that
        # is not a TSV column.
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    with_thr = sum(1 for r in rows if r["ba1_threshold"] or r["bs1_threshold"])
    print(f"specs parsed: {n_specs} | gene rows: {len(rows)} | with BA1/BS1: {with_thr}")
    print(f"written → {args.out}")


def _to_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# CLI override field → TSV column.
_OVERRIDE_FIELDS = {
    "ba1": "ba1_threshold",
    "bs1": "bs1_threshold",
    "bs1_strength": "bs1_strength",
    "bs1_exclude": "bs1_exclude",
    "af_basis": "af_basis",
    "pm2_threshold": "pm2_threshold",
    "pm2_strength": "pm2_strength",
    "pm2_basis": "pm2_basis",
    "pm2_subpop": "pm2_subpop",
    "pm2_zygosity": "pm2_zygosity",
    "pm4": "pm4",
    "pm4_supporting_max_aa": "pm4_supporting_max_aa",
    "pp2": "pp2",
    "pp2_requires": "pp2_requires",
    "pm5_grantham": "pm5_grantham",
    "pm5_excludes": "pm5_excludes",
    "pm5_max": "pm5_max",
    "pm5_lp": "pm5_lp",
    "bs2": "bs2",
    "bs2_count": "bs2_count",
    "bs2_strength": "bs2_strength",
    "bs2_female_only": "bs2_female_only",
    "bs2_hom_only": "bs2_hom_only",
    "pvs1": "pvs1",
    "ps1": "ps1",
    "ps1_splice": "ps1_splice",
    "ps1_max": "ps1_max",
    "ps1_paralog_group": "ps1_paralog_group",
    "ps1_paralog_strength": "ps1_paralog_strength",
    "bp1": "bp1",
    "bp1_target": "bp1_target",
    "bp1_exclude": "bp1_exclude",
    "bp1_strength": "bp1_strength",
    "bp1_no_splice": "bp1_no_splice",
    "bp3": "bp3",
    "bp3_regions": "bp3_regions",
    "bp7_phylop": "bp7_phylop",
    "bp7_intronic": "bp7_intronic",
    "bp4_splice_cutoff": "bp4_splice_cutoff",
    "bp7_splice_cutoff": "bp7_splice_cutoff",
    "revel_pp3_supporting": "revel_pp3_supporting",
    "revel_pp3_moderate": "revel_pp3_moderate",
    "revel_pp3_strong": "revel_pp3_strong",
    "revel_bp4_supporting": "revel_bp4_supporting",
    "revel_bp4_moderate": "revel_bp4_moderate",
    "revel_bp4_strong": "revel_bp4_strong",
    "inheritance": "inheritance",
}


def _parse_overrides(specs: list[str]) -> dict[str, dict[str, str]]:
    """Parse ``--override GENE:field=val[,field=val]`` strings into a
    {gene: {column: value}} map. Raises ValueError on a malformed spec or an
    unknown field so a typo fails loudly rather than being silently ignored."""
    out: dict[str, dict[str, str]] = {}
    for spec in specs:
        gene, sep, kvs = spec.partition(":")
        gene = gene.strip()
        if not sep or not gene:
            raise ValueError(f"--override '{spec}': expected GENE:field=value[,...]")
        fields: dict[str, str] = {}
        for kv in kvs.split(","):
            if not kv.strip():
                continue
            key, eq, val = kv.partition("=")
            key = key.strip().lower()
            if not eq or key not in _OVERRIDE_FIELDS:
                raise ValueError(
                    f"--override '{spec}': bad field '{kv.strip()}' "
                    f"(allowed: {', '.join(sorted(_OVERRIDE_FIELDS))})"
                )
            fields[_OVERRIDE_FIELDS[key]] = val.strip()
        out.setdefault(gene, {}).update(fields)
    return out


def _apply_overrides(by_gene: dict[str, dict], overrides: dict[str, dict[str, str]]) -> None:
    """Apply manual per-gene overrides in place (last word on a value).

    For multi-spec genes whose ambiguity cannot be auto-resolved (e.g. RYR1,
    ACTA1 — same specificity across distinct diseases), this lets the operator
    pin the disease-appropriate cutoff. A gene not produced by any spec is
    created so an override can also add a brand-new row."""
    for gene, fields in overrides.items():
        row = by_gene.get(gene)
        if row is None:
            row = {c: "" for c in COLUMNS}
            row["gene_symbol"] = gene
            by_gene[gene] = row
        row.update(fields)
        row["notes"] = (row.get("notes", "") + "; manual override").strip("; ")


if __name__ == "__main__":
    main()
