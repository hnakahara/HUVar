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
        --out resources/clingen/disease_prevalence.tsv

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
    "penetrance", "bs1_threshold", "bs1_strength", "ba1_threshold", "af_basis",
    "pm2_threshold", "pm2_strength", "pm2_basis", "pm4", "pp2",
    "pp2_requires", "pm5_grantham", "pm5_excludes", "pm5_max", "pm5_lp", "bs2",
    "bs2_count",
    "ps1", "ps1_splice", "bp1", "bp1_target", "bp1_exclude", "bp1_strength",
    "bp1_no_splice", "bp3", "bp3_regions",
    "source_vcep", "cspec_url", "notes",
]

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
    """"ge"/"gt"/"" — the PM5 Grantham gate operator for this rule set's gene(s).

    Returns "" when no applicable PM5 strength mentions Grantham. When the
    description carries an "equal" clause the gate is inclusive (``ge``); a
    strict "higher/less than" phrasing with no "equal" yields ``gt``. If both
    appear across strengths the inclusive ``ge`` wins (the spec's primary rule)."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PM5":
            continue
        op = ""
        for es in _applicable_strengths(code):
            desc = es.get("description", "") or ""
            if "grantham" not in desc.lower():
                continue
            if "equal" in desc.lower():
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


# BS2 (observed in a healthy individual) gene-level applicability. A VCEP may
# carry a BS2 code but decline general-population/gnomAD data for it (RASopathy
# GN004: "general population data should not be used for this criterion" — due
# to variable expressivity). Since our BS2 evaluator is gnomAD-based, treat that
# as not_applicable.
_BS2_NO_POPDATA = re.compile(
    r"population data should not be used|should not be used for this criterion",
    re.IGNORECASE,
)


def _bs2_applicability(rule_set: dict) -> str:
    """"applicable" / "not_applicable" / "" for the rule set's BS2 code.

    Applicable when the BS2 code has an Applicable strength whose description
    does not bar population data; not_applicable when the VCEP carries a BS2
    code it declined (or barred population data); "" when no BS2 code exists."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BS2":
            continue
        applic = _applicable_strengths(code)
        if not applic:
            return "not_applicable"
        if _BS2_NO_POPDATA.search(applic[0].get("description", "") or ""):
            return "not_applicable"
        return "applicable"
    return ""


# BS2 minimum observation count. Most VCEPs fire BS2 on 1-2 homozygotes (handled
# by the global default), but cancer/strict panels demand many more (CDH1 >=10
# individuals, TP53 >=8, PDHA1 >=16 hemizygotes) — under-counting there would
# FALSELY mark a pathogenic variant benign. Only operator-anchored integers tied
# to an observation noun are taken (never "allele" counts or "20x coverage"),
# the MAX is kept (the strictest gnomAD-applicable bar), and the value is
# sanity-bounded. ">N" means N+1 (">1 homozygote" → 2).
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


def _bs2_count(rule_set: dict) -> str:
    """The VCEP's minimum BS2 observation count for the gene, or "" (use the
    global default). The strictest (max) operator-anchored count wins."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "BS2":
            continue
        strs = _applicable_strengths(code)
        if not strs:
            return ""
        vals: list[int] = []
        for op, num in _BS2_COUNT.findall(strs[0].get("description", "")):
            n = int(num) + (1 if op == ">" else 0)
            if 1 <= n <= _BS2_COUNT_MAX:
                vals.append(n)
        return str(max(vals)) if vals else ""
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
        # No numeric cutoff: an "absent from controls" rule means threshold 0.
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
        bs1, bs1_note = _criterion_threshold(rs, "BS1")
        ba1, ba1_note = _criterion_threshold(rs, "BA1")
        bs1_strength = _bs1_strength(rs)
        af_basis = _af_basis(rs)
        pm2_strength = _pm2_strength(rs)
        pm2_basis = _pm2_basis(rs)
        pm4 = _pm4_applicability(rs)
        pp2_map = _pp2_applicability(rs)
        pp2_req = _pp2_requires(rs)
        pm5_op = _pm5_grantham_op(rs)
        pm5_excl = _pm5_excludes(rs)
        pm5_max = _pm5_max(rs)
        pm5_lp = _pm5_lp_comparator(rs)
        bs2 = _bs2_applicability(rs)
        bs2_count = _bs2_count(rs)
        ps1 = _ps1_applicability(rs)
        ps1_splice = _ps1_splice(rs)
        bp1, bp1_target = _bp1_applicability(rs)
        bp1_exclude = _bp1_exclude(rs)
        bp1_strength = _bp1_strength(rs)
        bp1_no_splice = _bp1_no_splice(rs)
        bp3 = _bp3_applicability(rs)
        bp3_regions = _bp3_regions(rs)
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
                "ba1_threshold": "" if ba1 is None else _fmt(ba1),
                "af_basis": af_basis,
                "pm2_threshold": "",
                "pm2_strength": "",
                "pm2_basis": "",
                "pm4": "",
                "pp2": "",
                "pp2_requires": "",
                "pm5_grantham": "",
                "pm5_excludes": "",
                "pm5_max": "",
                "pm5_lp": "",
                "bs2": "",
                "bs2_count": "",
                "ps1": "",
                "ps1_splice": "",
                "bp1": "",
                "bp1_target": "",
                "bp1_exclude": "",
                "bp1_strength": "",
                "bp1_no_splice": "",
                "bp3": "",
                "bp3_regions": "",
                "source_vcep": vcep,
                "cspec_url": url,
                "notes": "; ".join(x for x in (f"{gn} {status}", notes) if x),
                # Transient (not TSV columns; dropped at write time via
                # extrasaction="ignore"): drive multi-spec resolution and the
                # cross-spec PP2/PM5 aggregation in main().
                "_specificity": n_spec_genes,
                "_pm2_threshold": pm2_threshold,
                "_pm2_strength": pm2_strength,
                "_pm2_basis": pm2_basis,
                "_pm4": pm4,
                "_pp2": pp2_map.get(sym, ""),
                "_pp2_requires": pp2_req,
                "_pm5_grantham": pm5_op,
                "_pm5_excludes": pm5_excl,
                "_pm5_max": pm5_max,
                "_pm5_lp": pm5_lp,
                "_bs2": bs2,
                "_bs2_count": bs2_count,
                "_ps1": ps1,
                "_ps1_splice": ps1_splice,
                "_bp1": bp1,
                "_bp1_target": bp1_target,
                "_bp1_exclude": bp1_exclude,
                "_bp1_strength": bp1_strength,
                "_bp1_no_splice": bp1_no_splice,
                "_bp3": bp3,
                "_bp3_regions": bp3_regions,
            })
    return rows


def _fmt(x: float) -> str:
    # Compact but exact-ish representation (avoids 0.001 -> 1e-3 surprises).
    return ("%.10f" % x).rstrip("0").rstrip(".") if x < 1 else str(x)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-dir", default="resources/clingen")
    ap.add_argument("--out", default="resources/clingen/disease_prevalence.tsv")
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
    pm5_rank = {"ge": 2, "gt": 1, "": 0}
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
    # BP1 / BP3 applicability, most-specific spec (like PP2); BP1 also carries the
    # target consequence (missense / truncating).
    bp1_choice: dict[str, tuple[int, str, str]] = {}
    bp1_fields_by_gene: dict[str, dict] = {}
    bp3_choice: dict[str, tuple[int, str, str]] = {}
    bp3_regions_by_gene: dict[str, str] = {}
    # PM2 threshold/strength/basis, resolved to the most gene-specific spec that
    # carries an applicable PM2 code (like PP2/BS2/PS1). pm2_choice[g] =
    # (specificity, threshold, strength, basis).
    pm2_choice: dict[str, tuple[int, str, str, str]] = {}
    # PM4 applicability, most-specific spec (like PP2/BS2); on a tie the
    # conservative not_applicable wins.
    pm4_choice: dict[str, tuple[int, str, str]] = {}
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
                cand = (row["_specificity"], row["_pm4"], "")
                if g not in pm4_choice or _pp2_more_specific(cand, pm4_choice[g]):
                    pm4_choice[g] = cand
            if row["_bs2"] in ("applicable", "not_applicable"):
                # Carry bs2_count in the 3rd slot so the count comes from the
                # same (most-specific) spec that decided applicability.
                cand = (row["_specificity"], row["_bs2"], row["_bs2_count"])
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
                        spec, row["_pm2_threshold"], row["_pm2_strength"], row["_pm2_basis"],
                    )
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
            if row["_ps1"] in ("applicable", "not_applicable"):
                cand = (row["_specificity"], row["_ps1"], "")
                if g not in ps1_choice or _pp2_more_specific(cand, ps1_choice[g]):
                    ps1_choice[g] = cand
                # Splice mode comes from a spec that actually carries a PS1 code.
                spec, mode = row["_specificity"], row["_ps1_splice"]
                cur = ps1_splice_choice.get(g)
                if (cur is None or spec < cur[0]
                        or (spec == cur[0] and mode and not cur[1])):
                    ps1_splice_choice[g] = (spec, mode)
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
        # Emit the per-gene BS2 count only when BS2 is applicable for the gene.
        row["bs2_count"] = bchoice[2] if (bchoice and bchoice[1] == "applicable") else ""
        pm4c = pm4_choice.get(g)
        row["pm4"] = pm4c[1] if pm4c else ""
        pm2c = pm2_choice.get(g)
        # Emit pm2_strength only when Moderate (Supporting is the global default
        # the evaluator already applies); always emit the resolved threshold/basis.
        row["pm2_threshold"] = pm2c[1] if pm2c else ""
        row["pm2_strength"] = pm2c[2] if (pm2c and pm2c[2] == "Moderate") else ""
        row["pm2_basis"] = pm2c[3] if pm2c else ""
        pchoice = ps1_choice.get(g)
        row["ps1"] = pchoice[1] if pchoice else ""
        row["ps1_splice"] = ps1_splice_choice.get(g, (0, ""))[1]
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

    _apply_overrides(by_gene, overrides)
    rows = [by_gene[g] for g in sorted(by_gene)]
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        # extrasaction="ignore": rows carry a transient "_specificity" key that
        # is not a TSV column.
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
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
    "af_basis": "af_basis",
    "pm2_threshold": "pm2_threshold",
    "pm2_strength": "pm2_strength",
    "pm2_basis": "pm2_basis",
    "pm4": "pm4",
    "pp2": "pp2",
    "pp2_requires": "pp2_requires",
    "pm5_grantham": "pm5_grantham",
    "pm5_excludes": "pm5_excludes",
    "pm5_max": "pm5_max",
    "pm5_lp": "pm5_lp",
    "bs2": "bs2",
    "bs2_count": "bs2_count",
    "ps1": "ps1",
    "ps1_splice": "ps1_splice",
    "bp1": "bp1",
    "bp1_target": "bp1_target",
    "bp1_exclude": "bp1_exclude",
    "bp1_strength": "bp1_strength",
    "bp1_no_splice": "bp1_no_splice",
    "bp3": "bp3",
    "bp3_regions": "bp3_regions",
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
