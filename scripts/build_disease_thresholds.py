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
    "penetrance", "bs1_threshold", "ba1_threshold", "af_basis", "pp2",
    "pp2_requires", "pm5_grantham", "pm5_excludes", "pm5_max", "pm5_lp",
    "source_vcep", "cspec_url", "notes",
]

# Genes whose BA1/BS1 spec defines the cutoff on the *male* (XY/hemizygous)
# allele frequency rather than the overall population FAF — detected from
# "in males" / "hemizygous" wording in the applicable BA1/BS1 descriptions
# (e.g. RPGR, RS1, ABCD1, SLC6A8, OTC, all X-linked). Emitted in the af_basis
# column so the BA1/BS1 evaluators compare against gnomAD AF_XY for these genes.
_MALES_BASIS = re.compile(r"\bin males\b|hemizyg", re.IGNORECASE)

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
    """Per-spec PM5 strength ceiling: "Supporting" (only Supporting applicable —
    ATM, CDH1, PALB2), "Moderate" (Moderate or higher applicable), or "" (no
    applicable PM5). Aggregated across specs in main(): "Moderate" outranks
    "Supporting" so a gene that any VCEP allows at Moderate keeps the default."""
    for code in rule_set.get("criteriaCodes", []):
        if code.get("label") != "PM5":
            continue
        ranks = [_STRENGTH_RANK[es["label"]] for es in _applicable_strengths(code)
                 if es.get("label") in _STRENGTH_RANK]
        if not ranks:
            return ""
        return "Supporting" if max(ranks) == _STRENGTH_RANK["Supporting"] else "Moderate"
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
        af_basis = _af_basis(rs)
        pp2_map = _pp2_applicability(rs)
        pp2_req = _pp2_requires(rs)
        pm5_op = _pm5_grantham_op(rs)
        pm5_excl = _pm5_excludes(rs)
        pm5_max = _pm5_max(rs)
        pm5_lp = _pm5_lp_comparator(rs)
        notes = "; ".join(n for n in (ba1_note, bs1_note) if n and "not applicable" not in n and "absent" not in n)
        for gene in rs.get("genes", []):
            sym = gene.get("label")
            if not sym:
                continue
            rows.append({
                "gene_symbol": sym,
                "inheritance": _moi(gene),
                "prevalence": "", "allelic_het": "", "genetic_het": "", "penetrance": "",
                "bs1_threshold": "" if bs1 is None else _fmt(bs1),
                "ba1_threshold": "" if ba1 is None else _fmt(ba1),
                "af_basis": af_basis,
                "pp2": "",
                "pp2_requires": "",
                "pm5_grantham": "",
                "pm5_excludes": "",
                "pm5_max": "",
                "pm5_lp": "",
                "source_vcep": vcep,
                "cspec_url": url,
                "notes": "; ".join(x for x in (f"{gn} {status}", notes) if x),
                # Transient (not TSV columns; dropped at write time via
                # extrasaction="ignore"): drive multi-spec resolution and the
                # cross-spec PP2/PM5 aggregation in main().
                "_specificity": n_spec_genes,
                "_pp2": pp2_map.get(sym, ""),
                "_pp2_requires": pp2_req,
                "_pm5_grantham": pm5_op,
                "_pm5_excludes": pm5_excl,
                "_pm5_max": pm5_max,
                "_pm5_lp": pm5_lp,
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
    pm5_max_rank = {"Moderate": 2, "Supporting": 1, "": 0}
    pm5_max_by_gene: dict[str, str] = {}
    # PM5 comparator-significance policy, resolved to the most gene-specific spec
    # (a single-gene VCEP supersedes a grouped panel), like PP2. Stores
    # (specificity, "yes"/"no"); a "no" gene requires a Pathogenic comparator.
    pm5_lp_choice: dict[str, tuple[int, str]] = {}
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
        row["pm5_max"] = "Supporting" if pm5_max_by_gene.get(g, "") == "Supporting" else ""
        # Only "no" (Pathogenic comparator required) is recorded; "yes"/none
        # leave the column blank (LP comparator accepted — the default).
        lp = pm5_lp_choice.get(g)
        row["pm5_lp"] = "no" if lp and lp[1] == "no" else ""

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
    "af_basis": "af_basis",
    "pp2": "pp2",
    "pp2_requires": "pp2_requires",
    "pm5_grantham": "pm5_grantham",
    "pm5_excludes": "pm5_excludes",
    "pm5_max": "pm5_max",
    "pm5_lp": "pm5_lp",
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
