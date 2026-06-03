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
    "penetrance", "bs1_threshold", "ba1_threshold", "source_vcep", "cspec_url",
    "notes",
]

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
                "source_vcep": vcep,
                "cspec_url": url,
                "notes": "; ".join(x for x in (f"{gn} {status}", notes) if x),
                # Transient (not a TSV column; dropped at write time via
                # extrasaction="ignore"): drives multi-spec resolution.
                "_specificity": n_spec_genes,
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
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.json_dir, "GN*.json")))
    by_gene: dict[str, dict] = {}
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
            if g not in by_gene:
                by_gene[g] = row
                continue
            # Duplicate gene across specs. Prefer the more gene-specific spec
            # (fewer genes in scope): a single-gene VCEP supersedes a grouped
            # panel for that gene. Only on a specificity tie do we fall back to
            # keeping the most conservative (highest) BA1.
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
                if new_ba1 is not None and (old_ba1 is None or new_ba1 > old_ba1):
                    row["notes"] = (row["notes"] + "; multiple specs (kept conservative)").strip("; ")
                    by_gene[g] = row
                else:
                    prev["notes"] = (prev["notes"] + "; multiple specs").strip("; ")

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


if __name__ == "__main__":
    main()
