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
    """
    def _pct(num: str, pct: str) -> Optional[float]:
        try:
            v = float(num)
        except ValueError:
            return None
        if pct:
            v /= 100.0
        return v if 0.0 < v < 1.0 else None

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

    if not cands:
        return None, False
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
            # Duplicate gene across specs — keep the most conservative BA1.
            prev = by_gene[g]
            new_ba1 = _to_float(row["ba1_threshold"])
            old_ba1 = _to_float(prev["ba1_threshold"])
            if new_ba1 is not None and (old_ba1 is None or new_ba1 > old_ba1):
                row["notes"] = (row["notes"] + "; multiple specs (kept conservative)").strip("; ")
                by_gene[g] = row
            else:
                prev["notes"] = (prev["notes"] + "; multiple specs").strip("; ")

    rows = [by_gene[g] for g in sorted(by_gene)]
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
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
