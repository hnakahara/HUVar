#!/usr/bin/env python
"""Compute, from analysis_code_f1.tsv, the sensitivity/specificity/precision/F1 (+TP/FP/FN/TN)
of each criterion (ACMG code) and write them to a TSV.

- erepo_count = actual positives (variants to which the expert panel applied the code)
- tool_count  = variants to which the tool applied the code
- TP = recall*erepo_count (= precision*tool_count)
- FP = tool_count - TP,  FN = erepo_count - TP,  TN = N - TP - FP - FN  (N = number of evaluated variants)
- Specificity = TN/(TN+FP)
- codes not implemented (tool_count 0 in both builds) are N/A for that method (not written)
"""
import csv
from pathlib import Path

CODE_SRC = Path("outputs/analysis_code_f1.tsv")
CAT_SRC = Path("outputs/analysis_category_metrics.tsv")
OUT = Path("outputs/analysis_criteria_metrics.tsv")

# tool name in analysis_code_f1.tsv -> display name (6 methods, in order)
METHOD_ORDER = [
    ("InterVar",        "InterVar"),
    ("BIAS-2015",       "BIAS-2015 (default)"),
    ("BIAS_noPP5BP6",   "BIAS-2015 (no PP5/BP6)"),
    ("HUHVar_default",  "HUVar ESM1b+OpenSpliceAI"),
    ("HUHVar_spliceai", "HUVar AlphaMissense+SpliceAI"),
    ("HUHVar_revel",    "HUVar REVEL/BayesDel/CADD+SpliceAI"),
]
CODE_ORDER = ["PVS1", "PS1", "PS2", "PS3", "PS4",
              "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
              "PP1", "PP2", "PP3", "PP4", "PP5",
              "BA1", "BS1", "BS2", "BS3", "BS4",
              "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]
BUILD_ORDER = ["hg19", "hg38"]


def to_int(s):
    s = (s or "").strip()
    return int(s) if s else 0


def to_float(s):
    s = (s or "").strip()
    return float(s) if s else None


# --- N (number of evaluated variants) per build ---
N = {}
with CAT_SRC.open(encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        N[r["build"]] = int(r["n_eval"])

# --- read code F1: (build, tool, code) -> row ---
raw = {}
tool_emits = {}  # method name -> set of codes (implemented if tool_count > 0 at least once)
with CODE_SRC.open(encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        raw[(r["build"], r["tool"], r["code"])] = r
        if to_int(r["tool_count"]) > 0:
            tool_emits.setdefault(r["tool"], set()).add(r["code"])


def fmt(x):
    return "" if x is None else f"{x:.4f}"


rows_out = []
for build in BUILD_ORDER:
    n = N[build]
    for src_method, disp in METHOD_ORDER:
        implemented = tool_emits.get(src_method, set())
        for code in CODE_ORDER:
            if code not in implemented:      # codes not implemented are not written
                continue
            r = raw.get((build, src_method, code))
            if r is None:
                continue
            E = to_int(r["erepo_count"])
            T = to_int(r["tool_count"])
            rec = to_float(r["recall"])
            prec = to_float(r["precision"])
            if rec is not None and E > 0:
                TP = round(rec * E)
            elif prec is not None and T > 0:
                TP = round(prec * T)
            else:
                TP = 0
            TP = min(TP, E, T) if (E > 0 and T > 0) else TP
            FP = T - TP
            FN = E - TP
            TN = n - TP - FP - FN
            sens = TP / E if E > 0 else None
            spec = TN / (TN + FP) if (TN + FP) > 0 else None
            preci = TP / T if T > 0 else None
            if sens and preci and (sens + preci) > 0:
                f1 = 2 * preci * sens / (preci + sens)
            elif E > 0 and T > 0:
                f1 = 0.0
            else:
                f1 = None
            rows_out.append([build, disp, code, TP, FP, FN, TN,
                             fmt(sens), fmt(spec), fmt(preci), fmt(f1)])

with OUT.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["build", "method", "code", "TP", "FP", "FN", "TN",
                "sensitivity", "specificity", "precision", "f1"])
    w.writerows(rows_out)

print(f"Wrote {OUT}  ({len(rows_out)} rows)")
