#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_benchmark.py

Evaluate each tool against the eRepo, using the truth coordinates derived from the VCF.
As in the BIAS-2015 paper, the following are computed:

  (A) Category-level performance for the three groups P/LP, B/LB and VUS:
        sensitivity, specificity, precision and F1 (one-vs-rest)
      + 3-group confusion matrix / 5-category confusion matrix / overall accuracy
  (B) Number of variants assigned each ACMG code (eRepo / InterVar / BIAS / HUVar)
  (C) F1 per ACMG code (presence-based, with the eRepo code as the truth)

Coordinate aggregation:
  When a tool has several rows for the same coordinate, the representative call is chosen by
  majority vote, ties broken towards pathogenic (P>LP>VUS>LB>B). The number of split coordinates is reported.

HUVar is evaluated with its two classifiers (2015 / Bayesian) as separate methods. The code columns do not
depend on the classifier, so code counts are treated as two sources (default / spliceai).

Usage:
    python3 analyze_benchmark.py --input <uploads> --outdir <out> [--no-figs]
"""

import argparse
import csv
import os
from collections import Counter, defaultdict

import erepo_parsers as P

# 5 categories -> 3 groups
GROUP = {"P": "P/LP", "LP": "P/LP", "VUS": "VUS", "LB": "B/LB", "B": "B/LB"}
GROUPS = ["P/LP", "B/LB", "VUS"]
# Tie-break for the representative call (towards pathogenic)
SEV = {"P": 5, "LP": 4, "VUS": 3, "LB": 2, "B": 1}


def methods_for(reg, build):
    p = reg["builds"][build]
    return [
        ("InterVar", lambda: P.load_intervar(p["intervar"], build)),
        ("BIAS-2015", lambda: P.load_bias(p["bias"], build)),
        ("BIAS-2015 (no PP5/BP6)",
         lambda: P.load_bias(p["bias_nopp5bp6"], build, "bias_nopp5bp6")),
        ("HUHVar default (2015)",
         lambda: P.load_huhvar(p["huhvar_default"], build, "huhvar_default", "2015")),
        ("HUHVar default (bayesian)",
         lambda: P.load_huhvar(p["huhvar_default"], build, "huhvar_default", "bayesian")),
        ("HUHVar REVEL (2015)",
         lambda: P.load_huhvar(p["huhvar_revel"], build, "huhvar_revel", "2015")),
        ("HUHVar REVEL (bayesian)",
         lambda: P.load_huhvar(p["huhvar_revel"], build, "huhvar_revel", "bayesian")),
        ("HUHVar spliceAI (2015)",
         lambda: P.load_huhvar(p["huhvar_spliceai"], build, "huhvar_spliceai", "2015")),
        ("HUHVar spliceAI (bayesian)",
         lambda: P.load_huhvar(p["huhvar_spliceai"], build, "huhvar_spliceai", "bayesian")),
    ]


def aggregate_predictions(recs, truth):
    """coord -> representative 5-category prediction. Also returns the number of split coordinates."""
    by = defaultdict(list)
    for r in recs:
        a, coord = P.resolve_truth(r, truth)
        if coord is not None and a is not None and r["cls"] is not None:
            by[coord].append(r["cls"])
    pred = {}
    conflicts = 0
    for coord, lst in by.items():
        if len(set(lst)) > 1:
            conflicts += 1
        cnt = Counter(lst)
        top = max(cnt.values())
        cand = [c for c in lst if cnt[c] == top]
        pred[coord] = max(cand, key=lambda c: SEV[c])
    return pred, conflicts


def category_metrics(truth_map, pred_map):
    coords = [c for c in truth_map if c in pred_map]
    n = len(coords)
    conf3 = defaultdict(int)   # (truth_group, pred_group)
    conf5 = defaultdict(int)   # (truth5, pred5)
    for c in coords:
        t5, p5 = truth_map[c], pred_map[c]
        conf3[(GROUP[t5], GROUP[p5])] += 1
        conf5[(t5, p5)] += 1
    per = {}
    for g in GROUPS:
        TP = conf3[(g, g)]
        FP = sum(conf3[(t, g)] for t in GROUPS if t != g)
        FN = sum(conf3[(g, p)] for p in GROUPS if p != g)
        TN = n - TP - FP - FN
        sens = TP / (TP + FN) if (TP + FN) else float("nan")
        spec = TN / (TN + FP) if (TN + FP) else float("nan")
        prec = TP / (TP + FP) if (TP + FP) else float("nan")
        f1 = (2 * prec * sens / (prec + sens)
              if (prec + sens) and prec == prec and sens == sens else float("nan"))
        per[g] = dict(TP=TP, FP=FP, FN=FN, TN=TN,
                      sens=sens, spec=spec, prec=prec, f1=f1)
    acc = sum(conf3[(g, g)] for g in GROUPS) / n if n else float("nan")
    return per, conf3, conf5, n, acc


def erepo_coord_codes(reg, build, orig):
    """Build the set of eRepo codes for each truth coordinate (via the VCF NAME)."""
    name2codes = {}
    for r in orig["records"]:
        if r["hgvsc_raw"]:
            name2codes.setdefault(f"{r['gene']}:{r['hgvsc_raw']}", set()).update(r["codes"])
    out = defaultdict(set)
    for chrom, pos, ref, alt, name in P.load_vcf(reg["builds"][build]["vcf"]):
        coord = f"{chrom}:{pos}:{ref}:{alt}"
        out[coord].update(name2codes.get(name, set()))
    return out


def tool_coord_codes(recs, truth):
    out = defaultdict(set)
    for r in recs:
        _a, coord = P.resolve_truth(r, truth)
        if coord is not None:
            out[coord].update(r.get("codes", []))
    return out


def fmt(x):
    return "n/a" if x != x else f"{100*x:.2f}%"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="inputs")
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--no-figs", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    reg = P.file_registry(args.input)
    orig = P.load_original(reg["original"])

    report = []

    def out(s=""):
        print(s)
        report.append(s)

    # Output files
    catf = open(os.path.join(args.outdir, "analysis_category_metrics.tsv"),
                "w", newline="", encoding="utf-8")
    cw = csv.writer(catf, delimiter="\t")
    cw.writerow(["build", "method", "category", "n_eval", "TP", "FP", "FN", "TN",
                 "sensitivity", "specificity", "precision", "f1"])

    conf_lines = []
    code_count_rows = []   # (build, code, source, count)
    code_f1_rows = []      # (build, tool, code, erepo_n, tool_n, precision, recall, f1)

    summary = defaultdict(dict)  # (build,method)->{cat:metrics}

    for build in ("hg19", "hg38"):
        truth = P.build_truth(reg["builds"][build]["vcf"], orig)
        truth_map = truth["coord2assert"]
        n_truth = len(truth_map)
        out("=" * 78)
        out(f"Build {build} : truth variants = {n_truth}  (conflicting excluded: {truth['stats']['n_conflict']})")
        tdist = Counter(truth_map.values())
        out("  truth distribution: " + ", ".join(f"{k}={tdist.get(k,0)}" for k in P.CANON_ORDER))
        out("=" * 78)

        # ---- (A) Per category ----
        loaded = {}
        for name, loader in methods_for(reg, build):
            recs = loader()
            loaded[name] = recs
            pred, conflicts = aggregate_predictions(recs, truth)
            per, conf3, conf5, n, acc = category_metrics(truth_map, pred)
            summary[(build, name)] = per
            out(f"\n[{name}]  evaluated={n}  overall accuracy (3 groups)={fmt(acc)}  split coordinates={conflicts}")
            out(f"  {'category':6s} {'sens':>8s} {'spec':>8s} {'prec':>8s} {'F1':>8s}"
                f"  (TP/FP/FN/TN)")
            for g in GROUPS:
                m = per[g]
                out(f"  {g:6s} {fmt(m['sens']):>8s} {fmt(m['spec']):>8s} "
                    f"{fmt(m['prec']):>8s} {fmt(m['f1']):>8s}"
                    f"  ({m['TP']}/{m['FP']}/{m['FN']}/{m['TN']})")
                cw.writerow([build, name, g, n, m["TP"], m["FP"], m["FN"], m["TN"],
                             f"{m['sens']:.4f}", f"{m['spec']:.4f}",
                             f"{m['prec']:.4f}", f"{m['f1']:.4f}"])
            # Confusion matrices to a separate file
            conf_lines.append(f"### {build} / {name}  (n={n})")
            conf_lines.append("3-group confusion matrix (rows = truth, columns = prediction):")
            conf_lines.append("        " + "".join(f"{g:>8s}" for g in GROUPS))
            for tg in GROUPS:
                conf_lines.append(f"{tg:>7s} " +
                                  "".join(f"{conf3.get((tg,pg),0):8d}" for pg in GROUPS))
            conf_lines.append("5-category confusion matrix (rows = truth, columns = prediction):")
            conf_lines.append("       " + "".join(f"{c:>7s}" for c in P.CANON_ORDER))
            for tc in P.CANON_ORDER:
                conf_lines.append(f"{tc:>6s} " +
                                  "".join(f"{conf5.get((tc,pc),0):7d}" for pc in P.CANON_ORDER))
            conf_lines.append("")

        # ---- (B)/(C) ACMG codes ----
        erepo_cc = erepo_coord_codes(reg, build, orig)
        code_sources = {
            "eRepo": erepo_cc,
            "InterVar": tool_coord_codes(loaded["InterVar"], truth),
            "BIAS-2015": tool_coord_codes(loaded["BIAS-2015"], truth),
            "BIAS_noPP5BP6": tool_coord_codes(loaded["BIAS-2015 (no PP5/BP6)"], truth),
            "HUHVar_default": tool_coord_codes(loaded["HUHVar default (2015)"], truth),
            "HUHVar_revel": tool_coord_codes(loaded["HUHVar REVEL (2015)"], truth),
            "HUHVar_spliceai": tool_coord_codes(loaded["HUHVar spliceAI (2015)"], truth),
        }
        coords = list(truth_map.keys())
        # (B) Assignment counts
        for code in P.ACMG_CODES:
            for src, cc in code_sources.items():
                cnt = sum(1 for c in coords if code in cc.get(c, ()))
                code_count_rows.append((build, code, src, cnt))
        # (C) F1 per code (eRepo as truth)
        for tool in ["InterVar", "BIAS-2015", "BIAS_noPP5BP6",
                     "HUHVar_default", "HUHVar_revel", "HUHVar_spliceai"]:
            cc = code_sources[tool]
            for code in P.ACMG_CODES:
                TP = FP = FN = 0
                en = 0
                for c in coords:
                    t = code in erepo_cc.get(c, ())
                    p = code in cc.get(c, ())
                    en += t
                    if t and p:
                        TP += 1
                    elif p and not t:
                        FP += 1
                    elif t and not p:
                        FN += 1
                prec = TP / (TP + FP) if (TP + FP) else float("nan")
                rec = TP / (TP + FN) if (TP + FN) else float("nan")
                f1 = (2 * prec * rec / (prec + rec)
                      if (prec == prec and rec == rec and (prec + rec)) else float("nan"))
                code_f1_rows.append((build, tool, code, en, TP + FP,
                                     prec, rec, f1))

    catf.close()
    with open(os.path.join(args.outdir, "analysis_confusion.txt"),
              "w", encoding="utf-8") as f:
        f.write("\n".join(conf_lines) + "\n")

    with open(os.path.join(args.outdir, "analysis_acmg_code_counts.tsv"),
              "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["build", "code", "source", "count"])
        w.writerows(code_count_rows)

    with open(os.path.join(args.outdir, "analysis_code_f1.tsv"),
              "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["build", "tool", "code", "erepo_count", "tool_count",
                    "precision", "recall", "f1"])
        for b, t, c, en, tn, pr, rc, f1 in code_f1_rows:
            w.writerow([b, t, c, en, tn,
                        "" if pr != pr else f"{pr:.4f}",
                        "" if rc != rc else f"{rc:.4f}",
                        "" if f1 != f1 else f"{f1:.4f}"])

    # ---- Figures ----
    if not args.no_figs:
        make_figures(args.outdir, summary, code_count_rows, code_f1_rows)
        out("\n[Figures] category_f1_*.png / acmg_code_counts_*.png / code_f1_*.png written")

    report_path = os.path.join(args.outdir, "analysis_report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report) + "\n")
    print(f"\n[written] {report_path}")


def make_figures(outdir, summary, code_count_rows, code_f1_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    methods = ["InterVar", "BIAS-2015", "BIAS-2015 (no PP5/BP6)",
               "HUHVar default (2015)", "HUHVar default (bayesian)",
               "HUHVar REVEL (2015)", "HUHVar REVEL (bayesian)",
               "HUHVar spliceAI (2015)", "HUHVar spliceAI (bayesian)"]
    cats = GROUPS

    # (A) F1 per category
    for build in ("hg19", "hg38"):
        fig, ax = plt.subplots(figsize=(12, 5))
        x = np.arange(len(cats))
        w = 0.8 / len(methods)
        off = (len(methods) - 1) / 2.0
        for i, m in enumerate(methods):
            vals = [summary[(build, m)][c]["f1"] for c in cats]
            ax.bar(x + (i - off) * w, vals, w, label=m)
        ax.set_xticks(x)
        ax.set_xticklabels(cats)
        ax.set_ylabel("F1")
        ax.set_ylim(0, 1)
        ax.set_title(f"Category F1 ({build})")
        ax.legend(fontsize=7, ncol=3)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"category_f1_{build}.png"), dpi=130)
        plt.close(fig)

    # (B) ACMG code assignment counts
    sources = ["eRepo", "InterVar", "BIAS-2015", "BIAS_noPP5BP6",
               "HUHVar_default", "HUHVar_revel", "HUHVar_spliceai"]
    for build in ("hg19", "hg38"):
        data = defaultdict(dict)
        for b, code, src, cnt in code_count_rows:
            if b == build:
                data[src][code] = cnt
        fig, ax = plt.subplots(figsize=(14, 5))
        x = np.arange(len(P.ACMG_CODES))
        w = 0.85 / len(sources)
        off = (len(sources) - 1) / 2.0
        for i, src in enumerate(sources):
            vals = [data[src].get(code, 0) for code in P.ACMG_CODES]
            ax.bar(x + (i - off) * w, vals, w, label=src)
        ax.set_xticks(x)
        ax.set_xticklabels(P.ACMG_CODES, rotation=90, fontsize=7)
        ax.set_ylabel("variants with code")
        ax.set_title(f"ACMG code assignment counts ({build})")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"acmg_code_counts_{build}.png"), dpi=130)
        plt.close(fig)

    # (C) F1 per code
    tools = ["InterVar", "BIAS-2015", "BIAS_noPP5BP6",
             "HUHVar_default", "HUHVar_revel", "HUHVar_spliceai"]
    for build in ("hg19", "hg38"):
        data = defaultdict(dict)
        for b, t, code, en, tn, pr, rc, f1 in code_f1_rows:
            if b == build:
                data[t][code] = (f1 if f1 == f1 else 0.0)
        fig, ax = plt.subplots(figsize=(14, 5))
        x = np.arange(len(P.ACMG_CODES))
        w = 0.85 / len(tools)
        off = (len(tools) - 1) / 2.0
        for i, t in enumerate(tools):
            vals = [data[t].get(code, 0.0) for code in P.ACMG_CODES]
            ax.bar(x + (i - off) * w, vals, w, label=t)
        ax.set_xticks(x)
        ax.set_xticklabels(P.ACMG_CODES, rotation=90, fontsize=7)
        ax.set_ylabel("F1 vs eRepo code presence")
        ax.set_ylim(0, 1)
        ax.set_title(f"Per-code F1 vs eRepo ({build})")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"code_f1_{build}.png"), dpi=130)
        plt.close(fig)


if __name__ == "__main__":
    main()
