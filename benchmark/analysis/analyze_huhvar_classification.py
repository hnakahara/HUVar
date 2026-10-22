#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_huhvar_classification.py

Compare the final classifications of each HUVar configuration (default / manual / REVEL / REVEL+manual / SpliceAI / SpliceAI+manual)
with the original (VCF-derived truth) and write, as figures and tables, the
category-level (one-vs-rest for P/LP, B/LB and VUS)
  sensitivity, specificity, precision and F1
and the overall accuracy (3 groups).

Each configuration is evaluated with its two classifiers (2015 / Bayesian) as separate methods.
hg19 and hg38 are written separately (the build is added to the file name).
All erepo_{build}_huhvar_*_results.tsv files in the input folder are detected automatically,
so every existing result file, including the manual configurations, is used (no hard-coded file names).
Builds without a VCF are skipped automatically.
The detected result files and any expected configuration that is missing are reported at run time.

Aggregation rule (same as analyze_benchmark):
  several rows for the same coordinate -> majority vote, ties towards pathogenic (P>LP>VUS>LB>B).

Output (--outdir, {build} = hg19 / hg38):
  huhvar_category_metrics_{build}.tsv   long format (TP/FP/FN/TN and metrics per method x category)
  huhvar_f1_wide_{build}.tsv            wide format (F1 per category and overall accuracy, one row per method)
  huhvar_report_{build}.txt             text summary + 3-group confusion matrices
  huhvar_metrics_{build}_{system}.png   2x2 panels (sensitivity/specificity/precision/F1) per classifier
  huhvar_accuracy_{build}.png           overall accuracy (3 groups) of all methods

Usage:
  python3 analyze_huhvar_classification.py --input ./uploads --outdir ./outputs
  python3 analyze_huhvar_classification.py --no-figs
"""

import argparse
import csv
import glob
import os
import re

import erepo_parsers as P
import analyze_benchmark as AB   # reuse the aggregation and metric logic (GROUPS / aggregate_predictions / category_metrics)

# Known source tokens -> display labels (unknown tokens are shown as is)
SOURCE_LABELS = {
    "default": "default",
    "default_manual": "default+manual",
    "revel": "REVEL",
    "revel_manual": "REVEL+manual",
    "spliceai": "spliceAI",
    "spliceai_manual": "spliceAI+manual",
}
# Display order (other tokens are appended in the order detected)
SOURCE_ORDER = ["default", "default_manual", "revel", "revel_manual",
                "spliceai", "spliceai_manual"]
SYSTEMS = ["2015", "bayesian"]
GROUPS = AB.GROUPS  # ["P/LP", "B/LB", "VUS"]


def discover_sources(input_dir, build, log):
    """Detect erepo_{build}_huhvar_<token>(_results).tsv in input_dir and
    return a list of (method_label, source_label, system, loader).
    Every existing file is used, independent of hard-coded file names
    (with or without a trailing _results / _result).
    """
    prefix = f"erepo_{build}_huhvar_"
    found = {}   # token -> path
    for path in sorted(glob.glob(os.path.join(input_dir, f"{prefix}*.tsv"))):
        base = os.path.basename(path)
        token = base[len(prefix):-len(".tsv")]          # e.g. default_manual_results
        token = re.sub(r"_results?$", "", token)        # strip a trailing _results/_result
        if token:
            found[token] = path

    if not found:
        log(f"  [note] {build}: no HUVar result files found "
            f"(erepo_{build}_huhvar_*_results.tsv)")
        return []

    ordered = [t for t in SOURCE_ORDER if t in found]
    ordered += [t for t in found if t not in SOURCE_ORDER]

    log(f"  HUVar result files detected ({len(found)}):")
    for t in ordered:
        log(f"    - {os.path.basename(found[t])}  -> configuration '{SOURCE_LABELS.get(t, t)}'")
    missing = [t for t in SOURCE_ORDER if t not in found]
    if missing:
        log("  [note] expected configurations not found: " +
            ", ".join(SOURCE_LABELS.get(t, t) for t in missing))

    methods = []
    for t in ordered:
        label = SOURCE_LABELS.get(t, t)
        path = found[t]
        for system in SYSTEMS:
            methods.append((
                f"{label} ({system})", label, system,
                (lambda p=path, k=t, s=system:
                 P.load_huhvar(p, build, f"huhvar_{k}", s)),
            ))
    return methods


def main():
    ap = argparse.ArgumentParser(
        description="Classification performance of all HUVar configurations, per build")
    ap.add_argument("--input", default="inputs")
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--no-figs", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    reg = P.file_registry(args.input)
    orig = P.load_original(reg["original"])

    any_build = False
    for build in ("hg19", "hg38"):
        lines = []

        def out(s=""):
            print(s)
            lines.append(s)

        vcf = reg["builds"][build]["vcf"]
        out("=" * 78)
        out(f"Build {build}")
        out("=" * 78)
        if not os.path.exists(vcf):
            out(f"  VCF ({os.path.basename(vcf)}) not found; skipped."
                f" Place the VCF in the input folder to evaluate this build.")
            out("")
            # keep the report anyway
            with open(os.path.join(args.outdir, f"huhvar_report_{build}.txt"),
                      "w", encoding="utf-8") as f:
                f.write("\n".join(lines) + "\n")
            continue

        any_build = True
        truth = P.build_truth(vcf, orig)
        truth_map = truth["coord2assert"]
        from collections import Counter
        tdist = Counter(truth_map.values())
        out(f"  truth coordinates = {len(truth_map)}  "
            f"(conflicting excluded: {truth['stats']['n_conflict']})")
        out("  truth distribution: " +
            ", ".join(f"{k}={tdist.get(k,0)}" for k in P.CANON_ORDER))
        out("")
        out("  available configurations:")
        methods = discover_sources(args.input, build, out)
        out("")

        # ---- Metrics ----
        # results[label] = dict(per, acc, n, conflicts, source, system)
        results = {}
        for label, src_label, system, loader in methods:
            recs = loader()
            pred, conflicts = AB.aggregate_predictions(recs, truth)
            per, conf3, _conf5, n, acc = AB.category_metrics(truth_map, pred)
            results[label] = dict(per=per, acc=acc, n=n, conflicts=conflicts,
                                  source=src_label, system=system, conf3=conf3)

        if not results:
            out("  no configuration could be evaluated.")
            with open(os.path.join(args.outdir, f"huhvar_report_{build}.txt"),
                      "w", encoding="utf-8") as f:
                f.write("\n".join(lines) + "\n")
            continue

        # ---- Tables (long / wide) ----
        long_path = os.path.join(args.outdir, f"huhvar_category_metrics_{build}.tsv")
        with open(long_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["build", "method", "source", "system", "category",
                        "n_eval", "TP", "FP", "FN", "TN",
                        "sensitivity", "specificity", "precision", "f1"])
            for label, r in results.items():
                for g in GROUPS:
                    m = r["per"][g]
                    w.writerow([build, label, r["source"], r["system"], g, r["n"],
                                m["TP"], m["FP"], m["FN"], m["TN"],
                                f"{m['sens']:.4f}", f"{m['spec']:.4f}",
                                f"{m['prec']:.4f}", f"{m['f1']:.4f}"])

        wide_path = os.path.join(args.outdir, f"huhvar_f1_wide_{build}.tsv")
        with open(wide_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["build", "method", "n_eval"] +
                       [f"F1_{g}" for g in GROUPS] + ["accuracy_3group", "conflicts"])
            for label, r in results.items():
                w.writerow([build, label, r["n"]] +
                           [f"{r['per'][g]['f1']:.4f}" for g in GROUPS] +
                           [f"{r['acc']:.4f}", r["conflicts"]])

        # ---- Report ----
        out("-" * 78)
        out(f"{'method':22s} {'n':>5s} {'acc3':>7s}  "
            f"{'F1(P/LP)':>9s} {'F1(B/LB)':>9s} {'F1(VUS)':>9s}")
        for label, r in results.items():
            out(f"{label:22s} {r['n']:5d} {AB.fmt(r['acc']):>7s}  "
                f"{AB.fmt(r['per']['P/LP']['f1']):>9s} "
                f"{AB.fmt(r['per']['B/LB']['f1']):>9s} "
                f"{AB.fmt(r['per']['VUS']['f1']):>9s}")
        out("")
        # details of each method + 3-group confusion matrix
        for label, r in results.items():
            out(f"### {label}  (evaluated={r['n']}, split coordinates={r['conflicts']})")
            out(f"  {'category':6s} {'sens':>8s} {'spec':>8s} "
                f"{'prec':>8s} {'F1':>8s}  (TP/FP/FN/TN)")
            for g in GROUPS:
                m = r["per"][g]
                out(f"  {g:6s} {AB.fmt(m['sens']):>8s} {AB.fmt(m['spec']):>8s} "
                    f"{AB.fmt(m['prec']):>8s} {AB.fmt(m['f1']):>8s}"
                    f"  ({m['TP']}/{m['FP']}/{m['FN']}/{m['TN']})")
            conf3 = r["conf3"]
            out("  3-group confusion matrix (rows = truth, columns = prediction):")
            out("        " + "".join(f"{g:>8s}" for g in GROUPS))
            for tg in GROUPS:
                out(f"  {tg:>5s} " +
                    "".join(f"{conf3.get((tg, pg), 0):8d}" for pg in GROUPS))
            out("")

        report_path = os.path.join(args.outdir, f"huhvar_report_{build}.txt")
        with open(report_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        print(f"[written] {long_path}")
        print(f"[written] {wide_path}")
        print(f"[written] {report_path}")

        # ---- Figures ----
        if not args.no_figs:
            make_figures(args.outdir, build, results)
            print(f"[figs]    huhvar_metrics_{build}_*.png / huhvar_accuracy_{build}.png")

    if not any_build:
        print("\nNo build could be evaluated (no VCF found).")


def make_figures(outdir, build, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    metrics = [("sens", "Sensitivity"), ("spec", "Specificity"),
               ("prec", "Precision"), ("f1", "F1")]

    # (1) 2x2 panels per classifier (2015 / Bayesian) (metric x category x source)
    for system in SYSTEMS:
        labels = [l for l, r in results.items() if r["system"] == system]
        if not labels:
            continue
        src_labels = [results[l]["source"] for l in labels]
        fig, axes = plt.subplots(2, 2, figsize=(13, 9))
        x = np.arange(len(GROUPS))
        w = 0.8 / max(len(labels), 1)
        off = (len(labels) - 1) / 2.0
        for ax, (mkey, mname) in zip(axes.ravel(), metrics):
            for i, label in enumerate(labels):
                per = results[label]["per"]
                vals = [per[g][mkey] if per[g][mkey] == per[g][mkey] else 0.0
                        for g in GROUPS]
                ax.bar(x + (i - off) * w, vals, w, label=src_labels[i])
            ax.set_xticks(x)
            ax.set_xticklabels(GROUPS)
            ax.set_ylim(0, 1)
            ax.set_title(mname)
            ax.set_ylabel(mname)
        axes.ravel()[0].legend(fontsize=8, ncol=2)
        fig.suptitle(f"HUHVar classification metrics ({build}, {system})",
                     fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        fig.savefig(os.path.join(outdir, f"huhvar_metrics_{build}_{system}.png"),
                    dpi=130)
        plt.close(fig)

    # (2) overall accuracy (3 groups) of all methods (horizontal bars)
    labels = list(results.keys())
    accs = [results[l]["acc"] if results[l]["acc"] == results[l]["acc"] else 0.0
            for l in labels]
    fig, ax = plt.subplots(figsize=(9, max(3, 0.45 * len(labels) + 1)))
    y = np.arange(len(labels))
    ax.barh(y, accs, color="#4c72b0")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Overall accuracy (3-group)")
    ax.set_title(f"HUHVar overall accuracy ({build})")
    for yi, a in zip(y, accs):
        ax.text(min(a + 0.01, 0.98), yi, f"{a*100:.1f}%",
                va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"huhvar_accuracy_{build}.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
