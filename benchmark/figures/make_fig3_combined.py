#!/usr/bin/env python
"""Fig. 3 (revised; Reviewer 3, Comment 23): combines the former Fig. 3 (A) and the former Fig. 5 (B) in one figure.

A: 4 metrics per category (P/LP, VUS, B/LB) x 9 classifiers ... drawn as in make_metrics_barplot.py
B: change in HUVar performance with curated evidence (Base vs +Manual) ... drawn as in make_manual_f1_barplot.py

Input: outputs/analysis_category_metrics.tsv, outputs/huhvar_category_metrics_{build}.tsv
Output: outputs/fig3_combined_{build}.png / .pdf (hg19, hg38)
"""
import colorsys
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 15


def brighten(hexc, s_mul=1.35, v_mul=1.12):
    hexc = hexc.lstrip("#")
    r, g, b = (int(hexc[i:i + 2], 16) / 255 for i in (0, 2, 4))
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    s, v = min(1.0, s * s_mul), min(1.0, v * v_mul)
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return "#{:02X}{:02X}{:02X}".format(round(r * 255), round(g * 255), round(b * 255))


SRC_A = Path("outputs/analysis_category_metrics.tsv")
SRC_B = "outputs/huhvar_category_metrics_{build}.tsv"
OUT_TMPL = "out/fig3_combined_{build}"
OUT_FORMATS = ["png", "pdf"]
BUILD_ORDER = ["hg19", "hg38"]

METRICS = [("sensitivity", "Sensitivity"), ("specificity", "Specificity"),
           ("precision", "Precision"), ("f1", "F1 Score")]
CATEGORIES = ["P/LP", "VUS", "B/LB"]

# ---- A (former Fig. 3) ----
METHOD_ORDER = [
    ("InterVar", "InterVar"),
    ("BIAS-2015", "BIAS-2015"),
    ("BIAS-2015 (no PP5/BP6)", "BIAS-2015\n(no PP5/BP6)"),
    ("HUHVar default (2015)", "HUVar ESM1b+\nOpenSpliceAI\n(2015)"),
    ("HUHVar default (bayesian)", "HUVar ESM1b+\nOpenSpliceAI\n(Bayes)"),
    ("HUHVar spliceAI (2015)", "HUVar AM+\nSpliceAI\n(2015)"),
    ("HUHVar spliceAI (bayesian)", "HUVar AM+\nSpliceAI\n(Bayes)"),
    ("HUHVar REVEL (2015)", "HUVar REVEL+\nSpliceAI\n(2015)"),
    ("HUHVar REVEL (bayesian)", "HUVar REVEL+\nSpliceAI\n(Bayes)"),
]
COLORS_A = [brighten(c) for c in ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]]

# ---- B (former Fig. 5) ----
CONFIGS = [
    ("default (2015)", "default+manual (2015)", "ESM1b+\nOpenSpliceAI\n(2015)"),
    ("default (bayesian)", "default+manual (bayesian)", "ESM1b+\nOpenSpliceAI\n(Bayes)"),
    ("spliceAI (2015)", "spliceAI+manual (2015)", "AM+\nSpliceAI\n(2015)"),
    ("spliceAI (bayesian)", "spliceAI+manual (bayesian)", "AM+\nSpliceAI\n(Bayes)"),
    ("REVEL (2015)", "REVEL+manual (2015)", "REVEL+\nSpliceAI\n(2015)"),
    ("REVEL (bayesian)", "REVEL+manual (bayesian)", "REVEL+\nSpliceAI\n(Bayes)"),
]
COL_BASE = "#B9C4D0"
COL_MANUAL = "#2E67C5"

data_a = {}
with SRC_A.open(encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        data_a[(r["build"], r["method"], r["category"])] = r


def panel_a(sf, build):
    axes = sf.subplots(len(CATEGORIES), 1, sharex=True)
    x = np.arange(len(METHOD_ORDER))
    bw = 0.8 / len(METRICS)
    for ax, cat in zip(axes, CATEGORIES):
        for mi, (key, mlabel) in enumerate(METRICS):
            vals = [float(data_a[(build, m, cat)][key]) for m, _ in METHOD_ORDER]
            bars = ax.bar(x + (mi - (len(METRICS) - 1) / 2) * bw, vals, bw, label=mlabel,
                          color=COLORS_A[mi])
            ax.bar_label(bars, fmt="%.2f", fontsize=8, padding=1)
        ax.set_title(cat, fontsize=18, fontweight="bold", loc="left")
        ax.set_ylim(0, 1.08)
        ax.set_ylabel("Score")
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        ax.set_axisbelow(True)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([lbl for _, lbl in METHOD_ORDER], fontsize=11)
    h, l = axes[0].get_legend_handles_labels()
    sf.legend(h, l, ncol=4, loc="outside lower center", frameon=False, fontsize=16)
    (lambda *a,**k: None)("Classification performance by category", fontsize=20, fontweight="bold")


def panel_b(sf, build):
    data = {}
    with open(SRC_B.format(build=build), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            data[(r["method"], r["category"])] = r
    axes = sf.subplots(len(METRICS), len(CATEGORIES), sharex=True)
    x = np.arange(len(CONFIGS))
    bw = 0.38
    for ri, (mkey, mname) in enumerate(METRICS):
        for ci, cat in enumerate(CATEGORIES):
            ax = axes[ri][ci]
            bv = [float(data[(b, cat)][mkey]) for b, _, _ in CONFIGS]
            mv = [float(data[(m, cat)][mkey]) for _, m, _ in CONFIGS]
            ax.bar(x - bw / 2, bv, bw, label="Base", color=COL_BASE)
            ax.bar(x + bw / 2, mv, bw, label="+Manual", color=COL_MANUAL)
            for xi, b_, m_ in zip(x, bv, mv):
                d = m_ - b_
                ax.annotate(f"{d:+.2f}", (xi + bw / 2, m_), textcoords="offset points",
                            xytext=(0, 3), ha="center", fontsize=12, fontweight="bold",
                            color="#1a7a1a" if d >= 0 else "#b00000")
            ax.set_ylim(0, 1.12)
            ax.grid(axis="y", linestyle="--", alpha=0.4)
            ax.set_axisbelow(True)
            if ri == 0:
                ax.set_title(cat, fontsize=19, fontweight="bold")
            if ci == 0:
                ax.set_ylabel(mname, fontsize=18, fontweight="bold")
    for ci in range(len(CATEGORIES)):
        axes[-1][ci].set_xticks(x)
        axes[-1][ci].set_xticklabels([lbl for _, _, lbl in CONFIGS], fontsize=10)
    h, l = axes[0][1].get_legend_handles_labels()
    sf.legend(h, l, ncol=2, loc="outside lower center", frameon=False, fontsize=17)
    (lambda *a,**k: None)("Effect of manual evidence on category metrics (HUVar)", fontsize=20,
                fontweight="bold")


for build in BUILD_ORDER:
    fig = plt.figure(figsize=(20, 29), layout="constrained")
    fig.get_layout_engine().set(h_pad=0.08, hspace=0.04)
    sfa, sfb = fig.subfigures(2, 1, height_ratios=[12.5, 16])
    panel_a(sfa, build)
    panel_b(sfb, build)
    for sf, lab in ((sfa, "A"), (sfb, "B")):
        sf.text(0.0, 1.0, lab, fontsize=34, fontweight="bold", va="top", ha="left")
    base = OUT_TMPL.format(build=build)
    for ext in OUT_FORMATS:
        fig.savefig(f"{base}.{ext}", dpi=150, bbox_inches="tight")
        print(f"Wrote {base}.{ext}")
    plt.close(fig)
