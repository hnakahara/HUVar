#!/usr/bin/env python
"""Bar plots of analysis_criteria_metrics.tsv, one figure per build (hg19/hg38).
Each figure: 4 subplots (sensitivity/specificity/precision/F1).
x axis = implemented codes, grouped bars for 6 methods. No bar for codes a method does not implement. PNG + PDF output."""
import colorsys
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42  # embed TrueType fonts in the PDF (keeps Arial)
plt.rcParams["font.size"] = 15     # default font enlarged (tick labels, y labels, etc.)


def brighten(hexc, s_mul=1.35, v_mul=1.12):
    """Increase saturation and lightness while keeping the hue."""
    hexc = hexc.lstrip("#")
    r, g, b = (int(hexc[i:i + 2], 16) / 255 for i in (0, 2, 4))
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    s, v = min(1.0, s * s_mul), min(1.0, v * v_mul)
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return "#{:02X}{:02X}{:02X}".format(round(r * 255), round(g * 255), round(b * 255))

SRC = Path("outputs/analysis_criteria_metrics.tsv")
OUT_TMPL = "out/criteria_barplot_{build}"
OUT_FORMATS = ["png", "pdf"]

METHOD_ORDER = [
    ("InterVar",                            "InterVar"),
    ("BIAS-2015 (default)",                 "BIAS-2015"),
    ("BIAS-2015 (no PP5/BP6)",              "BIAS-2015 (no PP5/BP6)"),
    ("HUVar ESM1b+OpenSpliceAI",            "HUVar ESM1b+OpenSpliceAI"),
    ("HUVar AlphaMissense+SpliceAI",        "HUVar AM+SpliceAI"),
    ("HUVar REVEL/BayesDel/CADD+SpliceAI",  "HUVar REVEL+SpliceAI"),
]
CODE_ORDER = ["PVS1", "PS1", "PS2", "PS3", "PS4",
              "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
              "PP1", "PP2", "PP3", "PP4", "PP5",
              "BA1", "BS1", "BS2", "BS3", "BS4",
              "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]
BUILD_ORDER = ["hg19", "hg38"]

# codes a tool does not evaluate ("na" is shown only for these combinations; other empty positions are left blank)
_IV_NE = {"PS2", "PS3", "PM3", "PM6", "PP1", "PP4", "BS3", "BS4", "BP2", "BP5"}
_BIAS_NE = {"PS2", "PM3", "PM6", "PP1", "PP4", "BS3", "BS4", "BP2", "BP5"}
_HUVAR_NE = {"PS2", "PS4", "PM3", "PM6", "PP1", "PP4", "BS4", "BP2", "BP5", "PP5", "BP6"}  # HUVar: codes not evaluated automatically (PS4/PP1 curated only, others manual only) + PP5/BP6 (disabled per SVI)
NOT_EVALUATED = {"InterVar": _IV_NE, "BIAS-2015 (default)": _BIAS_NE, "BIAS-2015 (no PP5/BP6)": _BIAS_NE | {"PP5", "BP6"},
                 "HUVar ESM1b+OpenSpliceAI": _HUVAR_NE | {"BS3"}, "HUVar AlphaMissense+SpliceAI": _HUVAR_NE, "HUVar REVEL/BayesDel/CADD+SpliceAI": _HUVAR_NE}
METRICS = [("sensitivity", "Sensitivity"), ("specificity", "Specificity"),
           ("precision", "Precision"), ("f1", "F1 Score")]
# colours of the 6 methods
COLORS = [brighten(c) for c in
          ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"]]

data = {}
with SRC.open(encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        data[(r["build"], r["method"], r["code"])] = r


def val(build, method, code, key):
    r = data.get((build, method, code))
    if r is None:
        return np.nan
    s = r[key].strip()
    return float(s) if s else np.nan


n_methods = len(METHOD_ORDER)
bar_w = 0.8 / n_methods

for build in BUILD_ORDER:
    # codes the eRepo does not use (PP5/BP6 etc.; sensitivity/F1 N/A for all methods) are excluded from the figure
    codes = [c for c in CODE_ORDER
             if any(not np.isnan(val(build, m, c, "f1")) or
                    not np.isnan(val(build, m, c, "sensitivity"))
                    for m, _ in METHOD_ORDER)]
    x = np.arange(len(codes))

    fig, axes = plt.subplots(len(METRICS), 1, figsize=(22, 16), sharex=True)
    for ax, (key, mlabel) in zip(axes, METRICS):
        for mi, (method, _) in enumerate(METHOD_ORDER):
            vals = [val(build, method, c, key) for c in codes]
            offset = (mi - (n_methods - 1) / 2) * bar_w
            ax.bar(x + offset, vals, bar_w,
                   label=METHOD_ORDER[mi][1], color=COLORS[mi])
            # write "na" vertically only for codes a tool does not evaluate
            for xi, c in zip(x, codes):
                if c in NOT_EVALUATED.get(method, set()):
                    ax.text(xi + offset, 0.02, "na", rotation=90, ha="center", va="bottom",
                            fontsize=9, color=COLORS[mi], fontweight="bold")
        ax.set_title(mlabel, fontsize=19, fontweight="bold", loc="left")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Score")
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        ax.set_axisbelow(True)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(codes, rotation=0, ha="center", fontsize=20)
    handles, labels = axes[0].get_legend_handles_labels()
    (lambda *a,**k: None)(f"Per-criterion performance ({build})",
                 fontsize=21, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    fig.legend(handles, labels, ncol=6, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=16)
    base = OUT_TMPL.format(build=build)
    for ext in OUT_FORMATS:
        out = f"{base}.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Wrote {out}")
    plt.close(fig)
