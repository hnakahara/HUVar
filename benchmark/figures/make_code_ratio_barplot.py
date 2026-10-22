#!/usr/bin/env python
"""Figure 4 (revised; Reviewer 2, Comment 10): the grouped bars of make_code_counts_barplot.py, with the
number of variants assigned each criterion normalised to the number in the eRepo (truth).
y = assignments by the tool / assignments in the eRepo (linear); 1 = same as the eRepo. Bars above YMAX are
truncated at the limit and labelled with the actual ratio.
No bar for 0 assignments. "na" is written vertically near the x axis only for codes a tool does not evaluate (NOT_EVALUATED). PP5/BP6, which the eRepo does not use,
have no defined ratio and are excluded; their counts are given in the note. The raw counts are in the figure of make_code_counts_barplot.py
(supplementary) and in outputs/analysis_acmg_code_counts.tsv (supplementary table). PNG + PDF output.

  python make_code_ratio_barplot.py
"""

import colorsys
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42  # embed TrueType fonts in the PDF (keeps Arial)
plt.rcParams["font.size"] = 15


def brighten(hexc, s_mul=1.35, v_mul=1.12):
    """Increase saturation and lightness while keeping the hue (same as make_code_counts_barplot.py)."""
    hexc = hexc.lstrip("#")
    r, g, b = (int(hexc[i : i + 2], 16) / 255 for i in (0, 2, 4))
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    s, v = min(1.0, s * s_mul), min(1.0, v * v_mul)
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return "#{:02X}{:02X}{:02X}".format(round(r * 255), round(g * 255), round(b * 255))


SRC = Path("outputs/analysis_acmg_code_counts.tsv")
OUT_TMPL = "out/code_ratio_barplot_{build}"
OUT_FORMATS = ["png", "pdf"]

# same display names and colours as make_code_counts_barplot.py (the eRepo is the reference, so it has no bar)
SOURCE_ORDER = [
    ("InterVar", "InterVar", "#4C72B0"),
    ("BIAS-2015", "BIAS-2015", "#DD8452"),
    ("BIAS_noPP5BP6", "BIAS-2015 (no PP5/BP6)", "#55A868"),
    ("HUHVar_default", "HUVar ESM1b+OpenSpliceAI", "#C44E52"),
    ("HUHVar_spliceai", "HUVar AM+SpliceAI", "#8172B3"),
    ("HUHVar_revel", "HUVar REVEL+SpliceAI", "#937860"),
]
SOURCE_ORDER = [(s, d, brighten(c)) for s, d, c in SOURCE_ORDER]
CODE_ORDER = ["PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
              "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4",
              "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]
BUILD_ORDER = ["hg19", "hg38"]

# codes a tool does not evaluate ("na" is shown only for these combinations; other empty positions are left blank)
_IV_NE = {"PS2", "PS3", "PM3", "PM6", "PP1", "PP4", "BS3", "BS4", "BP2", "BP5"}
_BIAS_NE = {"PS2", "PM3", "PM6", "PP1", "PP4", "BS3", "BS4", "BP2", "BP5"}
_HUVAR_NE = {"PS2", "PS4", "PM3", "PM6", "PP1", "PP4", "BS4", "BP2", "BP5", "PP5", "BP6"}  # HUVar: codes not evaluated automatically (PS4/PP1 curated only, others manual only) + PP5/BP6 (disabled per SVI)
NOT_EVALUATED = {"InterVar": _IV_NE, "BIAS-2015": _BIAS_NE, "BIAS_noPP5BP6": _BIAS_NE | {"PP5", "BP6"},
                 "HUHVar_default": _HUVAR_NE | {"BS3"}, "HUHVar_spliceai": _HUVAR_NE, "HUHVar_revel": _HUVAR_NE}

counts = {}
with SRC.open(encoding="utf-8") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        counts[(r["build"], r["code"], r["source"])] = int(r["count"])

n_src = len(SOURCE_ORDER)
bar_w = 0.85 / n_src
YMAX = 4.0  # ratios above this are truncated and labelled

for build in BUILD_ORDER:
    codes = [c for c in CODE_ORDER if counts.get((build, c, "eRepo"), 0) > 0]
    no_ref = [c for c in CODE_ORDER if c not in codes]
    x = np.arange(len(codes))
    fig, ax = plt.subplots(figsize=(22, 9))
    for si, (src, disp, color) in enumerate(SOURCE_ORDER):
        vals, zeros = [], []
        for j, code in enumerate(codes):
            ref = counts[(build, code, "eRepo")]
            c = counts.get((build, code, src), 0)
            if c > 0:
                v = c / ref
                if v > YMAX:
                    ax.text(x[j] + (si - (n_src - 1) / 2) * bar_w, YMAX * 1.005, f"{v:.1f}",
                            ha="center", va="bottom", fontsize=9, color=color, rotation=90,
                            fontweight="bold")
                vals.append(min(v, YMAX))
            else:
                vals.append(np.nan)  # no bar for 0
                zeros.append(j)
        offset = (si - (n_src - 1) / 2) * bar_w
        ax.bar(x + offset, vals, bar_w, label=disp, color=color)
        for j in zeros:
            if codes[j] in NOT_EVALUATED.get(src, set()):
                ax.text(x[j] + offset, 0.03, "na", rotation=90, ha="center", va="bottom",
                        fontsize=9, color=color, fontweight="bold")
    ax.axhline(1, color="#333333", linewidth=1.2)
    ax.set_ylim(0, YMAX * 1.12)
    ax.set_yticks(np.arange(0, YMAX + 0.01, 0.5))
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c}\n({counts[(build, c, 'eRepo')]:,})" for c in codes],
                       rotation=0, ha="center", fontsize=12)
    ax.set_xlim(-0.6, len(codes) - 0.4)
    ax.set_ylabel("Variants assigned / eRepo")
    ax.grid(axis="y", which="major", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    handles, labels = ax.get_legend_handles_labels()
    (lambda *a,**k: None)(f"ACMG criterion assignment relative to the eRepo ({build})",
                 fontsize=21, fontweight="bold", y=0.99)
    fig.legend(handles, labels, ncol=6, loc="lower center", bbox_to_anchor=(0.5, 0.06),
               frameon=False, fontsize=16)
    note = ("Numbers under the criteria: variants assigned the criterion in the eRepo. "
            "Ratio 1 (black line) = same as the eRepo; bars above the axis limit are truncated and "
            "labelled with their ratio. na: not evaluated by the tool (HUVar: criteria applied only from curated or manual evidence, and PP5/BP6).\nNot used by the eRepo (ratio undefined): "
            + "; ".join(f"{c} (" + ", ".join(f"{d} {counts.get((build, c, s), 0):,}"
                                              for s, d, _ in SOURCE_ORDER
                                              if counts.get((build, c, s), 0) > 0) + ")"
                        for c in no_ref) + ".")
    (lambda *a,**k: None)(0.05, 0.015, note, fontsize=11, color="#444444")
    fig.subplots_adjust(top=0.93, left=0.06, right=0.99, bottom=0.25)
    base = OUT_TMPL.format(build=build)
    for ext in OUT_FORMATS:
        out = f"{base}.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Wrote {out}")
    plt.close(fig)
