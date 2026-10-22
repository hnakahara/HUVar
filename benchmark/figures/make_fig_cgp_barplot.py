#!/usr/bin/env python
"""Figure 5 (CGP cohort; redrawn former Fig. 6, Reviewer 3 Comment 19)

Performance of HUVar in the comprehensive genomic profiling (CGP) cohort,
drawn with the same grouped bars as Fig. 3.
  rows : P/LP, VUS, B/LB (one-versus-rest)
  x axis: 3 prediction configurations (ESM1b+OpenSpliceAI / AM+SpliceAI / REVEL+SpliceAI) x
        gene groups (All / With VCEP / Without VCEP)
  bars : Sensitivity, Specificity, Precision, F1 Score
Reference = institutional expert panel ("Expert panel" in the CGP supplementary table),
HUVar = Bayesian classification. Undefined values (denominator 0) have no bar and are marked "na".

  python make_fig_cgp_barplot.py --s3 "Supplementary Tables.xls" \
      --results-dir <folder with f1_acmg_hg38_{default,alphamissense,revel}_results.tsv> \
      [--out outputs/cgp_barplot_hg38]
Output: <out>.png / <out>.pdf / <out>_metrics.tsv
"""
import argparse
import colorsys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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


# genes in the CGP cohort whose VCEP specification is implemented in HUVar
VCEP_GENES = {"APC", "BRCA1", "BRCA2", "HNF1A", "MLH1", "MSH2", "MSH6", "PALB2",
              "PMS2", "PTEN", "TP53", "VHL"}
TRUTH3 = {"Pathogenic/Likely_pathogenic": "P/LP", "VUS": "VUS", "Benign/Likely_benign": "B/LB"}
PRED3 = {"Pathogenic": "P/LP", "LikelyPathogenic": "P/LP", "VUS": "VUS",
         "UncertainSignificance": "VUS", "LikelyBenign": "B/LB", "Benign": "B/LB"}
CONFIGS = [("default", "ESM1b+\nOpenSpliceAI"), ("alphamissense", "AM+\nSpliceAI"),
           ("revel", "REVEL+\nSpliceAI")]
GROUPS = ["All", "With VCEP", "Without VCEP"]
CATEGORIES = ["P/LP", "VUS", "B/LB"]
METRICS = [("sensitivity", "Sensitivity"), ("specificity", "Specificity"),
           ("precision", "Precision"), ("f1", "F1 Score")]
COLORS = [brighten(c) for c in ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]]


def ovr(t, p, c):
    tp = int(((t == c) & (p == c)).sum()); fp = int(((t != c) & (p == c)).sum())
    fn = int(((t == c) & (p != c)).sum()); tn = int(((t != c) & (p != c)).sum())
    se = tp / (tp + fn) if tp + fn else np.nan
    sp = tn / (tn + fp) if tn + fp else np.nan
    pr = tp / (tp + fp) if tp + fp else np.nan
    f1 = 2 * tp / (2 * tp + fp + fn) if tp + fp + fn else np.nan
    return dict(n_true=tp + fn, TP=tp, FP=fp, FN=fn, TN=tn, sensitivity=se, specificity=sp,
                precision=pr, f1=f1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s3", type=Path, required=True)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--out", default="outputs/cgp_barplot_hg38")
    a = ap.parse_args()

    s3 = pd.read_excel(a.s3, sheet_name="Supplementary Table S3", header=2)
    s3["variant_id"] = (s3.Chromosome.astype(str) + ":" + s3.Location.astype(str) + ":"
                        + s3.Ref.astype(str) + ":" + s3.Alt.astype(str))
    s3["truth"] = s3["Expert panel"].map(TRUTH3)
    s3["group"] = np.where(s3.Gene.isin(VCEP_GENES), "With VCEP", "Without VCEP")

    recs = []
    for cfg, clab in CONFIGS:
        r = pd.read_csv(a.results_dir / f"f1_acmg_hg38_{cfg}_results.tsv", sep="\t", dtype=str,
                        comment="#", usecols=["variant_id", "classification_bayesian"])
        r = r.drop_duplicates("variant_id").set_index("variant_id").classification_bayesian
        d = s3.assign(pred=s3.variant_id.map(r).map(PRED3)).dropna(subset=["truth", "pred"])
        for g in GROUPS:
            sub = d if g == "All" else d[d.group == g]
            acc = float((sub.truth == sub.pred).mean()) if len(sub) else np.nan
            for c in CATEGORIES:
                recs.append(dict(config=cfg, group=g, n=len(sub), genes=sub.Gene.nunique(),
                                 accuracy_3class=acc, category=c, **ovr(sub.truth, sub.pred, c)))
    met = pd.DataFrame(recs)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    met.to_csv(f"{a.out}_metrics.tsv", sep="\t", index=False)

    # x positions: 3 groups per configuration, with a gap between configurations
    xs, xlabels = [], []
    for ci, (cfg, clab) in enumerate(CONFIGS):
        for gi, g in enumerate(GROUPS):
            xs.append(ci * (len(GROUPS) + 0.6) + gi)
            row = met[(met.config == cfg) & (met.group == g)].iloc[0]
            xlabels.append(f"{g}\n(n = {row.n:,})")
    xs = np.array(xs)
    bw = 0.8 / len(METRICS)

    fig, axes = plt.subplots(len(CATEGORIES), 1, figsize=(16, 13), sharex=True)
    for ax, cat in zip(axes, CATEGORIES):
        for mi, (key, mlabel) in enumerate(METRICS):
            vals = []
            for cfg, _ in CONFIGS:
                for g in GROUPS:
                    vals.append(float(met[(met.config == cfg) & (met.group == g)
                                          & (met.category == cat)][key].iloc[0]))
            off = (mi - (len(METRICS) - 1) / 2) * bw
            bars = ax.bar(xs + off, vals, bw, label=mlabel, color=COLORS[mi])
            ax.bar_label(bars, labels=["" if np.isnan(v) else f"{v:.2f}" for v in vals],
                         fontsize=7.5, padding=1)
            for xi, v in zip(xs, vals):
                if np.isnan(v):
                    ax.text(xi + off, 0.02, "na", rotation=90, ha="center", va="bottom",
                            fontsize=9, color=COLORS[mi], fontweight="bold")
        ax.set_title(cat, fontsize=18, fontweight="bold", loc="left")
        ax.set_ylim(0, 1.1)
        ax.set_ylabel("Score")
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        ax.set_axisbelow(True)
        for ci in range(1, len(CONFIGS)):
            ax.axvline(ci * (len(GROUPS) + 0.6) - 0.8, color="#999999", linewidth=0.8)
    axes[-1].set_xticks(xs)
    axes[-1].set_xticklabels(xlabels, fontsize=10.5)
    for ci, (cfg, clab) in enumerate(CONFIGS):
        axes[-1].text(ci * (len(GROUPS) + 0.6) + 1, -0.30, clab.replace("+\n", " + "),
                      transform=axes[-1].get_xaxis_transform(), ha="center", va="top",
                      fontsize=14, fontweight="bold")
    h, l = axes[0].get_legend_handles_labels()
    fig.tight_layout(rect=[0, 0.05, 1, 1.0])
    fig.legend(h, l, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.0), frameon=False,
               fontsize=15)
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=150, bbox_inches="tight")
        print(f"Wrote {a.out}.{ext}")
    plt.close(fig)
    print(met[met.category == "P/LP"][["config", "group", "n", "genes", "sensitivity", "specificity",
                                       "precision", "f1", "accuracy_3class"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
