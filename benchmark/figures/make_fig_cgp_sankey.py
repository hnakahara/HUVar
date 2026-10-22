#!/usr/bin/env python
"""CGP cohort figure (redrawn former Fig. 6; Reviewer 3 Comment 19) - Sankey diagram.

Shows, in the same style as Supplementary Figure S6, how the 3-class HUVar (Bayesian) classifications correspond to those of the institutional expert panel.
Rows = prediction configurations (ESM1b+OpenSpliceAI / AM+SpliceAI / REVEL+SpliceAI),
columns = all genes / genes with a VCEP specification / genes without. The height of each panel is normalised to the number of variants
(comparison of proportions). Flows are coloured by the left-hand (HUVar) class.

With --f1cdx, three stages are drawn as in S6: F1CDx -> HUVar -> Expert panel
(TSV with variant_id and class columns; class is P/LP, VUS, B/LB or the long form).

  python make_fig_cgp_sankey.py --s3 "Supplementary Tables.xls" \
      --results-dir <folder with f1_acmg_hg38_{default,alphamissense,revel}_results.tsv> \
      [--f1cdx f1cdx_classes.tsv] [--groups all vcep] [--out outputs/cgp_sankey_hg38]
Output: <out>.png / <out>.pdf / <out>_flows.tsv
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 11

VCEP_GENES = {"APC", "BRCA1", "BRCA2", "HNF1A", "MLH1", "MSH2", "MSH6", "PALB2",
              "PMS2", "PTEN", "TP53", "VHL"}
NORM = {"Pathogenic/Likely_pathogenic": "P/LP", "Pathogenic": "P/LP", "LikelyPathogenic": "P/LP",
        "Likely_pathogenic": "P/LP", "Likely P": "P/LP", "P": "P/LP", "LP": "P/LP", "P/LP": "P/LP",
        "VUS": "VUS", "UncertainSignificance": "VUS", "Uncertain_significance": "VUS",
        "Benign/Likely_benign": "B/LB", "Benign": "B/LB", "LikelyBenign": "B/LB",
        "Likely_benign": "B/LB", "B": "B/LB", "LB": "B/LB", "B/LB": "B/LB"}
CATS = ["P/LP", "VUS", "B/LB"]
NODE_COL = {"P/LP": "#E0565B", "VUS": "#ADA59E", "B/LB": "#4A76A8"}
FLOW_COL = {"P/LP": "#EBA0A3", "VUS": "#D6D1CC", "B/LB": "#9DB6D3"}
CONFIGS = [("default", "ESM1b + OpenSpliceAI"), ("alphamissense", "AM + SpliceAI"),
           ("revel", "REVEL + SpliceAI")]
GROUPS = {"all": ("All genes", None), "with": ("Genes with VCEP specifications", True),
          "without": ("Genes without VCEP specifications", False)}

GAP = 0.04      # gap between class nodes (fraction of the panel height 1)
NODE_W = 0.025  # node width (x units)


def layout(counts, total):
    """(bottom, top) of each class node, from the top: P/LP, VUS, B/LB. Height = proportion."""
    present = [c for c in CATS if counts.get(c, 0) > 0]
    usable = 1.0 - GAP * (len(present) - 1)
    pos, y = {}, 1.0
    for c in CATS:
        h = usable * counts.get(c, 0) / total if total else 0
        pos[c] = (y - h, y)
        if counts.get(c, 0) > 0:
            y -= h + GAP
    return pos


def ribbon(ax, x0, x1, y0a, y0b, y1a, y1b, color):
    """Smooth band from [y0a, y0b] at x0 to [y1a, y1b] at x1."""
    t = np.linspace(0, 1, 60)
    s = 3 * t ** 2 - 2 * t ** 3
    xs = x0 + (x1 - x0) * t
    lo = y0a + (y1a - y0a) * s
    hi = y0b + (y1b - y0b) * s
    ax.fill_between(xs, lo, hi, color=color, alpha=0.75, linewidth=0)


def draw(ax, df, stages, title):
    """stages: list of column names (left to right). Each column of df is P/LP / VUS / B/LB."""
    n = len(df)
    xs = np.linspace(0, 1, len(stages))
    pos = {}
    for st in stages:
        pos[st] = layout(df[st].value_counts().to_dict(), n)
    usable = {st: 1.0 - GAP * (sum(df[st].value_counts().get(c, 0) > 0 for c in CATS) - 1)
              for st in stages}
    rows = []
    for k in range(len(stages) - 1):
        a, b = stages[k], stages[k + 1]
        ct = pd.crosstab(df[a], df[b]).reindex(index=CATS, columns=CATS, fill_value=0)
        out_off = {c: pos[a][c][1] for c in CATS}   # current position within the source node (from the top)
        in_off = {c: pos[b][c][1] for c in CATS}
        for sa in CATS:
            for sb in CATS:
                v = int(ct.at[sa, sb])
                if v == 0:
                    continue
                ha = usable[a] * v / n
                hb = usable[b] * v / n
                y0b, y0a = out_off[sa], out_off[sa] - ha
                y1b, y1a = in_off[sb], in_off[sb] - hb
                out_off[sa] -= ha
                in_off[sb] -= hb
                ribbon(ax, xs[k] + NODE_W / 2, xs[k + 1] - NODE_W / 2, y0a, y0b, y1a, y1b,
                       FLOW_COL[sa])
                rows.append(dict(panel=title, source_stage=a, source=sa, target_stage=b,
                                 target=sb, n=v))
    for i, st in enumerate(stages):
        vc = df[st].value_counts()
        for c in CATS:
            if vc.get(c, 0) == 0:
                continue
            lo, hi = pos[st][c]
            ax.add_patch(plt.Rectangle((xs[i] - NODE_W / 2, lo), NODE_W, hi - lo,
                                       color=NODE_COL[c], zorder=3, linewidth=0))
            ha, xt = ("left", xs[i] + NODE_W / 2 + 0.01) if i < len(stages) - 1 else \
                     ("right", xs[i] - NODE_W / 2 - 0.01)
            ax.text(xt, (lo + hi) / 2, f"{c} ({vc[c]:,})", ha=ha, va="center", fontsize=9,
                    zorder=4, bbox=dict(facecolor="white", alpha=0.7, edgecolor="none",
                                        pad=1.2))
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.02, 1.02)
    ax.axis("off")
    ax.set_title(title, fontsize=10.5, pad=14)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s3", type=Path, required=True)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--f1cdx", type=Path, default=None)
    ap.add_argument("--groups", nargs="+", default=["all", "with", "without"],
                    choices=list(GROUPS))
    ap.add_argument("--out", default="outputs/cgp_sankey_hg38")
    a = ap.parse_args()

    s3 = pd.read_excel(a.s3, sheet_name="Supplementary Table S3", header=2)
    s3["variant_id"] = (s3.Chromosome.astype(str) + ":" + s3.Location.astype(str) + ":"
                        + s3.Ref.astype(str) + ":" + s3.Alt.astype(str))
    s3["Expert panel"] = s3["Expert panel"].map(NORM)
    s3["vcep"] = s3.Gene.isin(VCEP_GENES)
    stages_base = ["HUVar", "Expert panel"]
    if a.f1cdx is not None:
        f = pd.read_csv(a.f1cdx, sep="\t", dtype=str)
        s3["F1CDx"] = s3.variant_id.map(f.drop_duplicates("variant_id").set_index("variant_id")["class"]
                                        .map(lambda x: NORM.get(str(x).strip(), np.nan)))
        stages_base = ["F1CDx"] + stages_base

    fig, axes = plt.subplots(len(CONFIGS), len(a.groups),
                             figsize=(4.6 * len(a.groups) + 0.6, 4.6 * len(CONFIGS)),
                             squeeze=False)
    rows = []
    for ri, (cfg, clab) in enumerate(CONFIGS):
        r = pd.read_csv(a.results_dir / f"f1_acmg_hg38_{cfg}_results.tsv", sep="\t", dtype=str,
                        comment="#", usecols=["variant_id", "classification_bayesian"])
        r = r.drop_duplicates("variant_id").set_index("variant_id").classification_bayesian
        d = s3.assign(HUVar=s3.variant_id.map(r).map(NORM)).dropna(subset=stages_base)
        for gi, g in enumerate(a.groups):
            glab, flag = GROUPS[g]
            sub = d if flag is None else d[d.vcep == flag]
            acc = (sub.HUVar == sub["Expert panel"]).mean()
            title = f"{clab}\nn = {len(sub):,} ({sub.Gene.nunique()} genes); " \
                    f"HUVar = expert panel {acc * 100:.1f}%"
            if "F1CDx" in sub:
                title += f"\nF1CDx = expert panel {(sub.F1CDx == sub['Expert panel']).mean() * 100:.1f}%"
            for rr in draw(axes[ri][gi], sub, stages_base, title):
                rows.append(dict(config=cfg, group=g, **rr))
    for ri in range(len(CONFIGS)):
        for gi in range(len(a.groups)):
            ax = axes[ri][gi]
            for k, st in enumerate(stages_base):
                x = np.linspace(0, 1, len(stages_base))[k]
                ax.text(x, 1.035, "HUVar" if st == "HUVar" else st, ha="center", va="bottom",
                        fontsize=9, color="#555555")
    fig.tight_layout(h_pad=2.6, w_pad=2.5, rect=[0, 0, 1, 0.965])
    for gi, g in enumerate(a.groups):
        bb = axes[0][gi].get_position()
        fig.text((bb.x0 + bb.x1) / 2, 0.985, GROUPS[g][0], ha="center", va="top", fontsize=13.5,
                 fontweight="bold")
    pd.DataFrame(rows).to_csv(f"{a.out}_flows.tsv", sep="\t", index=False)
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=200, bbox_inches="tight", facecolor="white")
        print(f"Wrote {a.out}.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    main()
