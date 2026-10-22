#!/usr/bin/env python
"""Supplementary Figure for Reviewer 3, Comment 15: concordance with the eRepo by
gene (exact / confidence difference / significant conflict / opposite), as
100% stacked horizontal bars, genes sorted by exact concordance.

Agreement categories use the blue arm (exact = dark, confidence = light) and
conflicts the orange arm (significant = light, opposite = dark) of a diverging
pair, so the agreement/conflict boundary reads at a glance.

  python make_fig_gene_concordance.py --summary confusion_excl/confusion_gene_summary.tsv \
      --config revel_manual --min-n 10 --out fig_gene_concordance
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "erepo_benchmark"))
from figstyle import INK, INK2, plt, save  # noqa: E402

CATS = [("exact", "Exact", "#184f95"), ("confidence", "Confidence (P↔LP, B↔LB)", "#86b6ef"),
        ("significant", "Significant (vs VUS)", "#f3b493"), ("opposite", "Opposite (P/LP↔B/LB)", "#c4441a")]
LABEL = {"revel_manual": "REVEL/CADD/BayesDel + SpliceAI with eRepo curated evidence",
         "revel": "REVEL/CADD/BayesDel + SpliceAI, automated only",
         "default_manual": "ESM1b + OpenSpliceAI with eRepo curated evidence",
         "default": "ESM1b + OpenSpliceAI, automated only"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--config", default="revel_manual")
    ap.add_argument("--min-n", type=int, default=10)
    ap.add_argument("--out", default="fig_gene_concordance")
    ap.add_argument("--no-title", action="store_true", help="omit the title above the plot")
    a = ap.parse_args()
    g = pd.read_csv(a.summary, sep="\t")
    g = g[(g.config == a.config) & (g.n >= a.min_n)].copy()
    g = g.sort_values(["exact %", "exact or confidence %"], ascending=[True, True])
    h = 0.16 * len(g) + 1.6
    fig, ax = plt.subplots(figsize=(7.2, h))
    left = pd.Series(0.0, index=g.index)
    for key, lab, col in CATS:
        ax.barh(range(len(g)), g[key + " %"], left=left, color=col, edgecolor="white",
                linewidth=0.6, height=0.82, label=lab)
        left += g[key + " %"]
    ax.set_yticks(range(len(g)))
    ax.set_yticklabels([f"{r.gene} ({r.n})" for r in g.itertuples()], fontsize=6.5, color=INK)
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.6, len(g) - 0.4)
    ax.set_xlabel("Variants (%)", color=INK)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    for x in (25, 50, 75):
        ax.axvline(x, color="white", linewidth=0.6, zorder=3)
    for s in ("right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.005), ncol=2, frameon=False, fontsize=7.5)
    med = g["exact %"].median()
    if not a.no_title:
        ax.set_title(f"{LABEL.get(a.config, a.config)}\n{len(g)} genes with ≥{a.min_n} variants; "
                     f"median exact {med:.1f}% (gene: number of variants)", fontsize=8.5, color=INK2,
                     loc="left", fontweight="normal", pad=28)
    save(fig, a.out)


if __name__ == "__main__":
    main()
