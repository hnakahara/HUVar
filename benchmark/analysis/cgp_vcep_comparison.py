#!/usr/bin/env python3
"""Reviewer 3, comment 19: HUVar performance in genes with vs without VCEP specifications
(CGP / FoundationOne CDx cohort, ACMG secondary-finding genes, GRCh38).

Truth = the "Expert panel" column of Supplementary Table S3 (P/LP, VUS, B/LB).
Predictions are read either from the HUVar columns of S3 (the submitted version)
or from new HUVar result TSVs (f1_acmg_hg38_<cfg>_results.tsv), matched on
chr:pos:ref:alt. Metrics: P/LP one-vs-rest sensitivity, specificity, precision,
F1, and three-class accuracy (Bayesian classification), per gene group and per
gene.

Usage:
  python cgp_vcep_comparison.py --s3 "Supplementary Tables.xls" \
      --results-dir /path/to/huvar_revision --out cgp_vcep
  # without --results-dir: recompute the submitted numbers from S3 itself
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

VCEP_GENES = {"APC", "BRCA1", "BRCA2", "HNF1A", "MLH1", "MSH2", "MSH6", "PALB2",
              "PMS2", "PTEN", "TP53", "VHL"}
TRUTH3 = {"Pathogenic/Likely_pathogenic": "P/LP", "VUS": "VUS",
          "Benign/Likely_benign": "B/LB"}
PRED3 = {"Pathogenic": "P/LP", "LikelyPathogenic": "P/LP", "VUS": "VUS",
         "UncertainSignificance": "VUS", "LikelyBenign": "B/LB", "Benign": "B/LB"}
CONFIGS = [("ESM1b+OpenSpliceAI", "default"), ("AlphaMissense+SpliceAI", "alphamissense"),
           ("REVEL+SpliceAI", "revel")]


def metrics(t: pd.Series, p: pd.Series) -> dict:
    tp = ((t == "P/LP") & (p == "P/LP")).sum()
    fp = ((t != "P/LP") & (p == "P/LP")).sum()
    fn = ((t == "P/LP") & (p != "P/LP")).sum()
    tn = ((t != "P/LP") & (p != "P/LP")).sum()
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")
    prec = tp / (tp + fp) if tp + fp else float("nan")
    f1 = 2 * prec * sens / (prec + sens) if prec + sens else float("nan")
    return dict(n=len(t), plp_truth=int((t == "P/LP").sum()), sensitivity=round(sens, 3),
                specificity=round(spec, 3), precision=round(prec, 3), f1=round(f1, 3),
                accuracy_3class=round((t == p).mean(), 3))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--s3", type=Path, required=True)
    ap.add_argument("--results-dir", type=Path, default=None)
    ap.add_argument("--out", default="cgp_vcep")
    a = ap.parse_args()

    s3 = pd.read_excel(a.s3, sheet_name="Supplementary Table S3", header=2)
    s3["variant_id"] = (s3.Chromosome.astype(str) + ":" + s3.Location.astype(str) + ":"
                        + s3.Ref.astype(str) + ":" + s3.Alt.astype(str))
    s3["truth"] = s3["Expert panel"].map(TRUTH3)
    s3["group"] = s3.Gene.map(lambda g: "With VCEP" if g in VCEP_GENES else "Without VCEP")

    preds = {}
    for col, cfg in CONFIGS:
        if a.results_dir:
            f = a.results_dir / f"f1_acmg_hg38_{cfg}_results.tsv"
            r = pd.read_csv(f, sep="\t", dtype=str, usecols=["variant_id", "classification_bayesian"])
            r = r.drop_duplicates("variant_id").set_index("variant_id")["classification_bayesian"]
            preds[col] = s3.variant_id.map(r).map(PRED3)
        else:
            preds[col] = s3[col].map(PRED3)

    rows, per_gene = [], []
    for col, _ in CONFIGS:
        s3["pred"] = preds[col]
        ok = s3.dropna(subset=["truth", "pred"])
        missing = len(s3) - len(ok)
        for grp, sub in ok.groupby("group"):
            rows.append(dict(group=grp, config=col, genes=sub.Gene.nunique(), **metrics(sub.truth, sub.pred)))
        rows.append(dict(group="All", config=col, genes=ok.Gene.nunique(), **metrics(ok.truth, ok.pred)))
        for gene, sub in ok.groupby("Gene"):
            per_gene.append(dict(gene=gene, group=sub.group.iloc[0], config=col, **metrics(sub.truth, sub.pred)))
        if missing:
            print(f"[{col}] {missing} variants without a truth or prediction were skipped")

    summ = pd.DataFrame(rows).sort_values(["group", "config"])
    pg = pd.DataFrame(per_gene).sort_values(["config", "group", "gene"])
    summ.to_csv(f"{a.out}_summary.tsv", sep="\t", index=False)
    pg.to_csv(f"{a.out}_per_gene.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))
    low = pg[(pg.config == CONFIGS[0][0]) & (pg.plp_truth >= 4)].sort_values("sensitivity").head(8)
    print("\nLowest P/LP sensitivity (ESM1b+OpenSpliceAI, >=4 P/LP):")
    print(low[["gene", "group", "plp_truth", "sensitivity", "precision"]].to_string(index=False))


if __name__ == "__main__":
    main()
