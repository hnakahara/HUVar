#!/usr/bin/env python3
"""Reviewer 3, comment 16: is the VCEP-evaluation-year effect explained by gene/VCEP?

Input: per_variant_dates.tsv written by erepo_date_analysis.py.
Older curations (2016-2019) come from a few early VCEPs, so a crude year trend in
HUVar-only PM5/PS1 calls or "HUVar more pathogenic" could reflect which genes were
curated early rather than how much ClinVar evidence accrued since. This script

1. lists the gene composition of each period (early <=2019, middle 2020-2023,
   recent >=2024);
2. compares early vs recent within genes curated in both periods (Mantel-Haenszel
   odds ratio stratified by gene, and recent rates directly standardised to the
   early gene mix);
3. fits logistic models with gene fixed effects (outcome ~ years since the VCEP
   evaluation, i.e. the time ClinVar had to accumulate new evidence, + gene);
4. shows which genes contribute the early-period excess.

Usage:
  python erepo_date_stratified.py per_variant_dates.tsv --out date_stratified
"""
from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import StratifiedTable

RANK = {"P": 2, "LP": 1, "VUS": 0, "LB": -1, "B": -2}
OUTCOMES = ["PM5_huvar_only", "PS1_huvar_only", "PS3_huvar_only", "PP1_huvar_only",
            "huvar_more_pathogenic", "exact", "significant_or_opposite"]


def period(y: float) -> str:
    return "early_le2019" if y <= 2019 else ("middle_2020_2023" if y <= 2023 else "recent_ge2024")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("tsv")
    ap.add_argument("--out", default="date_stratified")
    ap.add_argument("--min-per-period", type=int, default=5,
                    help="genes need this many variants in both early and recent periods")
    a = ap.parse_args()

    d = pd.read_csv(a.tsv, sep="\t", dtype={"Gene": str})
    d = d.dropna(subset=["eval_year"]).copy()
    for c in ("PM5_huvar_only", "PS1_huvar_only", "PS3_huvar_only", "PP1_huvar_only"):
        d[c] = d[c].astype(str).str.lower().eq("true").astype(int)
    d["huvar_more_pathogenic"] = (d.pred.map(RANK) > d.truth.map(RANK)).astype(int)
    d["exact"] = d.cat.eq("exact").astype(int)
    d["significant_or_opposite"] = d.cat.isin(["significant", "opposite"]).astype(int)
    d["period"] = d.eval_year.map(period)
    out = []

    # 1. gene composition
    comp = (d.groupby(["period", "Gene"]).size().rename("n").reset_index()
              .sort_values(["period", "n"], ascending=[True, False]))
    comp["pct_of_period"] = (100 * comp.n / comp.groupby("period").n.transform("sum")).round(1)
    comp.to_csv(f"{a.out}_gene_composition.tsv", sep="\t", index=False)
    out.append("## 1. Gene composition by period (top 10 genes)")
    for p, sub in comp.groupby("period"):
        top = ", ".join(f"{r.Gene} {r.n} ({r.pct_of_period}%)" for r in sub.head(10).itertuples())
        out.append(f"{p}: n={sub.n.sum()}, genes={len(sub)}; {top}")
    out.append("")

    # crude rates by period
    crude = d.groupby("period")[OUTCOMES].mean().mul(100).round(1)
    crude.insert(0, "n", d.groupby("period").size())
    out.append("## Crude rates by period (%)")
    out.append(crude.to_string())
    out.append("")

    # 2. within-gene early vs recent
    er = d[d.period != "middle_2020_2023"]
    cnt = er.groupby(["Gene", "period"]).size().unstack(fill_value=0)
    both = cnt[(cnt.get("early_le2019", 0) >= a.min_per_period)
               & (cnt.get("recent_ge2024", 0) >= a.min_per_period)].index
    w = er[er.Gene.isin(both)]
    out.append(f"## 2. Early (<=2019) vs recent (>=2024) within genes curated in both periods "
               f"(>= {a.min_per_period} variants each): {len(both)} genes, "
               f"early n={int((w.period == 'early_le2019').sum())}, "
               f"recent n={int((w.period == 'recent_ge2024').sum())}")
    out.append("genes: " + ", ".join(f"{g} ({cnt.loc[g, 'early_le2019']}/{cnt.loc[g, 'recent_ge2024']})"
                                     for g in both))
    early_mix = w[w.period == "early_le2019"].Gene.value_counts(normalize=True)
    rows = []
    for o in OUTCOMES:
        e = w[w.period == "early_le2019"][o].mean()
        r = w[w.period == "recent_ge2024"][o].mean()
        rg = w[w.period == "recent_ge2024"].groupby("Gene")[o].mean()
        r_std = float((rg.reindex(early_mix.index) * early_mix).sum())
        tables = []
        for g, sub in w.groupby("Gene"):
            x = sub.period.eq("early_le2019")
            t = np.array([[((sub[o] == 1) & x).sum(), ((sub[o] == 0) & x).sum()],
                          [((sub[o] == 1) & ~x).sum(), ((sub[o] == 0) & ~x).sum()]], dtype=float)
            if t[:, 0].sum() > 0 and t[:, 1].sum() > 0:  # informative strata only
                tables.append(t)
        if tables:
            st = StratifiedTable(tables)
            mh, (lo, hi), p = st.oddsratio_pooled, st.oddsratio_pooled_confint(), st.test_null_odds().pvalue
        else:
            mh = lo = hi = p = float("nan")
        rows.append(dict(outcome=o, early_pct=round(100 * e, 1), recent_pct=round(100 * r, 1),
                         recent_std_to_early_mix_pct=round(100 * r_std, 1),
                         MH_OR_early_vs_recent=round(mh, 2), ci_low=round(lo, 2), ci_high=round(hi, 2),
                         p=float(f"{p:.3g}"), informative_genes=len(tables)))
    wt = pd.DataFrame(rows)
    wt.to_csv(f"{a.out}_within_gene.tsv", sep="\t", index=False)
    out.append(wt.to_string(index=False))
    out.append("")

    # 3. conditional logistic: per year since VCEP evaluation, stratified by gene
    out.append("## 3. Logistic regression with gene fixed effects, all periods")
    out.append("OR per additional year between the VCEP evaluation and the ClinVar release "
               "(crude logistic OR for comparison)")
    import statsmodels.api as sm
    mrows = []
    for o in OUTCOMES:
        sub = d[d.groupby("Gene")[o].transform(lambda s: 0 < s.sum() < len(s))]
        X = pd.get_dummies(sub.Gene, drop_first=True, dtype=float)
        X.insert(0, "years", sub.years_before_clinvar.astype(float))
        X = sm.add_constant(X)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cl = sm.GLM(sub[o].values, X, family=sm.families.Binomial()).fit()
            cr = sm.Logit(d[o].values, sm.add_constant(d[["years_before_clinvar"]].astype(float))).fit(disp=0)
        b, se, pv = cl.params["years"], cl.bse["years"], cl.pvalues["years"]
        mrows.append(dict(outcome=o, crude_OR_per_year=round(np.exp(cr.params.iloc[1]), 3),
                          crude_p=float(f"{cr.pvalues.iloc[1]:.3g}"),
                          gene_adj_OR_per_year=round(np.exp(b), 3),
                          ci_low=round(np.exp(b - 1.96 * se), 3), ci_high=round(np.exp(b + 1.96 * se), 3),
                          p=float(f"{pv:.3g}"), n=len(sub), genes=sub.Gene.nunique()))
    mt = pd.DataFrame(mrows)
    mt.to_csv(f"{a.out}_gene_fixed_effects.tsv", sep="\t", index=False)
    out.append(mt.to_string(index=False))
    out.append("(n/genes: variants in genes where the outcome varies; other genes carry no "
               "within-gene information)")
    out.append("")

    # 4. which genes carry the early excess
    out.append("## 4. Early-period (<=2019) HUVar-only PM5 and HUVar-more-pathogenic by gene "
               "(genes with >=10 early variants)")
    eg = d[d.period == "early_le2019"].groupby("Gene").agg(
        n=("exact", "size"), PM5_only=("PM5_huvar_only", "sum"),
        more_path=("huvar_more_pathogenic", "sum"))
    rg = d[d.period == "recent_ge2024"].groupby("Gene").agg(
        n_recent=("exact", "size"), PM5_only_recent_pct=("PM5_huvar_only", "mean"),
        more_path_recent_pct=("huvar_more_pathogenic", "mean"))
    eg = eg[eg.n >= 10].join(rg)
    eg["PM5_only_pct"] = (100 * eg.PM5_only / eg.n).round(1)
    eg["more_path_pct"] = (100 * eg.more_path / eg.n).round(1)
    eg["PM5_only_recent_pct"] = (100 * eg.PM5_only_recent_pct).round(1)
    eg["more_path_recent_pct"] = (100 * eg.more_path_recent_pct).round(1)
    eg = eg.sort_values("PM5_only", ascending=False)
    eg.to_csv(f"{a.out}_early_genes.tsv", sep="\t")
    out.append(eg[["n", "PM5_only", "PM5_only_pct", "more_path", "more_path_pct", "n_recent",
                   "PM5_only_recent_pct", "more_path_recent_pct"]].to_string())

    txt = "\n".join(out) + "\n"
    with open(f"{a.out}_summary.txt", "w", encoding="utf-8") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
