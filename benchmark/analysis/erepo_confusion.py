#!/usr/bin/env python3
"""Reviewer 3, Comment 15: classification-level agreement with the eRepo.

For each eRepo class (P / LP / VUS / LB / B), the number and proportion of
variants that HUVar classified as P / LP / VUS / LB / B (5 x 5 matrix), the
exact / confidence / significant / opposite breakdown, the same per gene, and a
test of between-gene heterogeneity.

Input: the per-variant files written by erepo_concordance.py
  concordance_hg38_<config>_per_variant.tsv   (truth, pred, Gene, cat)
and optionally diag_<config>_variants.tsv (erepo_codes_cls = the class that
HUVar's Bayesian point rules give to the eRepo's own codes at their strengths).

  python erepo_confusion.py --dir <outputs> [--configs revel revel_manual] \
      [--diag diag_revel_manual_variants.tsv] [--min-gene-n 10] [--out <dir>]
Writes confusion_<config>.tsv, confusion_per_gene.tsv,
confusion_gene_summary.tsv, confusion_summary.txt and
Supplementary_Table_eRepo_confusion.xlsx.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import chi2_contingency

CLS = ["P", "LP", "VUS", "LB", "B"]
CATS = ["exact", "confidence", "significant", "opposite"]
LABEL = {"default": "Default (ESM1b + OpenSpliceAI), automated only",
         "default_manual": "Default + eRepo curated evidence",
         "revel": "VCEP-oriented (REVEL + SpliceAI + BayesDel + CADD), automated only",
         "revel_manual": "VCEP-oriented (REVEL + SpliceAI + BayesDel + CADD) + eRepo curated evidence",
         "spliceai": "AlphaMissense + SpliceAI, automated only",
         "spliceai_manual": "AlphaMissense + SpliceAI + eRepo curated evidence",
         "ceiling": "Reference: eRepo applied codes combined with HUVar Bayesian rules"}
RANK = {"P": 2, "LP": 1, "VUS": 0, "LB": -1, "B": -2}


def cat(t: str, p: str) -> str:
    if t == p:
        return "exact"
    a, b = RANK[t], RANK[p]
    return "confidence" if a * b > 0 else ("significant" if a == 0 or b == 0 else "opposite")


def matrix(d: pd.DataFrame, col: str = "pred") -> pd.DataFrame:
    """Rows = eRepo class, columns = HUVar class; 'n (row %)' plus totals."""
    ct = pd.crosstab(d.truth, d[col]).reindex(index=CLS, columns=CLS, fill_value=0)
    rows = []
    for t in CLS:
        n = ct.loc[t].sum()
        r = {"eRepo": t, "n": int(n)}
        for p in CLS:
            r[f"HUVar {p}"] = int(ct.at[t, p])
            r[f"HUVar {p} %"] = round(100 * ct.at[t, p] / n, 1) if n else None
        sub = d[d.truth == t]
        for c in CATS:
            r[c + " %"] = round(100 * (sub["_cat"] == c).mean(), 1) if len(sub) else None
        rows.append(r)
    tot = {"eRepo": "All", "n": int(ct.values.sum())}
    for p in CLS:
        tot[f"HUVar {p}"] = int(ct[p].sum())
        tot[f"HUVar {p} %"] = round(100 * ct[p].sum() / ct.values.sum(), 1)
    for c in CATS:
        tot[c + " %"] = round(100 * (d["_cat"] == c).mean(), 1)
    rows.append(tot)
    return pd.DataFrame(rows)


def pretty(m: pd.DataFrame) -> str:
    lines = [f"{'eRepo':6s}{'n':>7s}  " + "".join(f"{'HUVar ' + p:>17s}" for p in CLS)
             + "   exact  conf  signif  opp"]
    for _, r in m.iterrows():
        cells = "".join(f"{r[f'HUVar {p}']:>8d} ({r[f'HUVar {p} %']:5.1f}%)" for p in CLS)
        lines.append(f"{r.eRepo:6s}{r.n:>7d}  {cells}   {r['exact %']:5.1f} {r['confidence %']:5.1f}"
                     f"  {r['significant %']:5.1f} {r['opposite %']:4.1f}")
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--configs", nargs="+", default=["revel", "revel_manual"])
    ap.add_argument("--diag", type=Path, default=None,
                    help="diag_<config>_variants.tsv (adds the eRepo-codes reference matrix)")
    ap.add_argument("--min-gene-n", type=int, default=10)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--exclude", type=Path, default=None,
                    help="TSV with a variant_id column (e.g. erepo_codes_reclass_variants.tsv: eRepo "
                         "records whose own codes do not reproduce the assertion) to leave out")
    ap.add_argument("--build", default="hg38", choices=["hg38", "hg19"],
                    help="assembly of the input files (concordance_<build>_*, erepo_<build>_huvar_*)")
    a = ap.parse_args()
    out = a.out or a.dir
    out.mkdir(parents=True, exist_ok=True)

    excl = set()
    if a.exclude is not None:
        excl = set(pd.read_csv(a.exclude, sep="\t", dtype=str).variant_id)
    data = {}
    for c in a.configs:
        f = a.dir / f"concordance_{a.build}_{c}_per_variant.tsv"
        if not f.exists():
            raise SystemExit(f"missing {f} (run erepo_concordance.py first)")
        d = pd.read_csv(f, sep="\t", dtype=str)
        d = d[d.truth.isin(CLS) & d.pred.isin(CLS)].copy()
        d = d[~d.variant_id.isin(excl)]
        d["_cat"] = [cat(t, p) for t, p in zip(d.truth, d.pred)]
        data[c] = d
    if a.diag is not None and a.diag.exists():
        g = pd.read_csv(a.diag, sep="\t", dtype=str).drop_duplicates("variant_id")
        g = g[g.truth.isin(CLS) & g.erepo_codes_cls.isin(CLS) & ~g.variant_id.isin(excl)].copy()
        g["pred"] = g.erepo_codes_cls
        g["_cat"] = [cat(t, p) for t, p in zip(g.truth, g.pred)]
        data["ceiling"] = g

    L = []
    p = L.append
    if excl:
        p(f"Excluded {len(excl)} eRepo records listed in {a.exclude.name}")
        p("")
    mats, per_gene, gene_sum = {}, [], []
    for c, d in data.items():
        m = matrix(d)
        mats[c] = m
        m.to_csv(out / f"confusion_{c}.tsv", sep="\t", index=False)
        p(f"===== {LABEL.get(c, c)} (n = {len(d)}) =====")
        p(pretty(m))
        p("")
        if c == "ceiling":
            continue
        # per gene: matrix in long form and summary
        for gene, s in d.groupby("Gene"):
            ct = pd.crosstab(s.truth, s.pred).reindex(index=CLS, columns=CLS, fill_value=0)
            for t in CLS:
                if ct.loc[t].sum():
                    rec = {"config": c, "gene": gene, "eRepo": t, "n": int(ct.loc[t].sum())}
                    rec.update({f"HUVar {q}": int(ct.at[t, q]) for q in CLS})
                    per_gene.append(rec)
            rec = {"config": c, "gene": gene, "n": len(s)}
            for k in CATS:
                rec[k] = int((s._cat == k).sum())
                rec[k + " %"] = round(100 * rec[k] / len(s), 1)
            rec["exact or confidence %"] = round(100 * (rec["exact"] + rec["confidence"]) / len(s), 1)
            for t in CLS:
                st = s[s.truth == t]
                rec[f"eRepo {t} n"] = len(st)
                rec[f"eRepo {t} exact %"] = round(100 * (st.pred == t).mean(), 1) if len(st) else None
            gene_sum.append(rec)
        gs = pd.DataFrame([r for r in gene_sum if r["config"] == c])
        big = gs[gs.n >= a.min_gene_n]
        tab = big[CATS].to_numpy()
        tab = tab[:, tab.sum(0) > 0]
        chi2, pv, dof, _ = chi2_contingency(tab)
        ex = big[["exact"]].assign(other=big.n - big.exact).to_numpy()
        chi2e, pve, dofe, _ = chi2_contingency(ex)
        q = big["exact %"].quantile([0.25, 0.5, 0.75]).round(1).tolist()
        p(f"  by gene (genes with >= {a.min_gene_n} variants: {len(big)} of {len(gs)}; "
          f"{big.n.sum()} variants): exact % median {q[1]} (IQR {q[0]}-{q[2]}), "
          f"range {big['exact %'].min()}-{big['exact %'].max()}")
        p(f"  genes with exact >= 90%: {(big['exact %'] >= 90).sum()}; exact-or-confidence >= 90%: "
          f"{(big['exact or confidence %'] >= 90).sum()}; significant > 20%: {(big['significant %'] > 20).sum()}")
        p(f"  heterogeneity between genes: exact vs not chi2 = {chi2e:.0f}, df = {dofe}, p = {pve:.1e}; "
          f"4 categories chi2 = {chi2:.0f}, df = {dof}, p = {pv:.1e}")
        lo = big.sort_values("exact %").head(10)
        p("  lowest exact (n >= 30): " + "; ".join(
            f"{r.gene} {r['exact %']}% (n={r.n}; conf {r['confidence %']}, signif {r['significant %']})"
            for _, r in big[big.n >= 30].sort_values("exact %").head(10).iterrows()))
        p("  highest exact (n >= 30): " + "; ".join(
            f"{r.gene} {r['exact %']}% (n={r.n})"
            for _, r in big[big.n >= 30].sort_values("exact %", ascending=False).head(10).iterrows()))
        p("")

    # ceiling cross-check: confidence/significant cells explained by the combining rules
    if "ceiling" in data:
        for c in a.configs:
            d = data[c].merge(data["ceiling"][["variant_id", "pred"]].rename(columns={"pred": "ceil"}),
                              on="variant_id", how="left")
            p(f"== {c}: discordant cells where the eRepo's own codes give the same HUVar class "
              f"(i.e. explained by the combining rules, not by the criteria) ==")
            disc = d[d.truth != d.pred]
            for (t, pr), s in disc.groupby(["truth", "pred"]):
                if len(s) >= 20:
                    p(f"  eRepo {t:3s} -> HUVar {pr:3s}: n = {len(s):5d}; eRepo codes also give {pr}: "
                      f"{(s.ceil == pr).sum():5d} ({100 * (s.ceil == pr).mean():.1f}%)")
            p("")

    txt = "\n".join(L)
    (out / "confusion_summary.txt").write_text(txt + "\n", encoding="utf-8")
    print(txt)
    pg = pd.DataFrame(per_gene)
    gsum = pd.DataFrame(gene_sum)
    pg.to_csv(out / "confusion_per_gene.tsv", sep="\t", index=False)
    gsum.to_csv(out / "confusion_gene_summary.tsv", sep="\t", index=False)

    # xlsx
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    wb = Workbook()
    ws = wb.active
    ws.title = "Overall"
    bold = Font(bold=True)
    diag_fill = PatternFill("solid", fgColor="DDEBF7")
    row = 1
    ws.cell(row, 1, "Classification of eRepo variants by HUVar (GRCh38, Bayesian point system). "
                    "Rows: eRepo (VCEP) classification; columns: HUVar classification, n (row %).").font = bold
    row += 2
    for c, m in mats.items():
        ws.cell(row, 1, LABEL.get(c, c) + f" — n = {int(m.iloc[-1].n)}").font = bold
        row += 1
        hdr = ["eRepo", "n"] + [f"HUVar {q}" for q in CLS] + ["Exact", "Confidence", "Significant", "Opposite"]
        for j, h in enumerate(hdr, 1):
            ws.cell(row, j, h).font = bold
        row += 1
        for _, r in m.iterrows():
            ws.cell(row, 1, r.eRepo)
            ws.cell(row, 2, int(r.n))
            for j, q in enumerate(CLS, 3):
                cell = ws.cell(row, j, f"{r[f'HUVar {q}']} ({r[f'HUVar {q} %']:.1f}%)")
                cell.alignment = Alignment(horizontal="right")
                if r.eRepo == q:
                    cell.fill = diag_fill
            for j, k in enumerate(CATS, 8):
                ws.cell(row, j, f"{r[k + ' %']:.1f}%").alignment = Alignment(horizontal="right")
            row += 1
        row += 1
    ws.cell(row, 1, "Exact = same class; confidence = P vs LP or B vs LB; significant = P/LP or B/LB vs VUS; "
                    "opposite = P/LP vs B/LB. The reference matrix classifies the eRepo's own applied "
                    "codes (at their stated strengths) with HUVar's Bayesian rules, i.e. the upper bound "
                    "attainable when every criterion is reproduced exactly.")
    for col, w in zip("ABCDEFGHIJK", [8, 8, 15, 15, 15, 15, 15, 11, 11, 11, 11]):
        ws.column_dimensions[col].width = w
    for name, df in [("By gene", gsum), ("By gene x eRepo class", pg)]:
        w = wb.create_sheet(name)
        for j, h in enumerate(df.columns, 1):
            w.cell(1, j, h).font = bold
        for i, r in enumerate(df.itertuples(index=False), 2):
            for j, v in enumerate(r, 1):
                w.cell(i, j, None if pd.isna(v) else v)
        w.freeze_panes = "C2"
    wb.save(out / "Supplementary_Table_eRepo_confusion.xlsx")


if __name__ == "__main__":
    main()
