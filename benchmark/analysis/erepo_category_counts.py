#!/usr/bin/env python3
"""Reviewer 3, Comment 15: exact / confidence / significant / opposite counts and
proportions by configuration, for all eRepo records ("all") and after excluding
the records whose own codes do not reproduce the eRepo class ("excl"), overall
and by eRepo class, with the P/LP F1 score.

  python erepo_category_counts.py --dir <concordance outputs> --build hg19 --out categories_counts.tsv
  (--exclude <list> additionally reports the set without those records)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

CLS = ["P", "LP", "VUS", "LB", "B"]
RANK = {"P": 2, "LP": 1, "VUS": 0, "LB": -1, "B": -2}
CATS = ["exact", "confidence", "significant", "opposite"]
CONFIGS = ["default", "revel", "spliceai", "default_manual", "revel_manual", "spliceai_manual"]


def cat(t: str, p: str) -> str:
    if t == p:
        return "exact"
    a, b = RANK[t], RANK[p]
    return "confidence" if a * b > 0 else ("significant" if a == 0 or b == 0 else "opposite")


def row(d: pd.DataFrame, **k) -> dict:
    c = d._cat.value_counts()
    r = dict(k, n=len(d))
    for x in CATS:
        r[x] = int(c.get(x, 0))
        r[x + "_%"] = round(100 * c.get(x, 0) / len(d), 1) if len(d) else 0.0
    return r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--build", default="hg38", choices=["hg38", "hg19"])
    ap.add_argument("--exclude", type=Path, default=None,
                    help="optional: records to report also as an excluded set (default: all records only)")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    excl = set(pd.read_csv(a.exclude, sep="\t", dtype=str).variant_id) if a.exclude else None
    rows = []
    for cfg in CONFIGS:
        f = a.dir / f"concordance_{a.build}_{cfg}_per_variant.tsv"
        if not f.exists():
            continue
        d0 = pd.read_csv(f, sep="\t", dtype=str)
        d0 = d0[d0.truth.isin(CLS) & d0.pred.isin(CLS)].copy()
        d0["_cat"] = [cat(t, p) for t, p in zip(d0.truth, d0.pred)]
        sets = [("all", d0)] if excl is None else [("excl", d0[~d0.variant_id.isin(excl)]), ("all", d0)]
        for s, d in sets:
            t, p = d.truth.isin(["P", "LP"]), d.pred.isin(["P", "LP"])
            tp, fp, fn = int((t & p).sum()), int((~t & p).sum()), int((t & ~p).sum())
            f1 = round(2 * tp / (2 * tp + fp + fn), 3) if tp else 0.0
            rows.append(dict(row(d, config=cfg, set=s, eRepo="All"), PLP_F1=f1))
            for k in CLS:
                rows.append(row(d[d.truth == k], config=cfg, set=s, eRepo=k))
    out = pd.DataFrame(rows)
    out.to_csv(a.out, sep="\t", index=False)
    print(out[out.eRepo == "All"].to_string(index=False))


if __name__ == "__main__":
    main()
