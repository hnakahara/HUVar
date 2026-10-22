#!/usr/bin/env python3
"""Reviewer 3, comment 11: BS2 performance on eRepo variants.

1. Actual BS2 calls of a HUVar run (automated-only) vs eRepo BS2, by route,
   and the number of eRepo P/LP variants that receive BS2.
2. Simulation of the non-VCEP dominant carrier route on eRepo variants in
   AD / AD-AR genes (gene_inheritance.tsv), as if no VCEP rule existed:
   gnomAD AC >= 3 (approximate carrier count) with and without LOEUF < 0.5.

Usage:
  python bs2_erepo_eval.py --results erepo_hg38_huhvar_default_results.tsv \
     --erepo HUVar/resources/GRCh38/erepo_classifications_hg38.tsv \
     --dp HUVar/resources/shared/disease_prevalence.tsv \
     --inh HUVar/resources/shared/gene_inheritance.tsv --out bs2_eval
"""
import argparse
import pandas as pd


def route(met, e):
    if not met:
        return "-"
    for k, v in (("recessive", "AR homozygotes"), ("X-linked", "XL hemizygotes"),
                 ("incomplete penetrance", "AD homozygotes (VCEP hom-only)"),
                 ("dominant", "AD carriers"), ("expert-panel", "expert panel")):
        if k in e:
            return v
    return "homozygotes/hemizygotes (mode unknown)"


def main():
    ap = argparse.ArgumentParser()
    for a in ("--results", "--erepo", "--dp", "--inh"):
        ap.add_argument(a, required=True)
    ap.add_argument("--out", default="bs2_eval")
    a = ap.parse_args()
    n = pd.read_csv(a.results, sep="\t", dtype=str)
    t = pd.read_csv(a.erepo, sep="\t", comment="#", dtype=str)
    t["bs2_truth"] = t.met_codes.fillna("").str.split(",").apply(
        lambda l: any(c.strip().split("_")[0] == "BS2" for c in l))
    t["cls"] = t.assertion.map(lambda x: "P/LP" if "Pathogenic" in str(x)
                               else ("B/LB" if "Benign" in str(x) else "VUS"))
    d = n.merge(t[["variant_id", "bs2_truth", "cls"]], on="variant_id")
    dp = pd.read_csv(a.dp, sep="\t", dtype=str)
    d["vcep_bs2"] = d.gene.map(dict(zip(dp.gene_symbol, dp.bs2.fillna("")))).fillna("no VCEP")
    d["met"] = d.BS2 == "1"
    d["route"] = [route(m, e) for m, e in zip(d.met, d.BS2_evidence.fillna(""))]
    out = []
    tp = int((d.met & d.bs2_truth).sum())
    out.append(f"variants {len(d)}; eRepo BS2 {int(d.bs2_truth.sum())}; HUVar BS2 {int(d.met.sum())}; "
               f"agree {tp}; precision {tp / max(d.met.sum(), 1):.3f}; recall {tp / max(d.bs2_truth.sum(), 1):.3f}")
    out.append("HUVar BS2 calls by eRepo class: " + str(d[d.met].cls.value_counts().to_dict()))
    out.append("eRepo BS2 by VCEP BS2 status: " + str(d[d.bs2_truth].vcep_bs2.value_counts().to_dict()))
    out.append("variants in genes without VCEP rules: " + str(int((d.vcep_bs2 == "no VCEP").sum())))
    g = d[d.met].groupby("route").agg(n=("met", "size"), eRepo_BS2=("bs2_truth", "sum"),
                                      PLP=("cls", lambda s: (s == "P/LP").sum()),
                                      BLB=("cls", lambda s: (s == "B/LB").sum()))
    out += ["", "## Actual calls by route", g.to_string()]

    inh = pd.read_csv(a.inh, sep="\t", dtype=str)
    d["cgd"] = d.gene.map(dict(zip(inh.iloc[:, 0], inh.iloc[:, 1]))).fillna("")
    d["ac"] = pd.to_numeric(d.gnomad_ac, errors="coerce").fillna(0)
    d["loeuf"] = pd.to_numeric(d.gnomad_loeuf, errors="coerce")
    ad = d[d.cgd.str.contains("AD")]
    out += ["", f"## Simulation: dominant carrier route on AD/AD-AR genes "
                f"({ad.gene.nunique()} genes, {len(ad)} variants), ignoring VCEP rules"]
    rows = []
    for name, m in (("carriers >= 3", ad.ac >= 3), ("carriers >= 3 and LOEUF < 0.5", (ad.ac >= 3) & (ad.loeuf < 0.5))):
        s = ad[m]
        rows.append(dict(rule=name, n=len(s), PLP=int((s.cls == "P/LP").sum()), VUS=int((s.cls == "VUS").sum()),
                         BLB=int((s.cls == "B/LB").sum()), PLP_pct=round(100 * (s.cls == "P/LP").mean(), 1),
                         eRepo_BS2=int(s.bs2_truth.sum()), genes=s.gene.nunique()))
        out.append(f"{name}: P/LP genes {s[s.cls == 'P/LP'].gene.value_counts().head(8).to_dict()}")
    out.append(pd.DataFrame(rows).to_string(index=False))
    txt = "\n".join(out) + "\n"
    open(f"{a.out}_summary.txt", "w", encoding="utf-8").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
