#!/usr/bin/env python3
"""Reviewer 3, Comment 16: does ClinVar data newer than the VCEP evaluation explain discordance?

Run on the server where the ClinVar SQLite lives.

Usage:
  python erepo_date_analysis.py \
      --variants  diag_revel_manual_variants.tsv \
      --results   erepo_hg38_huhvar_revel_manual_results.tsv \
      --sqlite    $DATA/GRCh38/clinvar/clinvar_ps1_pm5_GRCh38.sqlite \
      --clinvar-release 2026-05-30 \
      --outdir    results/dates

Inputs
  --variants : per-variant concordance from erepo_concordance_diag.py
               (variant_id, Gene, truth, pred, cat, Applied Evidence Codes (Met), ...)
  --results  : HUVar result TSV for the same configuration (criteria flags + *_evidence)
  --sqlite   : HUVar ClinVar SQLite (table `variants`, column `last_evaluated`)

Outputs (in --outdir)
  per_variant_dates.tsv        VCEP evaluation date (ClinVar expert-panel record) per eRepo variant
  by_year.tsv                  concordance and HUVar-only PS1/PM5/PS3/PP1 rates by VCEP evaluation year
  ps1_pm5_comparators.tsv      HUVar-only PS1/PM5 calls with comparator last-evaluated dates
  summary.txt
"""
import argparse, re, sqlite3
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--variants", required=True)
ap.add_argument("--results", required=True)
ap.add_argument("--sqlite", required=True)
ap.add_argument("--clinvar-release", default="2026-05-30")
ap.add_argument("--outdir", default="results/dates")
a = ap.parse_args()
out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

v = pd.read_csv(a.variants, sep="\t", dtype=str)
cols = ["variant_id"] + [f"{c}{s}" for c in ("PS1", "PM5", "PS3", "PP1") for s in ("", "_evidence")]
r = pd.read_csv(a.results, sep="\t", dtype=str, usecols=cols)
m = v.merge(r, on="variant_id", how="left")

# immutable=1: the SQLite lives on a network mount where read-only locking fails
con = sqlite3.connect(f"file:{a.sqlite}?mode=ro&immutable=1", uri=True)

def ep_date(vid):
    chrom, pos, ref, alt = vid.split(":")
    chrom = chrom.replace("chr", "")
    row = con.execute(
        "SELECT MAX(last_evaluated) FROM variants WHERE chrom=? AND pos=? AND ref=? AND alt=? "
        "AND star_rating>=3", (chrom, int(pos), ref, alt)).fetchone()
    return row[0] if row and row[0] else None

m["vcep_eval_date"] = pd.to_datetime([ep_date(x) for x in m.variant_id], errors="coerce")
m["eval_year"] = m["vcep_eval_date"].dt.year
rel = pd.Timestamp(a.clinvar_release)
m["years_before_clinvar"] = ((rel - m["vcep_eval_date"]).dt.days / 365.25).round(2)

def has_code(codes, c):
    return any(re.match(rf"{c}(_|$)", t.strip()) for t in str(codes).split(","))

for c in ("PS1", "PM5", "PS3", "PP1"):
    m[f"{c}_huvar_only"] = (m[c] == "1") & ~m["Applied Evidence Codes (Met)"].map(lambda s, c=c: has_code(s, c))

m.to_csv(out / "per_variant_dates.tsv", sep="\t", index=False)

g = m.dropna(subset=["eval_year"]).groupby("eval_year")
by = pd.DataFrame({
    "n": g.size(),
    "exact_pct": g.apply(lambda x: 100 * (x["cat"] == "exact").mean()).round(1),
    "significant_or_opposite_pct": g.apply(lambda x: 100 * x["cat"].isin(["significant", "opposite"]).mean()).round(1),
    "huvar_more_pathogenic_pct": g.apply(lambda x: 100 * (
        x["pred"].map({"P": 2, "LP": 1, "VUS": 0, "LB": -1, "B": -2}) >
        x["truth"].map({"P": 2, "LP": 1, "VUS": 0, "LB": -1, "B": -2})).mean()).round(1),
    **{f"{c}_huvar_only_pct": g.apply(lambda x, c=c: 100 * x[f"{c}_huvar_only"].mean()).round(1)
       for c in ("PS1", "PM5", "PS3", "PP1")},
})
by.to_csv(out / "by_year.tsv", sep="\t")

# PS1/PM5 comparators: were they last evaluated after the VCEP evaluation?
rows = []
for _, x in m[(m.PS1_huvar_only) | (m.PM5_huvar_only)].iterrows():
    for c in ("PS1", "PM5"):
        if not x[f"{c}_huvar_only"]:
            continue
        ids = re.findall(r"\b(\d{3,9})\b", str(x[f"{c}_evidence"]).split(":", 1)[-1])
        for vid in ids[:3]:
            d = con.execute("SELECT MIN(last_evaluated), MAX(star_rating) FROM variants WHERE variation_id=?",
                            (vid,)).fetchone()
            rows.append(dict(variant_id=x.variant_id, gene=x.Gene, criterion=c, comparator_vcv=vid,
                             comparator_last_eval=d[0] if d else None,
                             comparator_stars=d[1] if d else None, vcep_eval_date=x.vcep_eval_date))
cp = pd.DataFrame(rows)
if len(cp):
    cp["comparator_last_eval"] = pd.to_datetime(cp["comparator_last_eval"], errors="coerce")
    cp["comparator_after_vcep"] = cp["comparator_last_eval"] > cp["vcep_eval_date"]
cp.to_csv(out / "ps1_pm5_comparators.tsv", sep="\t", index=False)

with open(out / "summary.txt", "w") as fh:
    fh.write(f"eRepo variants: {len(m)}; with expert-panel date: {m.vcep_eval_date.notna().sum()}\n")
    fh.write(f"VCEP evaluation date range: {m.vcep_eval_date.min()} - {m.vcep_eval_date.max()}\n")
    fh.write(f"Median years between VCEP evaluation and ClinVar release ({a.clinvar_release}): "
             f"{m.years_before_clinvar.median()}\n")
    if len(cp):
        per_var = cp.groupby(["variant_id", "criterion"])["comparator_after_vcep"].any()
        fh.write(f"HUVar-only PS1/PM5 calls: {len(per_var)}; with >=1 comparator last evaluated after "
                 f"the VCEP evaluation: {per_var.sum()} ({100*per_var.mean():.1f}%)\n")
    fh.write("\nBy VCEP evaluation year:\n" + by.to_string() + "\n")
print(open(out / "summary.txt").read())
