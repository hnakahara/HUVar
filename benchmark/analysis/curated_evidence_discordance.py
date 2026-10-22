#!/usr/bin/env python3
"""Contribution of curated evidence to the discordance with the eRepo (Reviewer 3, Comment 18).

For each HUVar configuration, compares the three-class classification (P/LP, VUS, B/LB) of the
automated-only run with that of the run with the eRepo curated-evidence supplement, and counts
  - variants to which the VCEP applied at least one criterion that requires curation,
  - discordant variants (automated only) and how many of them had such a criterion,
  - discordant variants resolved by adding the curated evidence, and newly discordant variants.

Input: the per-variant concordance tables written by erepo_concordance.py
(concordance_<build>_<cfg>_per_variant.tsv and concordance_<build>_<cfg>_manual_per_variant.tsv).

Usage:
  python curated_evidence_discordance.py --dir results/concordance/hg38 --build hg38 \\
      --disease-prevalence <HUVar>/resources/shared/disease_prevalence.tsv
"""
import argparse
import pandas as pd

CURATED = {"PS2", "PS3", "PS4", "PM3", "PM6", "PP1", "PP4", "BS3", "BS4", "BP2", "BP5"}
THREE = {"P": "P/LP", "LP": "P/LP", "VUS": "VUS", "LB": "B/LB", "B": "B/LB"}


def load(d, build, cfg):
    return (pd.read_csv(f"{d}/concordance_{build}_{cfg}_per_variant.tsv", sep="\t", dtype=str)
            .drop_duplicates("variant_id").set_index("variant_id"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--build", default="hg38")
    ap.add_argument("--disease-prevalence", required=True,
                    help="HUVar resources/shared/disease_prevalence.tsv (genes whose VCEP bars gnomAD-based BS2)")
    ap.add_argument("--configs", nargs="+", default=["default", "revel", "spliceai"])
    a = ap.parse_args()
    dp = pd.read_csv(a.disease_prevalence, sep="\t", dtype=str)
    bs2_curated = set(dp.loc[dp.bs2.fillna("") == "not_applicable", "gene_symbol"])
    print("Curated criteria:", ",".join(sorted(CURATED)),
          "+ BS2 in genes whose VCEP bars gnomAD-based BS2 (Bayesian, 3 classes)")
    for cfg in a.configs:
        auto, man = load(a.dir, a.build, cfg), load(a.dir, a.build, cfg + "_manual")
        idx = auto.index.intersection(man.index)
        auto, man = auto.loc[idx], man.loc[idx]
        truth, p_auto, p_man = auto.truth.map(THREE), auto.pred.map(THREE), man.pred.map(THREE)
        has = pd.Series([any((k := c.strip().split("_")[0]) in CURATED or (k == "BS2" and g in bs2_curated)
                             for c in str(s).split(",") if c.strip())
                         for s, g in zip(auto["Applied Evidence Codes (Met)"].fillna(""), auto["Gene"])],
                        index=auto.index)
        disc = p_auto != truth
        resolved = disc & (p_man == truth)
        new = (~disc) & (p_man != truth)
        print(f"{cfg:8} n={len(idx)} with_curated={has.sum()} ({100 * has.mean():.1f}%) "
              f"discordant={disc.sum()} ({100 * disc.mean():.1f}%) "
              f"discordant_with_curated={(disc & has).sum()} ({100 * (disc & has).sum() / disc.sum():.1f}%) "
              f"resolved={resolved.sum()} ({100 * resolved.sum() / disc.sum():.1f}%) newly_discordant={new.sum()}")
        for k in ["P/LP", "VUS", "B/LB"]:
            dk = disc & (truth == k)
            print(f"      truth {k}: discordant={dk.sum()} resolved={(dk & resolved).sum()} "
                  f"({100 * (dk & resolved).sum() / max(dk.sum(), 1):.1f}%)")


if __name__ == "__main__":
    main()
