#!/usr/bin/env python3
"""Reviewer 2, Comment 6: compare HUVar classifications of the eRepo RYR1 / ACTA1 /
VWF variants under each disease-specific CSpec (run_cspec_test.sh) with the
default merged specification and with the eRepo (VCEP) classification.

  python compare_cspec.py --dir <revision_cspec_YYYYMMDD> \
      [--erepo <HUVar>/resources/GRCh38/erepo_classifications_hg38.tsv] \
      [--erepo-report erepo.report.json]
--erepo-report (the eRepo report download, with the "MOI" column) assigns the
Congenital Myopathies VCEP classifications to the AD or AR specification by their
mode of inheritance; without it both myopathy specifications are "matching".
Writes cspec_summary.txt, cspec_per_variant.tsv (classification and met criteria
under every specification) and cspec_changed.tsv (variants whose classification
depends on the specification) into --dir.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

RUNS = {  # run -> {gene: cspec}
    "merged": {},
    "A": {"RYR1": "malignant_hyperthermia", "ACTA1": "congenital_myopathy_ad", "VWF": "vwd_ad"},
    "B": {"RYR1": "congenital_myopathy_ad", "ACTA1": "congenital_myopathy_ar", "VWF": "vwd_ar"},
    "C": {"RYR1": "congenital_myopathy_ar"},
}
SPECS = {"RYR1": ["malignant_hyperthermia", "congenital_myopathy_ad", "congenital_myopathy_ar"],
         "ACTA1": ["congenital_myopathy_ad", "congenital_myopathy_ar"],
         "VWF": ["vwd_ad", "vwd_ar"]}
CODES = ["PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1",
         "PP2", "PP3", "PP4", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4",
         "BP5", "BP7"]
HU = {"Pathogenic": "P", "LikelyPathogenic": "LP", "VUS": "VUS", "LikelyBenign": "LB", "Benign": "B"}
ER = {"Pathogenic": "P", "Likely Pathogenic": "LP", "Uncertain Significance": "VUS",
      "Likely Benign": "LB", "Benign": "B"}
RANK = {"P": 2, "LP": 1, "VUS": 0, "LB": -1, "B": -2}
SHORT = {"Supporting": "_P", "Moderate": "_M", "Strong": "_S", "VeryStrong": "_VS", "StandAlone": ""}


def _cm_specs(moi: str) -> list[str]:
    m = str(moi).lower()
    if "dominant" in m:
        return ["congenital_myopathy_ad"]
    if "recessive" in m:
        return ["congenital_myopathy_ar"]
    return ["congenital_myopathy_ad", "congenital_myopathy_ar"]


def matching_specs(gene: str, panel: str, disease: str, moi: str = "") -> list[str]:
    """Specifications corresponding to the VCEP that classified the variant (the
    Congenital Myopathies VCEP AD/AR specification from the eRepo MOI, if known)."""
    p, d = str(panel), str(disease).lower()
    if gene == "RYR1":
        out = []
        if "Malignant Hyperthermia" in p:
            out.append("malignant_hyperthermia")
        if "Congenital Myopathies" in p:
            out += _cm_specs(moi)
        return out
    if gene == "ACTA1":
        return _cm_specs(moi)
    if gene == "VWF":
        if "type 2" in d or "disease 2" in d:
            return ["vwd_ad"]
        if "type 3" in d:
            return ["vwd_ar"]
        return ["vwd_ad", "vwd_ar"]
    return []


def met(row) -> str:
    out = []
    for c in CODES:
        if str(row.get(c)) == "1":
            s = str(row.get(c + "_strength", ""))
            out.append(c + SHORT.get(s, "_" + s))
    return ",".join(out)


def cat(t: str, p: str) -> str:
    if t == p:
        return "exact"
    a, b = RANK[t], RANK[p]
    return "confidence" if a * b > 0 else ("significant" if a == 0 or b == 0 else "opposite")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--erepo", type=Path,
                    default=Path.home() / "tools/HUVar/resources/GRCh38/erepo_classifications_hg38.tsv")
    ap.add_argument("--erepo-report", type=Path, default=None,
                    help="eRepo report JSON ({'header':…, 'data':[…]}) with the MOI column")
    a = ap.parse_args()
    moi_map: dict[tuple[str, str], str] = {}
    if a.erepo_report:
        import json
        import re
        for x in json.load(open(a.erepo_report, encoding="utf-8"))["data"]:
            t = str(x.get("Preferred Variant Title", ""))
            mg = re.search(r"\(([A-Z0-9]+)\):(c\.[^ ]+)", t)
            if mg and mg.group(1) in SPECS:
                key = (mg.group(1), mg.group(2), str(x.get("Expert Panel", "")))
                moi_map.setdefault(key, str(x.get("MOI", "")))
    er = pd.read_csv(a.erepo, sep="\t", comment="#", dtype=str)
    er = er[er.gene.isin(SPECS)].drop_duplicates("variant_id").set_index("variant_id")
    er["truth"] = er.assertion.map(ER)

    L, rows = [], []
    p = L.append
    cfgs = sorted({f.name.split("_")[1] + ("_manual" if "_manual_" in f.name else "")
                   for f in a.dir.glob("cspec_*_merged.tsv")})
    if not cfgs:
        raise SystemExit(f"no cspec_<config>_merged.tsv in {a.dir} — pass the output folder of "
                         "run_cspec_test.sh (revision_cspec_YYYYMMDD)")
    for cfg in cfgs:
        res = {}
        for run in RUNS:
            f = a.dir / f"cspec_{cfg}_{run}.tsv"
            if not f.exists():
                continue
            d = pd.read_csv(f, sep="\t", dtype=str, comment="#")
            d = d[d.gene.isin(SPECS) & d.variant_id.isin(er.index)].drop_duplicates("variant_id").set_index("variant_id")
            res[run] = d
        if "merged" not in res:
            continue
        m = res["merged"]
        for vid, r in m.iterrows():
            g = r.gene
            cdot = str(r.hgvs_c).split(":")[-1]
            moi = (moi_map.get((g, cdot, "Congenital Myopathies VCEP"), "")
                   if "Congenital Myopathies" in str(er.at[vid, "expert_panel"]) else "")
            rec = {"config": cfg, "variant_id": vid, "gene": g, "hgvs_c": r.hgvs_c, "hgvs_p": r.hgvs_p,
                   "erepo": er.at[vid, "truth"], "expert_panel": er.at[vid, "expert_panel"],
                   "disease": er.at[vid, "disease"], "erepo_codes": er.at[vid, "met_codes"],
                   "moi": moi,
                   "matching_specs": ",".join(matching_specs(g, er.at[vid, "expert_panel"], er.at[vid, "disease"], moi)),
                   "merged": HU.get(r.classification_bayesian, r.classification_bayesian),
                   "merged_score": r.bayesian_score, "merged_criteria": met(r)}
            for run, sel in RUNS.items():
                if run == "merged" or g not in sel or run not in res or vid not in res[run].index:
                    continue
                x = res[run].loc[vid]
                sp = sel[g]
                rec[sp] = HU.get(x.classification_bayesian, x.classification_bayesian)
                rec[sp + "_score"] = x.bayesian_score
                rec[sp + "_criteria"] = met(x)
            rows.append(rec)

        df = pd.DataFrame([r for r in rows if r["config"] == cfg])
        p(f"===== {cfg} =====")
        for g, specs in SPECS.items():
            sub = df[df.gene == g]
            if sub.empty:
                continue
            p(f"-- {g} (n={len(sub)}; eRepo panels: {sub.expert_panel.value_counts().to_dict()})")
            p(f"   {'specification':26s} exact  P/LP-agree  changed_vs_merged  [matching-VCEP subset: n exact]")
            for sp in ["merged"] + specs:
                if sp not in sub:
                    continue
                ok = sub[sp].notna()
                s = sub[ok]
                ex = (s[sp] == s.erepo).mean()
                pl = ((s[sp].isin(["P", "LP"])) == (s.erepo.isin(["P", "LP"]))).mean()
                ch = "" if sp == "merged" else f"{(s[sp] != s.merged).sum():5d}"
                mt = s if sp == "merged" else s[s.matching_specs.str.contains(sp)]
                mex = (mt[sp] == mt.erepo).sum() if len(mt) else 0
                p(f"   {sp:26s} {ex:5.3f}  {pl:10.3f}  {ch:>17s}  [{len(mt)}: {mex}]")
            # best: classification under the VCEP-matching specification vs merged
            sub = sub.copy()
            def under_match(r):
                ms = [m for m in str(r.matching_specs).split(",") if m and isinstance(r.get(m), str)]
                return r[ms[0]] if len(ms) == 1 else None
            sub["matched"] = sub.apply(under_match, axis=1)
            t = sub[sub.matched.notna()]
            if len(t):
                p(f"   single matching spec (n={len(t)}): exact merged {(t.merged == t.erepo).mean():.3f}"
                  f" -> matched {(t.matched == t.erepo).mean():.3f}; categories merged "
                  f"{pd.Series([cat(a_, b_) for a_, b_ in zip(t.erepo, t.merged)]).value_counts().to_dict()}"
                  f" matched {pd.Series([cat(a_, b_) for a_, b_ in zip(t.erepo, t.matched)]).value_counts().to_dict()}")
            # Congenital Myopathies VCEP by mode of inheritance (eRepo MOI)
            if a.erepo_report and g in ("RYR1", "ACTA1"):
                cm = sub[sub.expert_panel.str.contains("Congenital Myopathies", na=False)]
                for moi, s2 in cm.groupby(cm.moi.replace("", "unknown")):
                    cols = [c for c in ["merged"] + specs if c in s2]
                    p(f"   CM VCEP, MOI {moi} (n={len(s2)}): " + ", ".join(
                        f"{c} {(s2[c] == s2.erepo).sum()} ({100 * (s2[c] == s2.erepo).mean():.1f}%)" for c in cols))
        p("")
    out = pd.DataFrame(rows)
    out.to_csv(a.dir / "cspec_per_variant.tsv", sep="\t", index=False)
    spec_cols = [s for g in SPECS for s in SPECS[g]]
    def varies(r):
        v = {r[s] for s in spec_cols if s in r and isinstance(r[s], str)} | {r["merged"]}
        return len(v) > 1
    ch = out[out.apply(varies, axis=1)]
    ch.to_csv(a.dir / "cspec_changed.tsv", sep="\t", index=False)
    p(f"variants whose classification depends on the specification: "
      f"{ch.groupby(['config', 'gene']).size().to_dict()}")
    txt = "\n".join(L)
    (a.dir / "cspec_summary.txt").write_text(txt + "\n", encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
