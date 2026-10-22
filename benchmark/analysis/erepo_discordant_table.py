#!/usr/bin/env python3
"""Reviewer 3, Comment 15: variant-level list of significant conflicts
(P/LP or B/LB vs VUS) and opposite classifications (P/LP vs B/LB) between HUVar
and the eRepo, for the supplementary table.

  python erepo_discordant_table.py --dir <eRepo outputs (concordance_hg38_<cfg>_per_variant.tsv)> \
      --results-dir <folder with erepo_hg38_huhvar_<cfg>_results.tsv> \
      --erepo erepo_classifications_hg38.tsv \
      --exclude erepo_codes_reclass_variants.tsv --out <dir> [--configs revel_manual revel]
Writes erepo_discordant_<cfg>.tsv and Supplementary_Table_eRepo_discordant.xlsx
(one sheet per configuration + summary).
"""
from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

import pandas as pd

CLS = ["P", "LP", "VUS", "LB", "B"]
RANK = {"P": 2, "LP": 1, "VUS": 0, "LB": -1, "B": -2}
CODES = ["PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2",
         "PP3", "PP4", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP7"]
SH = {"VeryStrong": "VS", "Strong": "S", "ThreePoint": "3pt", "Moderate": "M", "Supporting": "P", "StandAlone": "SA"}
PTS = {"VeryStrong": 8, "Strong": 4, "ThreePoint": 3, "Moderate": 2, "Supporting": 1, "StandAlone": 8}
# eRepo suffix -> HUVar strength name; unsuffixed = ACMG 2015 default
ER_ST = {"": None, "very": "VeryStrong", "verystrong": "VeryStrong", "very strong": "VeryStrong",
         "very_strong": "VeryStrong", "strong": "Strong", "moderate": "Moderate", "supporting": "Supporting",
         "stand alone": "StandAlone", "standalone": "StandAlone", "stand": "StandAlone", "stand_alone": "StandAlone"}
DEFAULT = {"PVS": "VeryStrong", "PS": "Strong", "PM": "Moderate", "PP": "Supporting", "BA": "StandAlone",
           "BS": "Strong", "BP": "Supporting"}
LABEL = {"revel_manual": "REVEL+SpliceAI+BayesDel+CADD + eRepo curated evidence",
         "revel": "REVEL+SpliceAI+BayesDel+CADD, automated only",
         "default_manual": "Default (ESM1b+OpenSpliceAI) + eRepo curated evidence",
         "default": "Default (ESM1b+OpenSpliceAI), automated only",
         "spliceai_manual": "AlphaMissense+SpliceAI + eRepo curated evidence",
         "spliceai": "AlphaMissense+SpliceAI, automated only"}


def cat(t: str, p: str) -> str:
    if t == p:
        return "exact"
    a, b = RANK[t], RANK[p]
    return "confidence" if a * b > 0 else ("significant" if a == 0 or b == 0 else "opposite")


def er_codes(s: str) -> dict:
    out = {}
    for t in str(s).split(","):
        t = t.strip()
        m = re.match(r"([A-Z]+\d)", t)
        if not m or t == "nan":
            continue
        b = m.group(1)
        sfx = t[len(b):].lstrip("_ ").lower()
        st = ER_ST.get(sfx, sfx)
        out[b] = st or DEFAULT.get(re.match(r"[A-Z]+", b).group(0), "")
    return out


def lab(met: dict) -> str:
    return ",".join(f"{c}_{SH.get(s, s)}({'-' if c.startswith('B') else '+'}{PTS.get(s, 0)})" for c, s in met.items())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--erepo", type=Path, required=True)
    ap.add_argument("--exclude", type=Path, default=None)
    ap.add_argument("--configs", nargs="+", default=["revel_manual", "revel"])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--build", default="hg38", choices=["hg38", "hg19"],
                    help="assembly of the input files (concordance_<build>_*, erepo_<build>_huvar_*)")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    excl = set(pd.read_csv(a.exclude, sep="\t", dtype=str).variant_id) if a.exclude else set()
    er = pd.read_csv(a.erepo, sep="\t", comment="#", dtype=str).drop_duplicates("variant_id").set_index("variant_id")

    sheets, summ = {}, []
    for cfg in a.configs:
        pv = pd.read_csv(a.dir / f"concordance_{a.build}_{cfg}_per_variant.tsv", sep="\t", dtype=str)
        pv = pv[pv.truth.isin(CLS) & pv.pred.isin(CLS) & ~pv.variant_id.isin(excl)]  # one row per eRepo record
        pv["category"] = [cat(t, p) for t, p in zip(pv.truth, pv.pred)]
        n_all = len(pv)
        x = pv[pv.category.isin(["significant", "opposite"])].copy()
        res = pd.read_csv(a.results_dir / f"erepo_{a.build}_huhvar_{cfg}_results.tsv", sep="\t", dtype=str, comment="#",
                          usecols=lambda c: c in ["variant_id", "hgvs_p", "bayesian_score"] + CODES
                          + [c_ + "_strength" for c_ in CODES])
        res = res.drop_duplicates("variant_id").set_index("variant_id")
        rows = []
        for r in x.itertuples(index=False):
            vid = r.variant_id
            ecodes = getattr(r, "_8")
            erd = er_codes(ecodes)
            met = {}
            if vid in res.index:
                rr = res.loc[vid]
                met = {c: str(rr.get(c + "_strength", "")) for c in CODES if str(rr.get(c)) == "1"}
            miss = [c for c in erd if c not in met]
            extra = [c for c in met if c not in erd]
            stdiff = [f"{c} {SH.get(erd[c], erd[c])}->{SH.get(met[c], met[c])}" for c in erd
                      if c in met and erd[c] != met[c] and c != "BA1"]
            direction = "HUVar more pathogenic" if RANK[r.pred] > RANK[r.truth] else "HUVar more benign"
            rows.append({
                "variant_id": vid, "gene": r.Gene, "erepo_record": r.key, "hgvs_c": er.hgvs_c.get(vid, ""),
                "hgvs_p": res.hgvs_p.get(vid, "") if "hgvs_p" in res else "",
                "expert_panel": er.expert_panel.get(vid, ""), "disease": er.disease.get(vid, ""),
                "approval_date": er.approval_date.get(vid, ""),
                "category": r.category, "eRepo": r.truth, "HUVar": r.pred, "direction": direction,
                "eRepo_codes": ecodes, "HUVar_score": res.bayesian_score.get(vid, "") if len(res) else "",
                "HUVar_criteria": lab(met), "eRepo_code_not_in_HUVar": ",".join(miss),
                "HUVar_code_not_in_eRepo": ",".join(extra), "strength_differs": "; ".join(stdiff),
                "erepo_url": er.erepo_url.get(vid, "") if "erepo_url" in er else ""})
        d = pd.DataFrame(rows)
        order = {"opposite": 0, "significant": 1}
        d = d.sort_values(["category", "eRepo", "HUVar", "gene"], key=lambda s: s.map(order) if s.name == "category"
                          else s.map(RANK).fillna(0) * -1 if s.name in ("eRepo", "HUVar") else s)
        d.to_csv(a.out / f"erepo_discordant_{cfg}.tsv", sep="\t", index=False)
        sheets[cfg] = d
        # summary
        cell = d.groupby(["category", "eRepo", "HUVar"]).size()
        miss = Counter(c for v in d.eRepo_code_not_in_HUVar for c in v.split(",") if c)
        extra = Counter(c for v in d.HUVar_code_not_in_eRepo for c in v.split(",") if c)
        summ.append((cfg, n_all, d, cell, miss, extra))

    L = []
    for cfg, n_all, d, cell, miss, extra in summ:
        L.append(f"===== {LABEL.get(cfg, cfg)} (n = {n_all}) =====")
        L.append(f"significant {int((d.category == 'significant').sum())}, opposite {int((d.category == 'opposite').sum())}; "
                 f"{d.direction.value_counts().to_dict()}")
        L.append(cell.to_string())
        L.append(f"eRepo code absent in HUVar: {dict(miss.most_common(12))}")
        L.append(f"HUVar code absent in eRepo: {dict(extra.most_common(12))}")
        L.append("by gene (top 15): " + str(d.gene.value_counts().head(15).to_dict()))
        L.append("")
    txt = "\n".join(L)
    (a.out / "erepo_discordant_summary.txt").write_text(txt + "\n", encoding="utf-8")
    print(txt)

    from openpyxl import Workbook
    from openpyxl.styles import Font
    bold = Font(bold=True)
    wb = Workbook()
    ws = wb.active
    ws.title = "Summary"
    ws.append(["Significant conflicts (P/LP or B/LB vs VUS) and opposite classifications (P/LP vs B/LB) between "
               "HUVar and the eRepo (GRCh38, Bayesian)."])
    ws["A1"].font = bold
    if excl:
        ws.append([f"{len(excl)} eRepo records whose applied codes do not reproduce the eRepo assertion were "
                   "excluded (listed separately)."])
    ws.append([])
    for cfg, n_all, d, cell, miss, extra in summ:
        ws.append([LABEL.get(cfg, cfg), f"n = {n_all}"])
        ws.cell(ws.max_row, 1).font = bold
        ws.append(["Category", "eRepo", "HUVar", "n", "% of all"])
        for c in ws[ws.max_row]:
            c.font = bold
        for (k, t, p), v in cell.items():
            ws.append([k, t, p, int(v), round(100 * v / n_all, 2)])
        ws.append(["Total", "", "", len(d), round(100 * len(d) / n_all, 1)])
        ws.append(["eRepo code absent in HUVar (top)", "; ".join(f"{k} {v}" for k, v in miss.most_common(10))])
        ws.append(["HUVar code absent in eRepo (top)", "; ".join(f"{k} {v}" for k, v in extra.most_common(10))])
        ws.append([])
    ws.column_dimensions["A"].width = 52
    ws.column_dimensions["B"].width = 14
    hdr = {"variant_id": "variant_id (GRCh38)", "gene": "gene", "erepo_record": "eRepo record (gene:HGVS)", "hgvs_c": "HGVS c.", "hgvs_p": "HGVS p.",
           "expert_panel": "expert panel", "disease": "disease", "approval_date": "eRepo approval date",
           "category": "category", "eRepo": "eRepo class", "HUVar": "HUVar class", "direction": "direction",
           "eRepo_codes": "eRepo applied codes", "HUVar_score": "HUVar Bayesian score",
           "HUVar_criteria": "HUVar criteria (points)", "eRepo_code_not_in_HUVar": "eRepo code not in HUVar",
           "HUVar_code_not_in_eRepo": "HUVar code not in eRepo",
           "strength_differs": "strength differs (eRepo->HUVar)", "erepo_url": "eRepo URL"}
    for cfg, d in sheets.items():
        w = wb.create_sheet(("Curated" if cfg.endswith("_manual") else "Automated") + f" ({cfg})"[:31])
        w.append([hdr[c] for c in d.columns])
        for c in w[1]:
            c.font = bold
        for r in d.itertuples(index=False):
            w.append([None if pd.isna(v) else v for v in r])
        w.freeze_panes = "C2"
    wb.save(a.out / "Supplementary_Table_eRepo_discordant.xlsx")


if __name__ == "__main__":
    main()
