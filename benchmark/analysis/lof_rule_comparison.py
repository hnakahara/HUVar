#!/usr/bin/env python3
"""Reviewer 2 C5 / Reviewer 3 C2: alternative rules for "LoF is an established
disease mechanism" (PVS1, genes without a VCEP decision), compared against
ClinGen haploinsufficiency (HI) scores and the VCEP PVS1 decisions.

Rules compared
  count>=3          current rule: >=3 ClinVar P/LP null variants
  frac>=X           null variants make up >=X of the gene's ClinVar P/LP variants
                    (X = 5, 10, 15, 20 %; the gene must have at least --min-plp
                    ClinVar P/LP variants in total, default 1, 3, 5, 10)
  LOEUF             gnomAD LOEUF below the gnomAD-recommended cut-off:
                    v2.1.1 < 0.35 (GRCh37 constraint table),
                    v4 < 0.45 (recommended for the v4.1.1 constraint release;
                    also applied to the v4.1 table used by HUVar)
  LOEUF | frac>=X   either
(HI score 3 and VCEP "applicable" are separate routes in HUVar; they are used
here only as references.)

References
  HI:   HI score 3 (positive) vs HI 0 / 1 / 40 (negative); HI 2 / 30 excluded.
        HI 30 (autosomal recessive) genes are reported separately as the share
        of recessive genes each rule recognises.
  VCEP: PVS1 applicable (positive) vs not applicable (negative).
Also: genes recognised by each rule among ClinGen GDV Definitive/Strong/Moderate
genes, by inheritance (gene_inheritance.tsv).

  python lof_rule_comparison.py --data-dir $DATA --out-dir lof_rules \
      [--constraint-v4 gnomad.v4.1.1.constraint_metrics.tsv] [--min-stars 1] [--min-plp 1 3 5 10]
Run in the HUVar environment (imports acmg_classifier).
"""
from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter
from pathlib import Path

from acmg_classifier.local_db.clingen_gene_db import load_gene_validity, load_hi_scores
from acmg_classifier.local_db.clinvar_sqlite import _get_conn, _is_null_change, _PP2_PATH
from acmg_classifier.local_db.gnomad_db import _load_constraint

FRACS = (0.05, 0.10, 0.15, 0.20)
SUPPORTED = {"definitive", "strong", "moderate"}
GDV_ORDER = ["definitive", "strong", "moderate", "limited", "disputed", "refuted",
             "no known disease relationship"]


def best_gdv(classes):
    alias = {"disputing": "disputed", "refuting": "refuted",
             "no reported evidence": "no known disease relationship"}
    norm = [alias.get(c.strip().lower(), c.strip().lower()) for c in classes]
    return next((c for c in GDV_ORDER if c in norm), norm[0] if norm else "")


def clinvar_counts(db: Path, min_stars: int) -> dict[str, tuple[int, int]]:
    """gene -> (P/LP null variants, all P/LP variants)."""
    rows = _get_conn(db).execute(
        "SELECT gene_symbol, hgvs_c, hgvs_p FROM variants "
        "WHERE star_rating >= ? AND clinical_significance IN (?,?,?)",
        (min_stars, *_PP2_PATH)).fetchall()
    tot: Counter = Counter(); nul: Counter = Counter()
    for g, hc, hp in rows:
        if not g:
            continue
        tot[g] += 1
        if _is_null_change(hc, hp):
            nul[g] += 1
    return {g: (nul[g], tot[g]) for g in tot}


def metrics(ref, pos, test):
    tp = sum(1 for r in ref if pos(r) and test(r)); fn = sum(1 for r in ref if pos(r) and not test(r))
    fp = sum(1 for r in ref if not pos(r) and test(r)); tn = sum(1 for r in ref if not pos(r) and not test(r))
    f = lambda x, y: x / (x + y) if x + y else float("nan")
    sens, spec, ppv = f(tp, fn), f(tn, fp), f(tp, fp)
    return sens, spec, ppv, (sens + spec) / 2 if sens == sens and spec == spec else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=Path("lof_rules"))
    ap.add_argument("--constraint-v4", type=Path, default=None,
                    help="gnomAD v4 constraint table (default: data/GRCh38/gnomad/gnomad_v4.1_constraint.tsv)")
    ap.add_argument("--constraint-v2", type=Path, default=None,
                    help="gnomAD v2.1.1 constraint table (default: data/GRCh37/gnomad/gnomad_v2.1.1_constraint.tsv)")
    ap.add_argument("--min-stars", type=int, default=1, help="ClinVar review stars (default 1)")
    ap.add_argument("--min-plp", type=int, nargs="+", default=[1, 3, 5, 10],
                    help="minimum numbers of ClinVar P/LP variants (all types) in the gene for the "
                         "fraction rule (default: 1 3 5 10)")
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    sh = a.data_dir / "shared"

    vcep = {}
    with (sh / "disease_prevalence.tsv").open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = (r.get("pvs1") or "").strip()
            if v in ("applicable", "not_applicable"):
                vcep[r["gene_symbol"].strip()] = v
    hi = load_hi_scores(sh / "ClinGen_gene_curation_list_GRCh38.tsv")
    gdv = {g: best_gdv(c) for g, c in load_gene_validity(sh / "clingen_gene_disease_validity.csv").items()}
    inh = {}
    if (sh / "gene_inheritance.tsv").exists():
        with (sh / "gene_inheritance.tsv").open(encoding="utf-8") as fh:
            inh = {r["gene"].strip(): r["inheritance"].strip() for r in csv.DictReader(fh, delimiter="\t")}
    cv = clinvar_counts(a.data_dir / "GRCh38" / "clinvar" / "clinvar_ps1_pm5_GRCh38.sqlite", a.min_stars)
    c4p = a.constraint_v4 or a.data_dir / "GRCh38" / "gnomad" / "gnomad_v4.1_constraint.tsv"
    c2p = a.constraint_v2 or a.data_dir / "GRCh37" / "gnomad" / "gnomad_v2.1.1_constraint.tsv"
    c4 = _load_constraint(c4p) if c4p.exists() else {}
    c2 = _load_constraint(c2p) if c2p.exists() else {}

    genes = sorted(set(vcep) | set(hi) | set(gdv) | set(cv))
    rows = []
    for g in genes:
        n, t = cv.get(g, (0, 0))
        rows.append({"gene": g, "vcep_pvs1": vcep.get(g, ""), "hi_score": (hi.get(g) or "").strip(),
                     "gdv_best": gdv.get(g, ""), "inheritance": inh.get(g, ""),
                     "clinvar_plp_null": n, "clinvar_plp_all": t,
                     "null_fraction": round(n / t, 4) if t else "",
                     "loeuf_v4": c4.get(g, (None, None, None))[1],
                     "loeuf_v2": c2.get(g, (None, None, None))[1]})
    with (a.out_dir / "lof_rules_per_gene.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    def frac(x, mn):
        return lambda r: r["clinvar_plp_all"] >= mn and r["clinvar_plp_null"] > 0 and \
            r["clinvar_plp_null"] / r["clinvar_plp_all"] >= x
    lo4 = lambda r: r["loeuf_v4"] is not None and r["loeuf_v4"] < 0.45
    lo2 = lambda r: r["loeuf_v2"] is not None and r["loeuf_v2"] < 0.35
    rules = [("count>=3 (current)", lambda r: r["clinvar_plp_null"] >= 3)]
    rules += [("LOEUF_v4<0.45", lo4), ("LOEUF_v2<0.35", lo2)]
    for mn in a.min_plp:
        sfx = f" (P/LP>={mn})"
        rules += [(f"frac>={int(x*100)}%{sfx}", frac(x, mn)) for x in FRACS]
        rules += [(f"LOEUF_v4<0.45 | frac>={int(x*100)}%{sfx}", (lambda f_: lambda r: lo4(r) or f_(r))(frac(x, mn)))
                  for x in FRACS]
        rules += [(f"LOEUF_v2<0.35 | frac>={int(x*100)}%{sfx}", (lambda f_: lambda r: lo2(r) or f_(r))(frac(x, mn)))
                  for x in FRACS]

    ref_hi = [r for r in rows if r["hi_score"] in ("3", "0", "1", "40")]
    ref_vc = [r for r in rows if r["vcep_pvs1"]]
    hi30 = [r for r in rows if r["hi_score"] == "30"]
    dsm = [r for r in rows if r["gdv_best"] in SUPPORTED]
    weak = [r for r in rows if r["gdv_best"] and r["gdv_best"] not in SUPPORTED]

    L = []; p = L.append
    p(f"ClinVar P/LP >= {a.min_stars} star; fraction rules require >= {a.min_plp} P/LP variants in the gene (P/LP>=)")
    p(f"gnomAD v4 constraint: {c4p.name if c4 else 'MISSING'} ({len(c4)} genes); "
      f"v2.1.1: {c2p.name if c2 else 'MISSING'} ({len(c2)} genes)")
    p(f"References: HI 3 = {sum(r['hi_score']=='3' for r in ref_hi)} vs HI 0/1/40 = "
      f"{sum(r['hi_score']!='3' for r in ref_hi)}; VCEP applicable {sum(r['vcep_pvs1']=='applicable' for r in ref_vc)}"
      f" vs not applicable {sum(r['vcep_pvs1']=='not_applicable' for r in ref_vc)}; HI 30 (AR) {len(hi30)};"
      f" GDV D/S/M {len(dsm)}, Limited or weaker {len(weak)}")
    p("")
    hdr = (f"{'rule':48s} | {'HI: sens':>8s} {'spec':>5s} {'PPV':>5s} {'BA':>5s} | {'VCEP: sens':>10s} {'spec':>5s}"
           f" {'PPV':>5s} {'BA':>5s} | {'HI30':>5s} | {'D/S/M':>12s} | {'Lim-':>11s}")
    p(hdr); p("-" * len(hdr))
    out = []
    for name, test in rules:
        h = metrics(ref_hi, lambda r: r["hi_score"] == "3", test)
        v = metrics(ref_vc, lambda r: r["vcep_pvs1"] == "applicable", test)
        a30 = sum(test(r) for r in hi30) / len(hi30) if hi30 else float("nan")
        nd = sum(test(r) for r in dsm); nw = sum(test(r) for r in weak)
        p(f"{name:48s} | {h[0]:8.3f} {h[1]:5.3f} {h[2]:5.3f} {h[3]:5.3f} | {v[0]:10.3f} {v[1]:5.3f} {v[2]:5.3f}"
          f" {v[3]:5.3f} | {a30:5.3f} | {nd:5d} ({100*nd/max(len(dsm),1):4.1f}%) | {nw:4d} ({100*nw/max(len(weak),1):4.1f}%)")
        out.append([name, *h, *v, a30, nd, nw])
    with (a.out_dir / "lof_rules_metrics.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["rule", "hi_sens", "hi_spec", "hi_ppv", "hi_ba", "vcep_sens", "vcep_spec", "vcep_ppv",
                    "vcep_ba", "hi30_recognised", "gdv_dsm_genes", "gdv_limited_or_weaker_genes"])
        w.writerows(out)

    p("")
    p("GDV Definitive/Strong/Moderate genes recognised, by inheritance (%):")
    modes = ("AD", "AR", "AD/AR", "XL")
    p(f"{'rule':48s} " + " ".join(f"{m:>7s}" for m in modes))
    p(f"{'(n genes)':48s} " + " ".join(f"{sum(r['inheritance']==m for r in dsm):7d}" for m in modes))
    for name, test in rules:
        cells = []
        for m in modes:
            sub = [r for r in dsm if r["inheritance"] == m]
            cells.append(f"{100*sum(test(r) for r in sub)/len(sub):6.1f}%" if sub else "     - ")
        p(f"{name:48s} " + " ".join(cells))
    p("")
    for m in modes:
        sub = [r["loeuf_v4"] for r in dsm if r["inheritance"] == m and r["loeuf_v4"] is not None]
        if sub:
            p(f"median LOEUF v4, {m:5s}: {statistics.median(sub):.2f} (n={len(sub)})")
    p("")
    p("Null fraction (ClinVar) of VCEP genes:")
    for lab in ("applicable", "not_applicable"):
        fr = sorted(r["clinvar_plp_null"] / r["clinvar_plp_all"] for r in ref_vc
                    if r["vcep_pvs1"] == lab and r["clinvar_plp_all"])
        if fr:
            q = lambda x: fr[min(len(fr) - 1, int(x * len(fr)))]
            p(f"  PVS1 {lab:15s} n={len(fr):3d}  median {statistics.median(fr):.3f}  "
              f"IQR {q(0.25):.3f}-{q(0.75):.3f}  <5%: {sum(f < 0.05 for f in fr)}  <20%: {sum(f < 0.20 for f in fr)}")
    txt = "\n".join(L)
    (a.out_dir / "lof_rules_summary.txt").write_text(txt + "\n", encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
