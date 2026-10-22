#!/usr/bin/env python3
"""Reviewer 2, Comment 5 (and Reviewer 1 Major 4 / Reviewer 3 Comment 4):
genes for which HUVar considers loss of function (LoF) an established disease
mechanism for PVS1, by route, and their distribution among ClinGen
gene-disease validity (GDV) genes. Also evaluates the ClinVar null-variant
threshold (1-10) against two references: ClinGen HI score 3, and the VCEP
PVS1 decision (applicable / not applicable).

Route priority (as in HUVar):
  VCEP applies PVS1  >  VCEP says not applicable (LoF NOT established)
  >  ClinGen HI score 3  >  >=3 ClinVar P/LP null variants (>=1 star)  >  not established

  python lof_mechanism_genes.py --data-dir $DATA --out-dir lof_out
Run in the HUVar environment (imports acmg_classifier). GRCh38 ClinVar SQLite.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from acmg_classifier.local_db.clingen_gene_db import load_gene_validity, load_hi_scores
from acmg_classifier.local_db.clinvar_sqlite import _get_conn, _is_null_change, _PP2_PATH
from acmg_classifier.local_db.gnomad_db import _load_constraint

GDV_ORDER = ["definitive", "strong", "moderate", "limited", "disputed", "refuted",
             "no known disease relationship"]
SUPPORTED = {"definitive", "strong", "moderate"}


def best_gdv(classes: list[str]) -> str:
    norm = []
    for c in classes:
        c = c.strip().lower()
        c = {"disputing": "disputed", "refuting": "refuted",
             "no reported evidence": "no known disease relationship"}.get(c, c)
        norm.append(c)
    for c in GDV_ORDER:
        if c in norm:
            return c
    return norm[0] if norm else ""


def null_counts(db: Path, min_stars: int = 1) -> dict[str, int]:
    con = _get_conn(db)
    rows = con.execute(
        "SELECT gene_symbol, hgvs_c, hgvs_p FROM variants "
        "WHERE star_rating >= ? AND clinical_significance IN (?,?,?)",
        (min_stars, *_PP2_PATH)).fetchall()
    out: Counter = Counter()
    for g, hc, hp in rows:
        if g and _is_null_change(hc, hp):
            out[g] += 1
    return dict(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--assembly", default="GRCh38")
    ap.add_argument("--out-dir", type=Path, default=Path("lof_out"))
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    shared = a.data_dir / "shared"

    vcep: dict[str, str] = {}
    with (shared / "disease_prevalence.tsv").open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = (r.get("pvs1") or "").strip()
            if v in ("applicable", "not_applicable"):
                vcep[r["gene_symbol"].strip()] = v
    hi = load_hi_scores(shared / "ClinGen_gene_curation_list_GRCh38.tsv")
    gdv = {g: best_gdv(c) for g, c in
           load_gene_validity(shared / "clingen_gene_disease_validity.csv").items()}
    db = a.data_dir / a.assembly / "clinvar" / f"clinvar_ps1_pm5_{a.assembly}.sqlite"
    nulls = null_counts(db)
    ver = "v4.1" if a.assembly == "GRCh38" else "v2.1.1"
    cons = _load_constraint(a.data_dir / a.assembly / "gnomad" / f"gnomad_{ver}_constraint.tsv")
    inh: dict[str, str] = {}
    ip = shared / "gene_inheritance.tsv"
    if ip.exists():
        with ip.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                inh[r["gene"].strip()] = r["inheritance"].strip()

    genes = sorted(set(vcep) | set(hi) | set(gdv) | set(nulls))
    rows = []
    for g in genes:
        n = nulls.get(g, 0)
        h = (hi.get(g) or "").strip()
        if vcep.get(g) == "applicable":
            route = "VCEP"
        elif vcep.get(g) == "not_applicable":
            route = "VCEP_not_applicable"
        elif h == "3":
            route = "HI3"
        elif n >= 3:
            route = "ClinVar>=3null"
        else:
            route = "not_established"
        rows.append({"gene": g, "route": route, "vcep_pvs1": vcep.get(g, ""),
                     "hi_score": h, "clinvar_plp_null": n, "gdv_best": gdv.get(g, ""),
                     "hi3": h == "3", "clinvar_ge3": n >= 3,
                     "loeuf": cons.get(g, (None, None, None))[1], "inheritance": inh.get(g, "")})
    with (a.out_dir / "lof_mechanism_per_gene.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    L = []
    p = L.append
    est = [r for r in rows if r["route"] in ("VCEP", "HI3", "ClinVar>=3null")]
    p(f"Genes considered (union of VCEP, ClinGen DS, ClinGen GDV, ClinVar P/LP null): {len(rows)}")
    p(f"ClinVar SQLite: {db}")
    p("")
    p("== (A/B) LoF established, by first route that applies ==")
    for k, v in Counter(r["route"] for r in rows).most_common():
        p(f"  {k:22s} {v}")
    p(f"  established total      {len(est)}")
    p(f"  HI3 genes also >=3 ClinVar null: {sum(1 for r in rows if r['hi3'] and r['clinvar_ge3'])}"
      f" / HI3 {sum(1 for r in rows if r['hi3'])}")
    p(f"  VCEP-applicable genes also HI3 or >=3 null: "
      f"{sum(1 for r in rows if r['route']=='VCEP' and (r['hi3'] or r['clinvar_ge3']))}"
      f" / {sum(1 for r in rows if r['route']=='VCEP')}")
    p("")
    p("== (C) Distribution among ClinGen GDV genes (best classification per gene) ==")
    p("  GDV class                       genes  VCEP  HI3  ClinVar>=3  established(%)  VCEP_NA")
    for cls in GDV_ORDER + ["(any supported: D/S/M)", "(all GDV genes)"]:
        if cls == "(any supported: D/S/M)":
            sub = [r for r in rows if r["gdv_best"] in SUPPORTED]
        elif cls == "(all GDV genes)":
            sub = [r for r in rows if r["gdv_best"]]
        else:
            sub = [r for r in rows if r["gdv_best"] == cls]
        if not sub:
            continue
        c = Counter(r["route"] for r in sub)
        e = c["VCEP"] + c["HI3"] + c["ClinVar>=3null"]
        p(f"  {cls:30s} {len(sub):6d} {c['VCEP']:5d} {c['HI3']:4d} {c['ClinVar>=3null']:11d}"
          f"  {e:6d} ({100*e/len(sub):4.1f}%)  {c['VCEP_not_applicable']:5d}")
    noGDV = [r for r in est if not r["gdv_best"]]
    p(f"  established genes without a GDV curation: {len(noGDV)}"
      f" (by route: {dict(Counter(r['route'] for r in noGDV))})")
    p("")
    p("== Threshold of ClinVar P/LP null variants (>=1 star) vs references ==")
    ref_hi = [r for r in rows if r["hi_score"] in ("3", "0", "1", "40")]
    p("  Reference 1: ClinGen HI score 3 (positive) vs HI 0/1/40 (negative; HI 2/30 excluded)")
    p(f"    positives {sum(r['hi_score']=='3' for r in ref_hi)}, negatives {sum(r['hi_score']!='3' for r in ref_hi)}")
    ref_vc = [r for r in rows if r["vcep_pvs1"]]
    p(f"  Reference 2: VCEP PVS1 applicable (positive {sum(r['vcep_pvs1']=='applicable' for r in ref_vc)})"
      f" vs not applicable (negative {sum(r['vcep_pvs1']=='not_applicable' for r in ref_vc)})")
    p("  thr | HI: sens  spec  PPV   | VCEP: sens  spec  PPV")
    for t in range(1, 11):
        def m(ref, pos):
            tp = sum(1 for r in ref if pos(r) and r["clinvar_plp_null"] >= t)
            fn = sum(1 for r in ref if pos(r) and r["clinvar_plp_null"] < t)
            fp = sum(1 for r in ref if not pos(r) and r["clinvar_plp_null"] >= t)
            tn = sum(1 for r in ref if not pos(r) and r["clinvar_plp_null"] < t)
            f = lambda x, y: x / (x + y) if x + y else float("nan")
            return f(tp, fn), f(tn, fp), f(tp, fp)
        s1 = m(ref_hi, lambda r: r["hi_score"] == "3")
        s2 = m(ref_vc, lambda r: r["vcep_pvs1"] == "applicable")
        p(f"  {t:3d} |     {s1[0]:.3f} {s1[1]:.3f} {s1[2]:.3f} |       {s2[0]:.3f} {s2[1]:.3f} {s2[2]:.3f}")
    p("")
    p("== gnomAD LOEUF as an alternative to the ClinVar route ==")
    p("  (genes with a LOEUF value; same references as above)")
    def mm(ref, pos, test):
        tp = sum(1 for r in ref if pos(r) and test(r)); fn = sum(1 for r in ref if pos(r) and not test(r))
        fp = sum(1 for r in ref if not pos(r) and test(r)); tn = sum(1 for r in ref if not pos(r) and not test(r))
        f = lambda x, y: x / (x + y) if x + y else float("nan")
        return f(tp, fn), f(tn, fp), f(tp, fp)
    hasl = lambda r: r["loeuf"] is not None
    rh = [r for r in ref_hi if hasl(r)]; rv = [r for r in ref_vc if hasl(r)]
    p(f"  n with LOEUF: HI reference {len(rh)}/{len(ref_hi)}, VCEP reference {len(rv)}/{len(ref_vc)}")
    p("  rule                    | HI: sens  spec  PPV   | VCEP: sens  spec  PPV")
    rules = [(f"LOEUF<{t}", (lambda t: lambda r: r["loeuf"] < t)(t)) for t in (0.35, 0.6, 0.8, 1.0, 1.2)]
    rules += [("ClinVar>=3", lambda r: r["clinvar_plp_null"] >= 3),
              ("ClinVar>=3 & LOEUF<1", lambda r: r["clinvar_plp_null"] >= 3 and r["loeuf"] < 1),
              ("ClinVar>=3 | LOEUF<0.6", lambda r: r["clinvar_plp_null"] >= 3 or r["loeuf"] < 0.6)]
    for name, test in rules:
        s1 = mm(rh, lambda r: r["hi_score"] == "3", test)
        s2 = mm(rv, lambda r: r["vcep_pvs1"] == "applicable", test)
        p(f"  {name:24s}|     {s1[0]:.3f} {s1[1]:.3f} {s1[2]:.3f} |       {s2[0]:.3f} {s2[1]:.3f} {s2[2]:.3f}")
    p("")
    p("  By inheritance (gene_inheritance.tsv), ClinGen GDV Definitive/Strong/Moderate genes with LOEUF:")
    p("  inheritance  genes  LOEUF<0.6  LOEUF<1  ClinVar>=3  HI3   median LOEUF")
    import statistics
    for mode in ("AD", "AR", "AD/AR", "XL"):
        sub = [r for r in rows if r["gdv_best"] in SUPPORTED and r["inheritance"] == mode and hasl(r)]
        if not sub:
            continue
        n_ = len(sub)
        f = lambda k: f"{k:5d} ({100*k/n_:4.1f}%)"
        p(f"  {mode:11s} {n_:6d}  {f(sum(r['loeuf'] < 0.6 for r in sub))}  {f(sum(r['loeuf'] < 1 for r in sub))}"
          f"  {f(sum(r['clinvar_ge3'] for r in sub))}  {sum(r['hi3'] for r in sub):4d}"
          f"   {statistics.median(r['loeuf'] for r in sub):.2f}")
    vc_ar = [r for r in rows if r["vcep_pvs1"] == "applicable" and r["inheritance"] == "AR" and hasl(r)]
    if vc_ar:
        p(f"  VCEP PVS1-applicable AR genes: {len(vc_ar)}, LOEUF<1: {sum(r['loeuf'] < 1 for r in vc_ar)}, "
          f"ClinVar>=3: {sum(r['clinvar_ge3'] for r in vc_ar)}")
    txt = "\n".join(L)
    (a.out_dir / "lof_mechanism_summary.txt").write_text(txt + "\n", encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
