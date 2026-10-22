#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_parsing.py  (VCF-based version)

Uses the TransVar-converted VCF as the bridge to the truth coordinates and checks that each tool can be matched by coordinate.
No performance metrics are computed; only the soundness of parsing and matching is checked.

Chain:
  original (GENE:HGVSc -> Assertion)
     └─ linked 100% through NAME=GENE:HGVSc in the VCF -> the assertion is attached to the VCF coordinate (= truth)
  each tool
     ├─ BIAS / HUVar : chr:pos:ref:alt matched directly to the truth coordinate
     └─ InterVar      : ANNOVAR coordinates reproduced deterministically from the VCF and matched

Usage:
    python3 verify_parsing.py --input <uploads> --outdir <out>
"""

import argparse
import csv
import os
from collections import Counter

import erepo_parsers as P


def build_methods(reg, build):
    paths = reg["builds"][build]
    return [
        ("InterVar", lambda: P.load_intervar(paths["intervar"], build)),
        ("BIAS-2015", lambda: P.load_bias(paths["bias"], build)),
        ("BIAS-2015 (no PP5/BP6)",
         lambda: P.load_bias(paths["bias_nopp5bp6"], build, "bias_nopp5bp6")),
        ("HUHVar default (2015)",
         lambda: P.load_huhvar(paths["huhvar_default"], build, "huhvar_default", "2015")),
        ("HUHVar default (bayesian)",
         lambda: P.load_huhvar(paths["huhvar_default"], build, "huhvar_default", "bayesian")),
        ("HUHVar REVEL (2015)",
         lambda: P.load_huhvar(paths["huhvar_revel"], build, "huhvar_revel", "2015")),
        ("HUHVar REVEL (bayesian)",
         lambda: P.load_huhvar(paths["huhvar_revel"], build, "huhvar_revel", "bayesian")),
        ("HUHVar spliceAI (2015)",
         lambda: P.load_huhvar(paths["huhvar_spliceai"], build, "huhvar_spliceai", "2015")),
        ("HUHVar spliceAI (bayesian)",
         lambda: P.load_huhvar(paths["huhvar_spliceai"], build, "huhvar_spliceai", "bayesian")),
    ]


def pct(n, d):
    return f"{(100.0 * n / d):.1f}%" if d else "n/a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="inputs")
    ap.add_argument("--outdir", default="outputs")
    ap.add_argument("--samples", type=int, default=20)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    reg = P.file_registry(args.input)
    orig = P.load_original(reg["original"])

    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 78)
    out("Reference file (erepo-original.csv)")
    out("=" * 78)
    orec = orig["records"]
    cls_c = Counter(r["cls"] for r in orec)
    out(f"Total rows: {len(orec)}")
    out("Assertion -> 5 categories: " +
        ", ".join(f"{k}={cls_c.get(k,0)}" for k in P.CANON_ORDER) +
        (f", [unrecognised={cls_c.get(None,0)}]" if cls_c.get(None) else ""))
    out(f"Unique GENE:HGVSc: "
        f"{len(set((r['gene'], r['hgvsc_raw']) for r in orec if r['hgvsc_raw']))}")
    out("")

    unmatched_path = os.path.join(args.outdir, "qc_unmatched_samples.tsv")
    uf = open(unmatched_path, "w", newline="", encoding="utf-8")
    uw = csv.writer(uf, delimiter="\t")
    uw.writerow(["build", "method", "gene", "transcript", "hgvsc_raw",
                 "coord_or_annovar", "cls_raw", "cls"])

    preview_path = os.path.join(args.outdir, "qc_parsed_preview.tsv")
    pf = open(preview_path, "w", newline="", encoding="utf-8")
    pw = csv.writer(pf, delimiter="\t")
    pw.writerow(["build", "method", "gene", "hgvsc_raw", "resolved_coord",
                 "cls_raw", "cls", "truth_assertion"])

    summary_rows = []

    for build in ("hg19", "hg38"):
        out("=" * 78)
        out(f"Build: {build}")
        out("=" * 78)

        truth = P.build_truth(reg["builds"][build]["vcf"], orig)
        st = truth["stats"]
        out("[VCF truth]")
        out(f"  VCF rows={st['n_vcf']}  linked to original={st['n_in_original']} "
            f"({pct(st['n_in_original'], st['n_vcf'])})")
        out(f"  unique coordinates with truth={len(truth['coord2assert'])}  "
            f"conflicting excluded={st['n_conflict']}")
        tdist = Counter(truth["coord2assert"].values())
        out("  truth distribution (5 cat): " +
            ", ".join(f"{k}={tdist.get(k,0)}" for k in P.CANON_ORDER))
        n_truth = len(truth["coord2assert"])
        out("")

        for name, loader in build_methods(reg, build):
            recs = loader()
            n = len(recs)
            cls_c = Counter(r["cls"] for r in recs)
            unknown = cls_c.get(None, 0)
            no_cls = sum(1 for r in recs if r["cls"] is None)

            resolved = 0
            with_truth = 0
            covered = set()
            sample_left = args.samples
            preview_left = 3
            for r in recs:
                a, coord = P.resolve_truth(r, truth)
                if coord is not None:
                    resolved += 1
                    if a is not None:
                        with_truth += 1
                        covered.add(coord)
                        if preview_left > 0:
                            pw.writerow([build, name, r["gene"], r["hgvsc_raw"],
                                         coord, r["cls_raw"], r["cls"], a])
                            preview_left -= 1
                else:
                    if sample_left > 0:
                        coa = (":".join(map(str, r["annovar"]))
                               if r["source"] == "intervar" else r.get("coord"))
                        uw.writerow([build, name, r["gene"], r["transcript"],
                                     r["hgvsc_raw"], coa, r["cls_raw"], r["cls"]])
                        sample_left -= 1

            out(f"--- {name} ---")
            out(f"  rows: {n}")
            out("  classification distribution (5 cat): " +
                ", ".join(f"{k}={cls_c.get(k,0)}" for k in P.CANON_ORDER) +
                (f", [unrecognised={unknown}]" if unknown else ""))
            if unknown:
                for s, c in Counter(r["cls_raw"] for r in recs
                                    if r["cls"] is None).most_common(5):
                    out(f"      [!] unrecognised label {s!r}: {c}")
            if no_cls and not unknown:
                out(f"  missing classification (blank): {no_cls}")
            out(f"  resolved to VCF coordinates: {resolved}/{n} ({pct(resolved,n)})")
            out(f"  with truth (evaluated): {with_truth}  "
                f"/ truth coverage {len(covered)}/{n_truth} ({pct(len(covered),n_truth)})")
            out("")
            summary_rows.append((build, name, n, with_truth,
                                 pct(len(covered), n_truth)))

    uf.close()
    pf.close()

    out("=" * 78)
    out("Number of evaluated variants (= variants that can be matched to the truth)")
    out("=" * 78)
    out(f"{'build':5s} {'method':28s} {'rows':>6s} {'eval':>6s} {'truth cov.':>10s}")
    for b, m, n, ev, cov in summary_rows:
        out(f"{b:5s} {m:28s} {n:6d} {ev:6d} {cov:>10s}")
    out("")
    out("Checks:")
    out(" - every classification label is mapped to one of the 5 categories (unrecognised = 0)")
    out(" - the resolution rate to VCF coordinates is close to 100% for each method")
    out(f" - unresolved examples: {unmatched_path}")
    out(f" - normalisation preview: {preview_path}")
    out("")
    out("If there are no problems, this VCF truth is used as the reference for")
    out("the performance (P/LP, B/LB, VUS) and per-code analyses.")

    report_path = os.path.join(args.outdir, "qc_report.txt")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n[written] {report_path}")


if __name__ == "__main__":
    main()
