#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_all.py  -- run the eRepo benchmark analysis in one command.

Steps:
  1. verify_parsing.py   : sanity check of parsing and VCF coordinate matching (QC; optional)
  2. analyze_benchmark.py: category-level performance / assignments per ACMG code / F1 per code + figures

Input files (all directly in the --input directory):
  erepo-original.csv
  erepo_hg19.vcf / erepo_hg38.vcf                         (TransVar conversion, NAME=GENE:HGVSc)
  erepo_{hg19,hg38}_{hg19,hg38}_multianno_txt.intervar    (ANNOVAR+InterVar)
  erepo_{hg19,hg38}_bias-2015.tsv                         (BIAS-2015)
  erepo_hg19_bias-2015_ignore_PP5_BP6.tsv                 (BIAS-2015 with PP5/BP6 disabled, hg19)
  erepo_hg38_bias-2015_PP5_BP6.tsv                        (BIAS-2015 with PP5/BP6 disabled, hg38; named differently from hg19)
  erepo_{hg19,hg38}_huhvar_default_results.tsv            (HUHVar default)
  erepo_{hg19,hg38}_huhvar_revel_results.tsv              (HUHVar REVEL)
  erepo_{hg19,hg38}_huhvar_spliceai_results.tsv           (HUHVar spliceAI)

Usage:
  python3 run_all.py --input ./uploads --outdir ./outputs
  python3 run_all.py --input inputs --outdir outputs
  python3 run_all.py --no-figs            # do not draw figures
  python3 run_all.py --skip-qc            # skip the QC step
"""

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, extra):
    cmd = [sys.executable, os.path.join(HERE, script)] + extra
    print(f"\n$ {' '.join(cmd)}\n" + "-" * 70)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit(f"[ERROR] {script} failed (code={r.returncode})")


def main():
    ap = argparse.ArgumentParser(description="Run the eRepo benchmark analysis")
    ap.add_argument("--input", default="inputs",
                    help="directory with the input files")
    ap.add_argument("--outdir", default="outputs",
                    help="output directory")
    ap.add_argument("--no-figs", action="store_true", help="do not draw figures")
    ap.add_argument("--skip-qc", action="store_true", help="skip QC")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    base = ["--input", args.input, "--outdir", args.outdir]
    if not args.skip_qc:
        run("verify_parsing.py", base)
    ana = base + (["--no-figs"] if args.no_figs else [])
    run("analyze_benchmark.py", ana)

    print("\n" + "=" * 70)
    print("Done. Main outputs:")
    for fn in ["qc_report.txt",
               "analysis_report.txt",
               "analysis_category_metrics.tsv",
               "analysis_acmg_code_counts.tsv",
               "analysis_code_f1.tsv",
               "analysis_confusion.txt",
               "category_f1_hg19.png", "category_f1_hg38.png",
               "acmg_code_counts_hg19.png", "acmg_code_counts_hg38.png",
               "code_f1_hg19.png", "code_f1_hg38.png"]:
        p = os.path.join(args.outdir, fn)
        mark = "OK " if os.path.exists(p) else "-- "
        print(f"  [{mark}] {fn}")


if __name__ == "__main__":
    main()
