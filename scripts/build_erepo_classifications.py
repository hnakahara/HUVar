"""Build ``resources/<assembly>/erepo_classifications_<hg>.tsv`` — the ClinGen
Evidence Repository (eRepo) classifications keyed by genomic coordinate.

Used (i) to restrict PS1/PM5 comparators to variants classified by the VCEP for
genes whose specification requires it, and (ii) by HUVar_app to show the existing
VCEP classification of a queried variant.

Inputs: the eRepo export converted to CSV (``Gene``, ``RefSeq (NM_)``, ``HGVSc``,
``Assertion``, ``Applied Evidence Codes (Met)``) and the TransVar-converted VCF
whose INFO/ID carries ``GENE:HGVSc`` (the same files used for the benchmark).
Variants with conflicting assertions across eRepo records are kept with all
assertions joined by ``|`` and flagged ``conflict=1``.

Usage::

    python scripts/build_erepo_classifications.py \
        --erepo erepo-original.csv --vcf erepo_hg38.vcf \
        --out resources/GRCh38/erepo_classifications_hg38.tsv --snapshot 2026-05-28
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

_ASSERT = {
    "pathogenic": "Pathogenic", "likely pathogenic": "Likely Pathogenic",
    "uncertain significance": "Uncertain Significance",
    "likely benign": "Likely Benign", "benign": "Benign",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--erepo", type=Path, required=True)
    ap.add_argument("--vcf", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--snapshot", default="")
    a = ap.parse_args()

    recs: dict[str, list[dict]] = defaultdict(list)
    with a.erepo.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            key = f"{r['Gene'].strip()}:{r['HGVSc'].strip()}"
            recs[key].append(r)

    rows = []
    with a.vcf.open(encoding="utf-8") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            for r in recs.get(vid, []):
                rows.append((f"{chrom}:{pos}:{ref}:{alt}", r))

    by_var: dict[str, list[dict]] = defaultdict(list)
    for v, r in rows:
        by_var[v].append(r)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", encoding="utf-8") as fh:
        fh.write(f"# ClinGen eRepo classifications (snapshot {a.snapshot})\n")
        fh.write("variant_id\tgene\ttranscript\thgvs_c\tassertion\tmet_codes\tconflict\n")
        for v in sorted(by_var):
            rs = by_var[v]
            assertions = sorted({_ASSERT.get(r["Assertion"].strip().lower(), r["Assertion"].strip())
                                 for r in rs})
            codes = sorted({c.strip() for r in rs
                            for c in (r.get("Applied Evidence Codes (Met)") or "").split(",") if c.strip()})
            r0 = rs[0]
            fh.write("\t".join([v, r0["Gene"], r0.get("RefSeq (NM_)", ""), r0["HGVSc"],
                                "|".join(assertions), ",".join(codes),
                                "1" if len(assertions) > 1 else "0"]) + "\n")
    print(f"wrote {len(by_var)} variants -> {a.out}")


if __name__ == "__main__":
    main()
