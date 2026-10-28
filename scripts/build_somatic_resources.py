#!/usr/bin/env python3
"""Build the bundled tables used by the somatic oncogenicity classifier.

1. ``resources/shared/somatic_gene_roles.tsv`` — cancer gene roles from the
   Network of Cancer Genes (NCG 7, Dressler et al. 2022), canonical drivers only
   (genes annotated by the Cancer Gene Census, Vogelstein et al. 2013 or Saito
   et al. 2020 lists inside NCG). The role is NCG's majority vote over those
   three lists (``NCG_oncogene`` / ``NCG_tsg``); canonical genes without a
   majority role are ``unresolved``. ``somatic_gene_role_overrides.tsv`` is
   applied last.

     python scripts/build_somatic_resources.py --ncg NCG_cancerdrivers_annotation_supporting_evidence.tsv

2. ``resources/shared/cancer_hotspots_v2.tsv`` — single-residue and in-frame
   indel hotspots of cancerhotspots.org v2 (Chang et al. 2018; 24,592 tumours),
   with the per-position and per-amino-acid-change sample counts used by the
   OS3 / OM3 / OP3 criteria. Data licence: ODbL 1.0 (see NOTICE).

     python scripts/build_somatic_resources.py --hotspots hotspots_v2.xls

``hotspots_v2.xls`` is in the cancerhotspots repository
(github.com/cBioPortal/cancerhotspots, webapp/src/main/resources/static/files/).
Reading the .xls needs ``xlrd``.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

_SHARED = Path(__file__).resolve().parents[1] / "resources" / "shared"


def _read_overrides(path: Path) -> dict[str, tuple[str, str]]:
    out: dict[str, tuple[str, str]] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            out[r["gene"].strip()] = (r["role"].strip(), r.get("rationale", "").strip())
    return out


def build_gene_roles(ncg: Path, overrides: Path, out: Path, version: str) -> int:
    by_gene: dict[str, list[dict]] = defaultdict(list)
    with ncg.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by_gene[r["symbol"].strip()].append(r)
    ov = _read_overrides(overrides)
    rows = []
    for gene in sorted(by_gene):
        rs = [r for r in by_gene[gene]
              if r["cgc_annotation"] or r["vogelstein_annotation"] or r["saito_annotation"]]
        if not rs:
            continue  # candidate driver: not used
        r = rs[0]

        def _clean(v: str) -> str:
            return re.sub(r"\s+", " ", (v or "").strip().strip('"'))

        cgc, vog, sai = _clean(r["cgc_annotation"]), _clean(r["vogelstein_annotation"]), _clean(r["saito_annotation"])
        onc, tsg = r["NCG_oncogene"] == "1", r["NCG_tsg"] == "1"
        role = "TSG" if tsg and not onc else "oncogene" if onc and not tsg else "unresolved"
        source = "NCG"
        note = ""
        cgc_roles = {x.strip().lower() for x in cgc.split(",")}
        if role == "TSG" and "oncogene" in cgc_roles:
            note = "CGC also lists an oncogene role"
        elif role == "oncogene" and "tsg" in cgc_roles:
            note = "CGC also lists a TSG role"
        if gene in ov:
            role, rationale = ov.pop(gene)
            source = "override"
            note = rationale
        rows.append([gene, role, source, cgc, vog, sai, note])
    for gene, (role, rationale) in sorted(ov.items()):  # override genes absent from NCG
        rows.append([gene, role, "override", "", "", "", rationale])
    with out.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# Cancer gene roles for somatic oncogenicity (OVS1 / OM2). Source: Network of "
                 f"Cancer Genes {version} (Dressler et al. 2022, canonical drivers; role = NCG majority "
                 "vote of CGC / Vogelstein 2013 / Saito 2020), corrected by "
                 "somatic_gene_role_overrides.tsv.\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "role", "source", "cgc", "vogelstein", "saito", "note"])
        w.writerows(rows)
    return len(rows)


def _count(field: str) -> tuple[str, int]:
    """'R:204' -> ('R', 204)."""
    k, _, n = str(field).rpartition(":")
    return k, int(float(n))


def build_hotspots(xls: Path, out: Path) -> int:
    import xlrd

    wb = xlrd.open_workbook(str(xls))
    rows = []
    for sheet, kind in (("SNV-hotspots", "snv"), ("INDEL-hotspots", "indel")):
        ws = wb.sheet_by_name(sheet)
        hdr = ws.row_values(0)
        for i in range(1, ws.nrows):
            r = dict(zip(hdr, ws.row_values(i)))
            pos = str(r["Amino_Acid_Position"]).strip()
            if pos.endswith(".0"):
                pos = pos[:-2]
            if "splice" in pos:
                continue
            gene = str(r["Hugo_Symbol"]).strip()
            n_pos = int(float(r["Mutation_Count"]))
            var, n_var = _count(r["Variant_Amino_Acid"])
            if kind == "snv":
                ref, _ = _count(str(r["Reference_Amino_Acid"]).split("|")[0])
            else:
                ref = ""
            rows.append([gene, kind, pos, ref, var, n_var, n_pos])
    rows.sort(key=lambda x: (x[0], x[1], int(re.match(r"\d+", x[2]).group()), -x[5]))
    with out.open("w", encoding="utf-8", newline="") as fh:
        fh.write("# cancerhotspots.org v2 (Chang et al. Cancer Discov 2018; 24,592 tumours). "
                 "Data licence: Open Database License (ODbL) 1.0. variant_count = samples with this "
                 "amino-acid change; position_count = samples with a mutation at the residue "
                 "(snv) or in the indel hotspot region (indel).\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "type", "position", "ref_aa", "variant", "variant_count", "position_count"])
        w.writerows(rows)
    return len(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ncg", type=Path, help="NCG_cancerdrivers_annotation_supporting_evidence.tsv")
    ap.add_argument("--ncg-version", default="7.1")
    ap.add_argument("--overrides", type=Path, default=_SHARED / "somatic_gene_role_overrides.tsv")
    ap.add_argument("--hotspots", type=Path, help="cancerhotspots hotspots_v2.xls")
    ap.add_argument("--out-dir", type=Path, default=_SHARED)
    a = ap.parse_args()
    if not a.ncg and not a.hotspots:
        ap.error("give --ncg and/or --hotspots")
    if a.ncg:
        n = build_gene_roles(a.ncg, a.overrides, a.out_dir / "somatic_gene_roles.tsv", a.ncg_version)
        print(f"{n} genes -> {a.out_dir / 'somatic_gene_roles.tsv'}")
    if a.hotspots:
        n = build_hotspots(a.hotspots, a.out_dir / "cancer_hotspots_v2.tsv")
        print(f"{n} hotspot variants -> {a.out_dir / 'cancer_hotspots_v2.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
