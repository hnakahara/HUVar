#!/usr/bin/env python3
"""Build the VCEP criteria provenance table.

For every gene the app applies VCEP-specific ACMG criteria to, record the
governing ClinGen specification (GN id), its version and status, the VCEP name,
and which criteria are implemented for that gene. Sources:
  * resources/shared/disease_prevalence.tsv — per-gene criterion columns +
    cspec_url (the GN id of the governing spec the build resolved for the gene).
  * resources/shared/pm1_hotspots.tsv       — PM1 hotspot rows.
  * resources/shared/pm4_regions.tsv        — PM4 region / gate rows.
  * src/acmg_classifier/pvs1/vcep_pvs1.py   — gene-specific PVS1 decision trees.
  * resources/clingen/cspec_json/GN*.json   — spec version / status / label.

Writes a TSV and a Markdown table to docs/.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DP = ROOT / "resources" / "shared" / "disease_prevalence.tsv"
MULTISPEC = ROOT / "resources" / "shared" / "disease_prevalence_multispec.tsv"
PM1 = ROOT / "resources" / "shared" / "pm1_hotspots.tsv"
PM4 = ROOT / "resources" / "shared" / "pm4_regions.tsv"
PVS1 = ROOT / "src" / "acmg_classifier" / "pvs1" / "vcep_pvs1.py"
JSON_DIR = ROOT / "resources" / "clingen" / "cspec_json"
OUT_TSV = ROOT / "docs" / "vcep_criteria_provenance.tsv"
OUT_MD = ROOT / "docs" / "vcep_criteria_provenance.md"

# A criterion is "implemented" for a gene when any of these columns is non-blank.
_COL_CRITERIA: dict[str, list[str]] = {
    "BA1": ["ba1_threshold"],
    "BS1": ["bs1_threshold", "bs1_strength", "bs1_exclude"],
    "BS2": ["bs2", "bs2_strength"],
    "PM2": ["pm2_threshold", "pm2_strength", "pm2_subpop", "pm2_zygosity"],
    "PM4": ["pm4", "pm4_supporting_max_aa"],
    "PM5": ["pm5_grantham", "pm5_excludes", "pm5_max", "pm5_lp"],
    "PP2": ["pp2", "pp2_requires"],
    "PP3": ["revel_pp3_supporting", "revel_pp3_moderate", "revel_pp3_strong"],
    "BP4": ["revel_bp4_supporting", "revel_bp4_moderate", "revel_bp4_strong",
            "bp4_splice_cutoff"],
    "BP7": ["bp7_phylop", "bp7_intronic", "bp7_splice_cutoff"],
    "PVS1": ["pvs1"],
    "PS1": ["ps1", "ps1_splice", "ps1_max", "ps1_paralog_group"],
    "BP1": ["bp1"],
    "BP3": ["bp3"],
}


def _gn_from_url(url: str) -> str:
    m = re.search(r"(GN\d+)", url or "")
    return m.group(1) if m else ""


def _spec_meta() -> dict[str, tuple[str, str, str]]:
    """GN id -> (version, status, vcep_label)."""
    out: dict[str, tuple[str, str, str]] = {}
    for f in JSON_DIR.glob("GN*.json"):
        try:
            d = json.load(open(f, encoding="utf-8"))
        except Exception:
            continue
        gn = f.stem
        status = d.get("currentStatus")
        status = status.get("label") if isinstance(status, dict) else (status or "")
        label = (d.get("label") or "").replace(
            "ClinGen ", "").replace(
            " Specifications to the ACMG/AMP Variant Interpretation Guidelines", "")
        # Sanitise for a Markdown/TSV cell: no newlines or pipes, collapse spaces.
        label = re.sub(r"\s+", " ", label.replace("|", "/")).strip()
        out[gn] = (d.get("version", ""), status or "", label[:70])
    return out


def _multispec_rows(meta: dict[str, tuple[str, str, str]],
                    primary_gn: dict[str, str]) -> list[dict]:
    """Alternative per-CSpec specs for the multi-disease genes (RYR1/ACTA1/VWF).

    Sourced from ``disease_prevalence_multispec.tsv`` (one row per gene×CSpec).
    The per-gene table only carries the single conservative spec the batch build
    resolved; these are the additional CSpecs the HUVar app can switch a variant
    to. Empty if the multispec table is absent.
    """
    if not MULTISPEC.exists():
        return []
    out: list[dict] = []
    for r in csv.DictReader(open(MULTISPEC, encoding="utf-8"), delimiter="\t"):
        gene = (r.get("gene_symbol") or "").strip()
        gn = (r.get("source_gn") or "").strip()
        ver, status, _ = meta.get(gn, ("", "", ""))
        out.append({
            "gene": gene,
            "disease": (r.get("disease_label") or "").strip(),
            "cspec_id": (r.get("cspec_id") or "").strip(),
            "GN_ID": gn,
            "version": ver,
            "status": status,
            "conservative": "✓" if primary_gn.get(gene) == gn else "",
        })
    return out


def main() -> None:
    meta = _spec_meta()

    pm1_genes = {r["gene_symbol"] for r in csv.DictReader(open(PM1, encoding="utf-8"), delimiter="\t")}
    pm4_genes = {r["gene_symbol"] for r in csv.DictReader(open(PM4, encoding="utf-8"), delimiter="\t")}
    pvs1_genes = set(re.findall(r'^\s*"([A-Z0-9]+)":\s*_GeneSpec\(', PVS1.read_text(encoding="utf-8"), re.M))

    rows = []
    for r in csv.DictReader(open(DP, encoding="utf-8"), delimiter="\t"):
        gene = r["gene_symbol"]
        crits = []
        for crit, cols in _COL_CRITERIA.items():
            if any((r.get(c) or "").strip() for c in cols):
                crits.append(crit)
        if gene in pm1_genes:
            crits.append("PM1")
        if gene in pm4_genes and "PM4" not in crits:
            crits.append("PM4")
        if gene in pvs1_genes and "PVS1" not in crits:
            crits.append("PVS1")
        if not crits:
            continue
        gn = _gn_from_url(r.get("cspec_url", ""))
        ver, status, label = meta.get(gn, ("", "", r.get("source_vcep", "")))
        rows.append({
            "gene": gene, "GN_ID": gn, "version": ver, "status": status,
            "VCEP": label or r.get("source_vcep", ""),
            "criteria": " ".join(sorted(set(crits))),
        })

    rows.sort(key=lambda x: x["gene"])
    multispec = _multispec_rows(meta, {r["gene"]: r["GN_ID"] for r in rows})
    fields = ["gene", "GN_ID", "version", "status", "VCEP", "criteria"]
    with open(OUT_TSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # Cross-spec provenance notes: criteria whose source spec differs from the
    # gene's primary cspec_url (curated, recorded in the build scripts / code).
    cross = [
        ("PS1 paralogue", "HRAS, NRAS, KRAS, MAP2K1, MAP2K2, SOS1, SOS2", "GN004",
         meta.get("GN004", ("", "", ""))[0], "RASopathy VCEP"),
        ("PS1 paralogue", "HBA2 (← HBA1)", "GN173", meta.get("GN173", ("", "", ""))[0],
         "Hemoglobinopathy VCEP"),
        ("PS1 paralogue (analogous-residue alignment)",
         "SCN1A, SCN2A, SCN3A, SCN8A", "GN067-070",
         meta.get("GN067", ("", "", ""))[0], "SCN epilepsy VCEPs"),
        ("PS1 paralogue (analogous-residue alignment, → KCNQ2, Moderate)",
         "KCNQ1", "GN112", meta.get("GN112", ("", "", ""))[0], "Long QT VCEP"),
    ]

    with open(OUT_MD, "w", encoding="utf-8") as fh:
        fh.write("# VCEP criteria provenance\n\n")
        fh.write("Per gene: the governing ClinGen specification (GN id + version) the app's "
                 "VCEP-specific ACMG criteria are based on, and which criteria are implemented "
                 "for that gene. Auto-generated by `scripts/build_criteria_provenance.py` from "
                 "the resource TSVs and `resources/clingen/cspec_json/`.\n\n")
        fh.write(f"Genes covered: **{len(rows)}**.\n\n")
        fh.write("| Gene | GN_ID | Version | Status | VCEP | Implemented criteria |\n")
        fh.write("|------|-------|---------|--------|------|----------------------|\n")
        for x in rows:
            fh.write(f"| {x['gene']} | {x['GN_ID']} | {x['version']} | {x['status']} "
                     f"| {x['VCEP']} | {x['criteria']} |\n")
        fh.write("\n## Cross-spec criteria (source differs from the gene's primary GN)\n\n")
        fh.write("| Criterion | Genes | GN_ID | Version | VCEP |\n")
        fh.write("|-----------|-------|-------|---------|------|\n")
        for crit, genes, gn, ver, vcep in cross:
            fh.write(f"| {crit} | {genes} | {gn} | {ver} | {vcep} |\n")

        if multispec:
            fh.write("\n## Alternative CSpecs (multi-disease genes)\n\n")
            fh.write("A few genes carry several Released, clinically distinct ClinGen "
                     "CSpecs. The per-gene table above lists only the single "
                     "**conservative** spec the batch `disease_prevalence.tsv` resolves to "
                     "(used by the CLI and batch runs, marked ✓ below); the HUVar app can "
                     "additionally evaluate these genes under any of the alternatives below, "
                     "via `resources/shared/disease_prevalence_multispec.tsv`.\n\n")
            fh.write("| Gene | Disease / mode | GN_ID | Version | Status | Conservative |\n")
            fh.write("|------|----------------|-------|---------|--------|--------------|\n")
            for x in multispec:
                fh.write(f"| {x['gene']} | {x['disease']} | {x['GN_ID']} | {x['version']} "
                         f"| {x['status']} | {x['conservative']} |\n")

    print(f"provenance rows: {len(rows)} | multispec rows: {len(multispec)} "
          f"| written → {OUT_TSV} and {OUT_MD}")


if __name__ == "__main__":
    main()
