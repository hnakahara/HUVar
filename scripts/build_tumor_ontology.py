#!/usr/bin/env python3
"""Build the tumour-type tables used for the AMP/ASCO/CAP 2017 tiers.

* ``resources/shared/oncotree_doid.tsv`` — OncoTree code -> Disease Ontology IDs.
  Order of mapping: curated overrides (``oncotree_doid_overrides.tsv``), the node's
  NCIt cross-reference, its UMLS CUI, its exact name (DO name or EXACT synonym),
  then the nearest mapped ancestor node (method ``parent:<CODE>``).
* ``resources/shared/doid_hierarchy.tsv`` — DOID, name, is_a parents (non-obsolete
  terms), used to decide whether a CIViC evidence disease is the patient's tumour
  type, a more general one (e.g. "cancer" for tumour-agnostic approvals) or a
  subtype.

Inputs: an OncoTree release JSON (github.com/cBioPortal/oncotree, trees/, CC BY 4.0)
and the Disease Ontology ``doid.obo`` (github.com/DiseaseOntology/HumanDiseaseOntology,
src/ontology/releases/, CC0).

  python scripts/build_tumor_ontology.py --oncotree oncotree_2025_10_03.json --doid doid.obo
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

_SHARED = Path(__file__).resolve().parents[1] / "resources" / "shared"


def read_doid(path: Path) -> tuple[dict, str]:
    text = path.read_text(encoding="utf-8")
    version = (re.search(r"^data-version: (.*)$", text, re.M) or [None, ""])[1]
    terms = {}
    for blk in text.split("\n[Term]")[1:]:
        m = re.search(r"^id: (DOID:\d+)", blk, re.M)
        if not m or "is_obsolete: true" in blk:
            continue
        terms[m.group(1)] = {
            "name": re.search(r"^name: (.*)$", blk, re.M).group(1).strip(),
            "is_a": re.findall(r"^is_a: (DOID:\d+)", blk, re.M),
            "xref": re.findall(r"^xref: (\S+)", blk, re.M),
            "syn": re.findall(r'^synonym: "([^"]+)" EXACT', blk, re.M),
        }
    return terms, version


def build(oncotree: Path, doid: Path, overrides: Path, out_dir: Path) -> dict:
    terms, do_version = read_doid(doid)
    idx: dict[str, set] = defaultdict(set)
    for k, t in terms.items():
        for x in t["xref"]:
            if x.startswith("NCI:"):
                idx["NCI:" + x[4:]].add(k)
            elif x.startswith("UMLS_CUI:"):
                idx["UMLS:" + x[9:]].add(k)
        for nm in [t["name"]] + t["syn"]:
            idx["NAME:" + nm.lower()].add(k)
    ov = {}
    if overrides.exists():
        with overrides.open(encoding="utf-8") as fh:
            for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
                ds = [d.strip() for d in r["doids"].split(",") if d.strip()]
                missing = [d for d in ds if d not in terms]
                if missing:
                    raise SystemExit(f"override {r['oncotree_code']}: unknown DOIDs {missing}")
                ov[r["oncotree_code"].strip()] = set(ds)
    tree = json.loads(oncotree.read_text(encoding="utf-8"))
    nodes: dict[str, tuple[dict, str | None]] = {}

    def walk(n, parent=None):
        nodes[n["code"]] = (n, parent)
        for c in (n.get("children") or {}).values():
            walk(c, n["code"])

    walk(tree["TISSUE"] if "TISSUE" in tree else next(iter(tree.values())))

    def direct(code):
        if code in ov:
            return ov[code], "curated"
        n = nodes[code][0]
        ext = n.get("externalReferences") or {}
        for src in ("NCI", "UMLS"):
            s = set().union(*(idx.get(f"{src}:{r}", set()) for r in ext.get(src, [])))
            if s:
                return s, src.lower()
        s = idx.get("NAME:" + n["name"].lower(), set())
        return (set(s), "name") if s else (set(), "")

    rows = []
    for code in nodes:
        s, method = direct(code)
        c = code
        while not s and nodes[c][1]:
            c = nodes[c][1]
            s, _ = direct(c)
            method = f"parent:{c}" if s else ""
        n, parent = nodes[code]
        rows.append([code, n["name"], n.get("mainType") or "", parent or "",
                     ",".join(sorted(s)), "; ".join(terms[d]["name"] for d in sorted(s)), method])
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "oncotree_doid.tsv").open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# OncoTree ({oncotree.name}; CC BY 4.0) -> Disease Ontology ({do_version}; CC0)\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["oncotree_code", "name", "main_type", "parent", "doids", "doid_names", "method"])
        w.writerows(sorted(rows))
    with (out_dir / "doid_hierarchy.tsv").open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# Disease Ontology {do_version} (CC0): non-obsolete terms and is_a parents\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["doid", "name", "parents", "synonyms"])
        for k in sorted(terms):
            t = terms[k]
            w.writerow([k, t["name"], ",".join(t["is_a"]), "|".join(t["syn"])])
    return {"oncotree": len(rows), "doid": len(terms), "do_version": do_version,
            "unmapped": sum(1 for r in rows if not r[4])}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--oncotree", type=Path, required=True)
    ap.add_argument("--doid", type=Path, required=True)
    ap.add_argument("--overrides", type=Path, default=_SHARED / "oncotree_doid_overrides.tsv")
    ap.add_argument("--out-dir", type=Path, default=_SHARED)
    a = ap.parse_args()
    info = build(a.oncotree, a.doid, a.overrides, a.out_dir)
    print(f"{info['oncotree']} OncoTree codes ({info['unmapped']} unmapped), "
          f"{info['doid']} DO terms ({info['do_version']}) -> {a.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
