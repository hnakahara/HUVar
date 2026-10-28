#!/usr/bin/env python3
"""Build the MaveDB functional-evidence table used for PS3/BS3.

Reads a MaveDB bulk dump (``mavedb-dump.<timestamp>.tar.gz`` from Zenodo, or the
extracted directory). The dump holds only published, CC0-licensed data.

Kept: score sets that have at least one *clinical* score calibration — not
research-use-only, with ACMG/AMP classifications (PS3/BS3 + evidence strength)
on its functional ranges — and are not superseded by another score set in the
dump. ExCALIBR and other research-use-only calibrations are dropped.

For every mapped variant (current MaveDB mapping) the score is placed in each
clinical calibration's functional ranges (numeric ranges, or class labels read
from the score set's own class column), giving one row per (variant,
calibration). The runtime (``criteria/functional_assays.py``) picks the
calibration (primary / most conservative) and re-derives the strength from the
range's OddsPath with the VCEP or Brnich et al. 2019 thresholds.

Outputs (``--out-dir``, default resources/shared):
  mavedb_functional.tsv.gz   gene, urn, level, key, score, cal_urn, primary,
                             criterion, mavedb_strength, oddspath, label
       level = genomic (key chrN:pos:ref:alt, GRCh38) or protein (key NP_xxx:p.Xxx123Yyy)
  mavedb_scoresets.tsv       one row per kept score set (metadata, calibrations)

Usage:
  python scripts/build_mavedb_functional.py --dump mavedb-dump.2026062418131.tar.gz
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import sys
import tarfile
from pathlib import Path

_NC_CHROM = {f"NC_0000{i:02d}": f"chr{i}" for i in range(1, 23)}
_NC_CHROM.update({"NC_000023": "chrX", "NC_000024": "chrY"})
_HGVS_G_SNV = re.compile(r"^(NC_\d{6})\.\d+:g\.(\d+)([ACGT])>([ACGT])$")
_HGVS_P = re.compile(r"^(N[PM]_\d+)(?:\.\d+)?:(p\.[A-Za-z]{3}\d+(?:[A-Za-z]{3}|=))$")


def _iter_dump(dump: Path):
    """Yield (relative path, bytes-reader) for the files of the dump (tar.gz or dir)."""
    if dump.is_dir():
        for p in sorted(dump.rglob("*")):
            if p.is_file() and not p.name.startswith("._"):
                yield str(p.relative_to(dump)), p.open("rb")
        return
    with tarfile.open(dump, "r|gz") as tf:
        for m in tf:
            name = m.name.split("/", 1)[-1]
            if not m.isfile() or Path(name).name.startswith("._"):
                continue
            yield name, tf.extractfile(m)


def _clinical(cal: dict) -> bool:
    return (not cal.get("researchUseOnly")
            and any(f.get("acmgClassification") for f in cal.get("functionalClassifications") or []))


def _in_range(v: float, f: dict) -> bool:
    lo, hi = (f.get("range") or [None, None])
    if lo is not None and (v < lo or (v == lo and not f.get("inclusiveLowerBound"))):
        return False
    if hi is not None and (v > hi or (v == hi and not f.get("inclusiveUpperBound"))):
        return False
    return True


def _class_column(cal: dict, rows: list[dict]) -> str | None:
    """Score column whose values are the calibration's class labels."""
    labels = {f.get("class") for f in cal["functionalClassifications"] if f.get("class")}
    if not labels or not rows:
        return None
    for col in rows[0]:
        if col.startswith("scores.") and labels & {r[col] for r in rows[:5000]}:
            return col
    return None


def _key(mapped: dict) -> tuple[str, str] | None:
    pm = mapped.get("postMapped") or {}
    if pm.get("type") != "Allele":
        return None
    for e in pm.get("expressions") or []:
        v = e.get("value") or ""
        m = _HGVS_G_SNV.match(v)
        if m and m.group(1) in _NC_CHROM:
            return "genomic", f"{_NC_CHROM[m.group(1)]}:{m.group(2)}:{m.group(3)}:{m.group(4)}"
        m = _HGVS_P.match(v)
        if m:
            return "protein", f"{m.group(1)}:{m.group(2)}"
    return None


def build(dump: Path, out_dir: Path) -> dict:
    meta = None
    wanted: dict[str, dict] = {}
    mapped: dict[str, list] = {}
    scores: dict[str, list] = {}
    for name, fh in _iter_dump(dump):
        if name == "main.json":
            meta = json.load(fh)
            sets = [s for es in meta["experimentSets"] for e in es["experiments"]
                    for s in e["scoreSets"]]
            urns = {s["urn"] for s in sets}
            for s in sets:
                if not any(_clinical(c) for c in s.get("scoreCalibrations") or []):
                    continue
                sup = (s.get("supersedingScoreSet") or {})
                sup_urn = sup.get("urn") if isinstance(sup, dict) else sup
                if sup_urn and sup_urn in urns:
                    continue
                if len(s.get("targetGenes") or []) != 1:
                    continue
                wanted[s["urn"].replace(":", "-")] = s
            continue
        if meta is None:
            continue
        base = Path(name).name
        prefix = base.split(".")[0]
        if prefix not in wanted:
            continue
        if base.endswith(".mapped-variants.json"):
            mapped[prefix] = json.load(fh)
        elif base.endswith(".scores.csv"):
            scores[prefix] = list(csv.DictReader(io.StringIO(fh.read().decode("utf-8"))))
    if meta is None:
        raise SystemExit("main.json not found in the dump")

    out_dir.mkdir(parents=True, exist_ok=True)
    rows_out, sets_out = [], []
    for prefix, s in sorted(wanted.items()):
        t = s["targetGenes"][0]
        gene = t.get("mappedHgncName") or t.get("name") or ""
        cals = [c for c in s["scoreCalibrations"] if _clinical(c)]
        score_rows = scores.get(prefix) or []
        by_acc = {r["accession"]: r for r in score_rows}
        class_cols = {c["urn"]: _class_column(c, score_rows) for c in cals
                      if any(f.get("range") is None for f in c["functionalClassifications"])}
        usable = [c for c in cals if c["urn"] not in class_cols or class_cols[c["urn"]]]
        n = 0
        levels = set()
        for m in mapped.get(prefix) or []:
            if not m.get("current"):
                continue
            k = _key(m)
            r = by_acc.get(m.get("variantUrn"))
            if k is None or r is None:
                continue
            try:
                v = float(r["scores.score"])
            except (TypeError, ValueError):
                v = None
            levels.add(k[0])
            for c in usable:
                col = class_cols.get(c["urn"])
                hit = None
                for f in c["functionalClassifications"]:
                    if col:
                        if f.get("class") and r.get(col) == f["class"]:
                            hit = f
                            break
                    elif v is not None and f.get("range") and _in_range(v, f):
                        hit = f
                        break
                if col is None and v is None:
                    continue
                a = (hit or {}).get("acmgClassification") or {}
                rows_out.append([
                    gene, s["urn"], k[0], k[1], "" if v is None else f"{v:.6g}",
                    c["urn"], "1" if c.get("primary") else "0",
                    a.get("criterion", ""), a.get("evidenceStrength", ""),
                    "" if (hit or {}).get("oddspathsRatio") is None else f"{hit['oddspathsRatio']:.6g}",
                    (hit or {}).get("label", "") or "",
                ])
            n += 1
        pub = (s.get("primaryPublicationIdentifiers") or [{}])[0]
        authors = pub.get("authors") or []
        first = next((a["name"] for a in authors if a.get("primary")), authors[0]["name"] if authors else "")
        sets_out.append([
            gene, s["urn"], s.get("title", ""), pub.get("identifier", ""), first.split(",")[0],
            str(pub.get("publicationYear") or ""), "/".join(sorted(levels)), str(n),
            " | ".join(f"{c['urn']}{'*' if c.get('primary') else ''} {c['title']}" for c in usable),
        ])

    with gzip.open(out_dir / "mavedb_functional.tsv.gz", "wt", encoding="utf-8", newline="") as fh:
        fh.write(f"# MaveDB functional evidence (CC0); dump asOf {meta.get('asOf')}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "urn", "level", "key", "score", "cal_urn", "primary",
                    "criterion", "mavedb_strength", "oddspath", "label"])
        w.writerows(rows_out)
    with (out_dir / "mavedb_scoresets.tsv").open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# MaveDB score sets with clinical calibrations; dump asOf {meta.get('asOf')}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "urn", "title", "publication", "first_author", "year", "level",
                    "n_mapped", "calibrations"])
        w.writerows(sets_out)
    return {"asOf": meta.get("asOf"), "score_sets": len(sets_out), "rows": len(rows_out)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dump", type=Path, required=True,
                    help="MaveDB dump (.tar.gz) or its extracted directory")
    ap.add_argument("--out-dir", type=Path, default=Path("resources/shared"))
    a = ap.parse_args()
    info = build(a.dump, a.out_dir)
    print(f"MaveDB dump {info['asOf']}: {info['score_sets']} score sets, {info['rows']} rows "
          f"-> {a.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
