#!/usr/bin/env python3
"""Validate the expert-panel criterion import (PS3/PS4/PP1/BS2) against the eRepo,
old rule vs revised rule, directly from the ClinVar RCV XML (no SQLite rebuild).

old: all comments of a >=3-star ClinVarSet (any submitter), whole-text negation
     patterns (HUVar as submitted / before 2026-10-09).
new: expert-panel SCV comments only, clause-level negation
     (clinvar_builder._ep_scv_texts / _mine_ep_criteria / _mine_bs2 in the
     installed HUVar).

  python ep_import_validate.py --xml $DATA/GRCh38/clinvar/ClinVarRCVRelease.xml.gz \
     --results erepo_benchmark/inputs/erepo_hg38_huhvar_revel_manual_results.tsv \
     --erepo <HUVar>/resources/GRCh38/erepo_classifications_hg38.tsv \
     --out ep_validate --workers 16

Reference: the eRepo met codes of the same variant. Note that the ClinVar record
and the eRepo snapshot can differ in version, so 100% agreement is not expected.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import xml.etree.ElementTree as ET
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

from acmg_classifier.setup.clinvar_builder import (
    _batched,
    _ep_scv_texts,
    _iter_clinvarset_chunks,
    _mine_bs2,
    _mine_ep_criteria,
    _rcv_classification,
    _star_rating,
)

CODES = ("PS3", "PS4", "PP1", "BS2")

# ---- old rule (verbatim from the submitted version) -------------------------
_OLD_BS2_POS = re.compile(r"\bBS2(?:_?(VeryStrong|Very[\s_]Strong|Strong|Moderate|Supporting))?\b", re.I)
_OLD_BS2_NEG = re.compile(
    r"BS2[\s_]*(?:is\s+|was\s+|were\s+|are\s+)?(?:not\s+(?:applicable|met|applied)|n/?a\b)|"
    r"(?:not\s+applicable|does\s+not\s+meet|did\s+not\s+meet|cannot\s+(?:be\s+)?"
    r"(?:applied|met)|unable\s+to\s+(?:apply|meet))[^.]{0,30}?BS2", re.I)


def _old_pos(code):
    return re.compile(rf"\b{code}(?:[_\s]?(VeryStrong|Very[\s_]Strong|Strong|Moderate|Supporting))?\b", re.I)


def _old_neg(code):
    return re.compile(
        rf"{code}[\s_]*(?:\w+\s+){{0,1}}(?:is\s+|was\s+|were\s+)?(?:not\s+(?:applicable|met|applied)|n/?a\b)|"
        rf"(?:not\s+(?:met|applied|applicable)|does\s+not\s+meet|did\s+not\s+meet)\s*:?[^.]{{0,40}}?\b{code}\b", re.I)


_OLD = {c: (_old_pos(c), _old_neg(c)) for c in ("PS3", "PS4", "PP1")}


def old_codes(text):
    out = set()
    for c, (pos, neg) in _OLD.items():
        if not neg.search(text) and pos.search(text):
            out.add(c)
    if not _OLD_BS2_NEG.search(text) and _OLD_BS2_POS.search(text):
        out.add("BS2")
    return out


def new_codes(text):
    out = {kv.split(":")[0] for kv in (_mine_ep_criteria(text) or "").split(";") if kv}
    if _mine_bs2(text)[0]:
        out.add("BS2")
    return out


IDS: set[str] = set()


def _init(ids):
    global IDS
    IDS = ids


def _work(chunks):
    res = []
    for ch in chunks:
        m = re.search(r'<MeasureSet[^>]*\bID="(\d+)"', ch)
        if not m or m.group(1) not in IDS:
            continue
        try:
            el = ET.fromstring(ch)
        except ET.ParseError:
            continue
        ref = el.find(".//ReferenceClinVarAssertion")
        if ref is None:
            continue
        _, rev, _ = _rcv_classification(ref)
        if _star_rating(rev) < 3:
            continue
        all_text = " ".join([c.text for c in el.findall(".//Comment") if c.text]
                            + [d.text for d in el.findall(".//Attribute[@Type='Description']") if d.text])
        res.append((m.group(1), old_codes(all_text), new_codes(" ".join(_ep_scv_texts(el)))))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--xml", type=Path, required=True)
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--erepo", type=Path, required=True)
    ap.add_argument("--out", default="ep_validate")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()

    truth = {}
    with a.erepo.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            truth[r["variant_id"]] = {c.strip().split("_")[0] for c in (r.get("met_codes") or "").split(",") if c.strip()}
    vid2var = {}
    with a.results.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = (r.get("clinvar_variation_id") or "").strip()
            if v and r["variant_id"] in truth:
                vid2var[str(int(float(v)))] = r["variant_id"]

    # a variant can appear in several ClinVarSets (one per condition): union
    old, new = {}, {}
    with gzip.open(a.xml, "rt", encoding="utf-8") as fh, Pool(a.workers, initializer=_init, initargs=(set(vid2var),)) as pool:
        for res in pool.imap_unordered(_work, _batched(_iter_clinvarset_chunks(fh), 200), chunksize=4):
            for vid, o, n in res:
                old.setdefault(vid, set()).update(o)
                new.setdefault(vid, set()).update(n)

    lines = [f"eRepo variants with a ClinVar variation ID: {len(vid2var)}; with a >=3-star ClinVar record: {len(old)}",
             "code  rule  imported  agree_eRepo  precision  eRepo_code(n, among >=3-star)  recall"]
    rows = []
    for c in CODES:
        tr = {v for v in old if c in truth[vid2var[v]]}
        for name, d in (("old", old), ("new", new)):
            imp = {v for v, s in d.items() if c in s}
            agree = imp & tr
            prec = len(agree) / len(imp) if imp else float("nan")
            rec = len(agree) / len(tr) if tr else float("nan")
            lines.append(f"{c:4}  {name:4}  {len(imp):8}  {len(agree):11}  {prec:9.3f}  {len(tr):29}  {rec:6.3f}")
            for v in sorted(imp - tr):
                rows.append(dict(rule=name, code=c, variation_id=v, variant_id=vid2var[v], kind="imported_not_in_eRepo"))
    with open(f"{a.out}_disagreements.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["rule", "code", "variation_id", "variant_id", "kind"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    txt = "\n".join(lines) + "\n"
    Path(f"{a.out}_summary.txt").write_text(txt, encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
