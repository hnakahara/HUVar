#!/usr/bin/env python3
"""How often does ClinVar carry structured family / segregation data? (PP1)

Reviewer 3 (comments 6/17): PP1 can no longer be text-mined from ClinVar. This
script checks whether the *structured* family data that some submitters provide
(``ObservedIn/Sample/FamilyData``: NumFamilies, NumFamiliesWithVariant,
NumFamiliesWithSegregationObserved, SegregationObserved, PedigreeID, and any
``ObservedData`` attribute whose type mentions families or segregation) could be
used instead.

It streams the same ClinVar RCV XML that setup_data.py uses and reports, for all
variants with GRCh38/GRCh37 coordinates and for the eRepo variants:

* how many variants / SCVs have any FamilyData, and how many state that
  segregation was observed (or give a number of families with segregation);
* the same for P/LP SCVs only, and excluding ClinGen expert-panel SCVs;
* for comparison, the free-text "segregat" mention used by the old PP1 rule;
* among eRepo variants: how the structured data relate to the PP1 the VCEP
  applied (eRepo met codes).

Run on the analysis server (needs the HUVar package for the XML helpers):

    python clinvar_familydata_survey.py \
        --xml $DATA/GRCh38/clinvar/ClinVarRCVRelease.xml.gz --assembly GRCh38 \
        --erepo <HUVar>/resources/GRCh38/erepo_classifications_hg38.tsv \
        --out familydata_GRCh38 --workers 16

Outputs ``<out>_summary.txt`` and ``<out>_variants.tsv`` (variants with any
FamilyData, any segregation text, or present in the eRepo table).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

from acmg_classifier.setup.clinvar_builder import (  # same parsing as the SQLite build
    _SEGREGATION_NEG,
    _SEGREGATION_POS,
    _batched,
    _iter_clinvarset_chunks,
    _rcv_classification,
    _scv_is_expert_panel,
    _scv_significance,
    _star_rating,
)

_FAM_ATTR_RE = re.compile(r"famil|segregat|meios|pedigree", re.IGNORECASE)
_INT_ATTRS = ("NumFamilies", "NumFamiliesWithVariant", "NumFamiliesWithSegregationObserved")


def _to_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _parse(chunk: str, assembly: str):
    try:
        elem = ET.fromstring(chunk)
    except ET.ParseError:
        return None
    ref = elem.find(".//ReferenceClinVarAssertion")
    if ref is None:
        return None
    loc = ref.find(f".//MeasureSet/Measure/SequenceLocation[@Assembly='{assembly}']")
    if loc is None or not loc.get("positionVCF"):
        return None
    key = (f"chr{loc.get('Chr')}:{loc.get('positionVCF')}:"
           f"{loc.get('referenceAlleleVCF')}:{loc.get('alternateAlleleVCF')}")
    _, rev, _ = _rcv_classification(ref)
    scvs = []
    for scv in elem.findall(".//ClinVarAssertion"):
        acc = scv.find(".//ClinVarAccession")
        acc_id = acc.get("Acc", "") if acc is not None else ""
        sig = _scv_significance(scv)
        rec = {
            "acc": acc_id,
            "plp": sig.startswith("pathogenic") or sig.startswith("likely pathogenic"),
            "ep": _scv_is_expert_panel(scv),
            "fam": False, "seg_yes": False, "seg_no": False,
            "n_fam": 0, "n_fam_var": 0, "n_fam_seg": 0, "attr_types": [],
            "text_seg": False,
        }
        for fd in scv.findall(".//FamilyData"):
            rec["fam"] = True
            so = (fd.get("SegregationObserved") or "").strip().lower()
            rec["seg_yes"] |= so == "yes"
            rec["seg_no"] |= so == "no"
            n = [_to_int(fd.get(a)) for a in _INT_ATTRS]
            rec["n_fam"] += n[0] or 0
            rec["n_fam_var"] += n[1] or 0
            rec["n_fam_seg"] += n[2] or 0
        for at in scv.findall(".//ObservedData/Attribute"):
            t = at.get("Type", "")
            if _FAM_ATTR_RE.search(t):
                rec["attr_types"].append(t)
                rec["fam"] = True
                v = _to_int(at.get("integerValue") or at.text)
                if "segregation" in t.lower() and v:
                    rec["n_fam_seg"] += v
        texts = [c.text for c in scv.findall(".//Comment") if c.text]
        texts += [d.text for d in scv.findall(".//ObservedData/Attribute[@Type='Description']") if d.text]
        combined = " ".join(texts)
        rec["text_seg"] = bool(combined and _SEGREGATION_POS.search(combined)
                               and not _SEGREGATION_NEG.search(combined))
        scvs.append(rec)
    return key, _star_rating(rev), scvs


def _work(args):
    chunks, assembly = args
    out = []
    for c in chunks:
        r = _parse(c, assembly)
        if r:
            out.append(r)
    return out


def load_erepo(path: Path | None) -> dict[str, dict]:
    if not path:
        return {}
    out = {}
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            codes = {c.strip().split("_")[0] for c in (r.get("met_codes") or "").split(",") if c.strip()}
            out[r["variant_id"]] = {"gene": r.get("gene", ""), "pp1": "PP1" in codes,
                                    "assertion": r.get("assertion", "")}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--xml", type=Path, required=True)
    ap.add_argument("--assembly", default="GRCh38", choices=["GRCh38", "GRCh37"])
    ap.add_argument("--erepo", type=Path, default=None)
    ap.add_argument("--out", default="familydata")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()

    erepo = load_erepo(a.erepo)
    # per variant: max stars, SCVs keyed by accession (a SCV is listed once per RCV)
    var_scv: dict[str, dict[str, dict]] = defaultdict(dict)
    var_stars: dict[str, int] = {}
    attr_types = Counter()
    with gzip.open(a.xml, "rt", encoding="utf-8") as fh, Pool(a.workers) as pool:
        batches = ((b, a.assembly) for b in _batched(_iter_clinvarset_chunks(fh), 200))
        for res in pool.imap_unordered(_work, batches, chunksize=4):
            for key, stars, scvs in res:
                var_stars[key] = max(stars, var_stars.get(key, 0))
                for s in scvs:
                    attr_types.update(s["attr_types"])
                    var_scv[key][s["acc"] or f"_{id(s)}"] = s

    def agg(scvs, plp_only=False, no_ep=False):
        sel = [s for s in scvs if (not plp_only or s["plp"]) and (not no_ep or not s["ep"])]
        return {
            "n_scv": len(sel),
            "fam": sum(s["fam"] for s in sel),
            "seg_yes": sum(s["seg_yes"] for s in sel),
            "seg_no": sum(s["seg_no"] for s in sel),
            "n_fam_seg": sum(s["n_fam_seg"] for s in sel),
            "text_seg": sum(s["text_seg"] for s in sel),
        }

    rows = []
    for key, d in var_scv.items():
        scvs = list(d.values())
        al, pl, pn = agg(scvs), agg(scvs, plp_only=True), agg(scvs, plp_only=True, no_ep=True)
        e = erepo.get(key)
        rows.append({
            "variant_id": key, "stars": var_stars.get(key, 0),
            "in_erepo": int(e is not None),
            "erepo_gene": e["gene"] if e else "", "erepo_assertion": e["assertion"] if e else "",
            "erepo_pp1": int(bool(e and e["pp1"])),
            **{f"all_{k}": v for k, v in al.items()},
            **{f"plp_{k}": v for k, v in pl.items()},
            **{f"plp_nonEP_{k}": v for k, v in pn.items()},
        })

    def pct(n, d):
        return f"{n:,} / {d:,} ({100 * n / d:.2f}%)" if d else f"{n} / 0"

    lines = [f"ClinVar FamilyData survey — {a.xml.name} ({a.assembly})", ""]
    for label, sub in (("All ClinVar variants with coordinates", rows),
                       ("ClinVar variants with a P/LP SCV", [r for r in rows if r["plp_n_scv"] > 0]),
                       ("eRepo variants found in ClinVar", [r for r in rows if r["in_erepo"]])):
        n = len(sub)
        lines.append(f"## {label}: n = {n:,}")
        lines.append("  any SCV with FamilyData/family attribute : " + pct(sum(r["all_fam"] > 0 for r in sub), n))
        lines.append("  P/LP SCV with FamilyData                 : " + pct(sum(r["plp_fam"] > 0 for r in sub), n))
        lines.append("  P/LP SCV, non-expert-panel, FamilyData   : " + pct(sum(r["plp_nonEP_fam"] > 0 for r in sub), n))
        lines.append("  P/LP SCV with SegregationObserved=yes    : " + pct(sum(r["plp_seg_yes"] > 0 for r in sub), n))
        lines.append("  P/LP SCV with >=1 family with segregation: " + pct(sum(r["plp_n_fam_seg"] > 0 for r in sub), n))
        lines.append("    ... non-expert-panel only              : " + pct(sum(r["plp_nonEP_seg_yes"] > 0 or r["plp_nonEP_n_fam_seg"] > 0 for r in sub), n))
        lines.append("  P/LP SCV, non-EP, 'segregat' in text     : " + pct(sum(r["plp_nonEP_text_seg"] > 0 for r in sub), n))
        lines.append("")
    er = [r for r in rows if r["in_erepo"]]
    if er:
        pp1 = [r for r in er if r["erepo_pp1"]]
        no = [r for r in er if not r["erepo_pp1"]]

        def structured(r):  # structured segregation evidence from non-EP P/LP SCVs
            return r["plp_nonEP_seg_yes"] > 0 or r["plp_nonEP_n_fam_seg"] > 0
        lines.append(f"## eRepo PP1 (VCEP applied PP1: {len(pp1):,}; not applied: {len(no):,})")
        lines.append("  structured segregation (non-EP P/LP), eRepo PP1     : " + pct(sum(structured(r) for r in pp1), len(pp1)))
        lines.append("  structured segregation (non-EP P/LP), eRepo not PP1 : " + pct(sum(structured(r) for r in no), len(no)))
        lines.append("  text 'segregat' (non-EP P/LP), eRepo PP1            : " + pct(sum(r["plp_nonEP_text_seg"] > 0 for r in pp1), len(pp1)))
        lines.append("  text 'segregat' (non-EP P/LP), eRepo not PP1        : " + pct(sum(r["plp_nonEP_text_seg"] > 0 for r in no), len(no)))
        lines.append(f"  eRepo variants not found in this ClinVar XML         : {len(erepo) - len(er):,}")
        lines.append("")
    lines.append("## ObservedData attribute types mentioning family/segregation/meiosis/pedigree")
    for t, c in attr_types.most_common(30):
        lines.append(f"  {c:>9,}  {t}")
    if not attr_types:
        lines.append("  (none)")

    Path(f"{a.out}_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    keep = [r for r in rows if r["in_erepo"] or r["all_fam"] or r["all_text_seg"]]
    with open(f"{a.out}_variants.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(keep, key=lambda r: r["variant_id"]))
    print("\n".join(lines))
    print(f"wrote {a.out}_summary.txt and {a.out}_variants.tsv ({len(keep):,} variants)")


if __name__ == "__main__":
    main()
