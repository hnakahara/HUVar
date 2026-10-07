"""Add VCEP details to ``resources/<assembly>/erepo_classifications_<hg>.tsv``.

HUVar_app shows the existing ClinGen VCEP classification of a queried variant
(Genome Medicine revision, Reviewer 3, comment 20). The classification table
built by ``build_erepo_classifications.py`` carries the assertion and met codes;
this script adds, from the tab-separated eRepo export, the disease, the expert
panel, the panel's CSpec affiliation page, the approval date and the eRepo record
link. Records are matched on (RefSeq accession without version, HGVS c.) against
the ``HGVS Expressions`` column. Several eRepo records for one variant are joined
with ``|`` in the same order. Retracted records are dropped.

Usage::

    python scripts/add_erepo_details.py --tabbed erepo-tabbed.tsv \
        resources/GRCh38/erepo_classifications_hg38.tsv \
        resources/GRCh37/erepo_classifications_hg19.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

DETAIL_COLS = ["disease", "expert_panel", "cspec_url", "approval_date", "erepo_url"]
_NM = re.compile(r"(NM_\d+)\.\d+:(c\..+)$")


def load_tabbed(path: Path) -> dict[tuple[str, str], list[dict]]:
    csv.field_size_limit(10**9)
    idx: dict[tuple[str, str], dict[str, dict]] = defaultdict(dict)
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("Retracted") or "").strip().lower() == "true":
                continue
            for h in (r.get("HGVS Expressions") or "").split(","):
                m = _NM.match(h.strip())
                if m:
                    idx[(m.group(1), m.group(2))][r["Uuid"]] = r
    return {k: list(v.values()) for k, v in idx.items()}


def enrich(table: Path, idx) -> tuple[int, int]:
    lines = table.read_text(encoding="utf-8").splitlines()
    head = [ln for ln in lines if ln.startswith("#")]
    body = [ln for ln in lines if not ln.startswith("#")]
    rows = list(csv.DictReader(body, delimiter="\t"))
    cols = [c for c in csv.reader([body[0]], delimiter="\t").__next__() if c not in DETAIL_COLS]
    n_hit = 0
    for r in rows:
        recs = idx.get((r["transcript"].split(".")[0], r["hgvs_c"]), [])
        recs.sort(key=lambda x: (x.get("Approval Date") or "", x.get("Uuid") or ""))
        if recs:
            n_hit += 1
        r["disease"] = "|".join((x.get("Disease") or "").strip() for x in recs)
        r["expert_panel"] = "|".join((x.get("Expert Panel") or "").strip() for x in recs)
        r["cspec_url"] = "|".join((x.get("Guideline") or "").strip() for x in recs)
        r["approval_date"] = "|".join((x.get("Approval Date") or "").strip() for x in recs)
        r["erepo_url"] = "|".join((x.get("Evidence Repo Link") or "").strip() for x in recs)
    with table.open("w", encoding="utf-8", newline="\n") as fh:
        for h in head:
            fh.write(h + "\n")
        w = csv.DictWriter(fh, fieldnames=cols + DETAIL_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return len(rows), n_hit


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--tabbed", type=Path, required=True, help="eRepo tab-separated export")
    ap.add_argument("tables", type=Path, nargs="+")
    a = ap.parse_args()
    idx = load_tabbed(a.tabbed)
    for t in a.tables:
        n, hit = enrich(t, idx)
        print(f"{t}: {hit}/{n} variants with eRepo details")


if __name__ == "__main__":
    main()
