"""Look up the existing ClinGen VCEP classification of a variant.

Reads the bundled eRepo snapshot ``resources/<assembly>/erepo_classifications_<hg>.tsv``
(built by ``scripts/build_erepo_classifications.py`` and enriched by
``scripts/add_erepo_details.py``). When a VCEP has already classified the queried
variant, HUVar reports it so that users can rely on the expert classification
rather than on the automated result (Genome Medicine revision, Reviewer 3,
comment 20).
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
from typing import Optional

_REPO = Path(__file__).resolve().parents[3]


def _strip_chr(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def _table(assembly: str) -> Path:
    hg = "hg38" if assembly == "GRCh38" else "hg19"
    return _REPO / "resources" / assembly / f"erepo_classifications_{hg}.tsv"


@lru_cache(maxsize=4)
def _load(path: Path) -> tuple[dict[str, dict], str]:
    out: dict[str, dict] = {}
    snapshot = ""
    if not path.exists():
        return out, snapshot
    with path.open(encoding="utf-8") as fh:
        lines = []
        for ln in fh:
            if ln.startswith("#"):
                if "snapshot" in ln:
                    snapshot = ln.split("snapshot", 1)[1].strip(" )\n")
                continue
            lines.append(ln)
    for r in csv.DictReader(lines, delimiter="\t"):
        chrom, pos, ref, alt = r["variant_id"].split(":")
        out[f"{_strip_chr(chrom)}:{pos}:{ref}:{alt}"] = r
    return out, snapshot


def lookup(assembly: str, chrom: str, pos: int, ref: str, alt: str) -> Optional[dict]:
    """The VCEP classification of the variant, or None.

    Keys: ``assertion``, ``met_codes``, ``conflict`` (bool), ``gene``,
    ``hgvs_c``, ``disease``, ``expert_panel``, ``cspec_url``, ``approval_date``,
    ``erepo_url``, ``snapshot``. Multi-record fields are ``|``-joined.
    """
    table, snapshot = _load(_table(assembly))
    r = table.get(f"{_strip_chr(str(chrom))}:{pos}:{ref}:{alt}")
    if r is None:
        return None
    return {
        "assertion": r.get("assertion", ""),
        "met_codes": r.get("met_codes", ""),
        "conflict": (r.get("conflict") or "0").strip() == "1",
        "gene": r.get("gene", ""),
        "transcript": r.get("transcript", ""),
        "hgvs_c": r.get("hgvs_c", ""),
        "disease": r.get("disease", ""),
        "expert_panel": r.get("expert_panel", ""),
        "cspec_url": r.get("cspec_url", ""),
        "approval_date": r.get("approval_date", ""),
        "erepo_url": r.get("erepo_url", ""),
        "snapshot": snapshot,
    }


def vcep_warning(assembly: str, chrom: str, pos: int, ref: str, alt: str) -> list[str]:
    """``["VCEP_CLASSIFIED: ..."]`` when a VCEP has classified the variant."""
    rec = lookup(assembly, chrom, pos, ref, alt)
    if rec is None:
        return []
    panel = rec["expert_panel"] or "a ClinGen VCEP"
    date = f" on {rec['approval_date']}" if rec["approval_date"] else ""
    url = rec["erepo_url"].split("|")[0] if rec["erepo_url"] else ""
    return [f"VCEP_CLASSIFIED: {panel} classified this variant as {rec['assertion']}{date}"
            f"{' (' + url + ')' if url else ''}; the expert classification should "
            "generally take precedence over this automated result"]
