"""CIViC (CC0) functional evidence for the somatic OS2 / SBS2 criteria.

Two kinds of CIViC records are used, both for single-variant molecular profiles:

* accepted **Oncogenic assertions** whose ClinGen/CGC/VICC codes include OS2 or
  SBS2 (CIViC's own SOP-based curation of the variant);
* accepted **Functional evidence items** (significance Gain of Function,
  Neomorphic, Loss of Function, Dominant Negative, Unaltered Function; direction
  Supports). OS2/SBS2 need "well-established" studies, so evidence items count
  only with a trust rating >= ``somatic_civic_min_rating`` from >=
  ``somatic_civic_min_sources`` distinct publications, in the same direction.

Records are matched by gene and protein change in CIViC's short notation
(``V600E``, ``E746_A750del``), i.e. any nucleotide change giving that protein
change. ``build_civic_table`` turns the CIViC monthly/nightly TSV downloads
(``*-AcceptedClinicalEvidenceSummaries.tsv``, ``*-VariantSummaries.tsv``,
``*-AcceptedAssertionSummaries.tsv``) into ``data/shared/civic_somatic.tsv``.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional

csv.field_size_limit(1 << 26)

FUNCTIONAL_SIGNIFICANCES = {"Gain of Function", "Neomorphic", "Loss of Function",
                            "Dominant Negative", "Unaltered Function"}
_ACTIVATING = {"Gain of Function", "Neomorphic"}
_INACTIVATING = {"Loss of Function", "Dominant Negative"}
_NEUTRAL = {"Unaltered Function"}
# CIViC short protein-change names: V600E, R248*, E746_A750del, A767_V769dup,
# L747_P753delinsS, V600_K601delinsE.
_SHORT_P = re.compile(r"^[A-Z*]\d+(?:_[A-Z*]\d+)?(?:[A-Z*]|del|dup|ins[A-Z*]+|delins[A-Z*]+)$")
_HGVS_P = re.compile(r"p\.\(?([A-Z][a-z]{2}\d+(?:_[A-Z][a-z]{2}\d+)?[A-Za-z*]*)\)?")

HEADER = ["gene", "protein_change", "record", "record_id", "variant_id", "molecular_profile_id",
          "evidence_type", "significance", "direction", "level", "rating", "codes",
          "citation_id", "disease", "variant_origin"]


def _protein_keys(variant_name: str, hgvs: str) -> list[str]:
    """Short protein changes for a CIViC variant: the HGVS p. change (HGVS numbering,
    e.g. H3-3A p.Lys28Met) first, then the CIViC name when it differs (legacy
    numbering such as K27M)."""
    from acmg_classifier.somatic.resources import short_protein_change
    keys: list[str] = []
    for h in (hgvs or "").split(","):
        m = _HGVS_P.search(h)
        if m:
            key = short_protein_change("p." + m.group(1))
            if key and _SHORT_P.match(key):
                keys.append(key)
                break
    name = (variant_name or "").strip()
    if _SHORT_P.match(name) and name not in keys:
        keys.append(name)
    return keys


def _rows(path: Path):
    with path.open(encoding="utf-8", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def build_civic_table(evidence_tsv: Path, variants_tsv: Path, assertions_tsv: Optional[Path],
                      out: Path, release: str = "") -> dict:
    """Compact table of the CIViC records usable for OS2 / SBS2 (see module doc)."""
    # single-variant MP id -> [(gene, protein change, variant id)]
    mp_variant: dict[str, list[tuple[str, str, str]]] = {}
    for r in _rows(variants_tsv):
        if (r.get("feature_type") or "Gene") not in ("Gene", ""):
            continue
        mp = (r.get("single_variant_molecular_profile_id") or "").strip()
        keys = _protein_keys(r.get("variant", ""), r.get("hgvs_descriptions", ""))
        gene = (r.get("gene") or r.get("feature_name") or "").strip()
        if mp and keys and gene:
            mp_variant[mp] = [(gene, k, r.get("variant_id", "")) for k in keys]
    out_rows = []
    n_ev = n_as = 0
    for r in _rows(evidence_tsv):
        if (r.get("evidence_type") != "Functional"
                or r.get("significance") not in FUNCTIONAL_SIGNIFICANCES
                or (r.get("evidence_status") or "accepted") != "accepted"
                or r.get("variant_origin") == "Common Germline"):
            continue
        hits = mp_variant.get((r.get("molecular_profile_id") or "").strip())
        if not hits:
            continue
        for gene, key, vid in hits:
            out_rows.append([gene, key, "evidence", r.get("evidence_id", ""), vid,
                             r.get("molecular_profile_id", ""), "Functional", r["significance"],
                             r.get("evidence_direction", ""), r.get("evidence_level", ""),
                             r.get("rating", ""), "", r.get("citation_id", ""),
                             r.get("disease", ""), r.get("variant_origin", "")])
        n_ev += 1
    if assertions_tsv is not None and assertions_tsv.exists():
        for r in _rows(assertions_tsv):
            codes = (r.get("clingen_codes") or "").replace(" ", "")
            if r.get("assertion_type") != "Oncogenic" or not re.search(r"\b(OS2|SBS2)\b", codes):
                continue
            hits = mp_variant.get((r.get("molecular_profile_id") or "").strip())
            if not hits:
                continue
            for gene, key, vid in hits:
                out_rows.append([gene, key, "assertion", r.get("assertion_id", ""), vid,
                                 r.get("molecular_profile_id", ""), "Oncogenic",
                                 r.get("significance", ""), r.get("assertion_direction", ""),
                                 r.get("amp_category", ""), "", codes, "", r.get("disease", ""),
                                 ""])
            n_as += 1
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# CIViC (CC0) functional evidence / oncogenic assertions for OS2-SBS2; "
                 f"release {release or 'unknown'}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(sorted(out_rows, key=lambda x: (x[0], x[1], x[2], x[3])))
    return {"evidence_items": n_ev, "assertions": n_as, "release": release}


@lru_cache(maxsize=4)
def load_civic(path: Path) -> dict[tuple[str, str], list[dict]]:
    out: dict[tuple[str, str], list[dict]] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            out.setdefault((r["gene"], r["protein_change"]), []).append(r)
    return out


@dataclass
class CivicCall:
    criterion: Optional[str]      # "OS2" | "SBS2" | None
    evidence: str


def _eid(r: dict) -> str:
    prefix = "AID" if r["record"] == "assertion" else "EID"
    return f"{prefix}{r['record_id']}"


def civic_call(records: list[dict], role: Optional[str], min_rating: int = 3,
               min_sources: int = 2) -> Optional[CivicCall]:
    """OS2 / SBS2 from the CIViC records of one protein change (None = no record).

    ``role`` is the gene role (TSG / oncogene / dual / unresolved / None)."""
    if not records:
        return None
    for r in records:  # an SOP-based CIViC assertion settles it
        if r["record"] != "assertion" or r.get("direction") != "Supports":
            continue
        crit = "SBS2" if re.search(r"\bSBS2\b", r["codes"]) else "OS2"
        return CivicCall(crit, f"CIViC Oncogenic assertion {_eid(r)} ({r['significance']}; "
                               f"ClinGen/CGC/VICC codes {r['codes']})")
    if role in ("oncogene", "unresolved"):
        onco = _ACTIVATING
    elif role == "TSG":
        onco = _INACTIVATING | {"Neomorphic"}
    elif role == "dual":
        onco = _ACTIVATING | _INACTIVATING
    else:
        return CivicCall(None, "CIViC functional evidence present, but the gene role is unknown "
                               "(not a canonical cancer gene) - not interpreted")

    def qualifying(sigs):
        rows = [r for r in records if r["record"] == "evidence" and r["significance"] in sigs
                and r.get("direction") == "Supports"]
        good = [r for r in rows if _int(r.get("rating")) >= min_rating]
        return rows, good, {r["citation_id"] for r in good if r.get("citation_id")}

    o_all, o_good, o_src = qualifying(onco)
    b_all, b_good, b_src = qualifying(_NEUTRAL)

    def describe(rows):
        return ", ".join(f"{_eid(r)} {r['significance']} (PMID {r['citation_id']}, rating "
                         f"{r.get('rating') or '-'})" for r in rows[:5])

    o_ok, b_ok = len(o_src) >= min_sources, len(b_src) >= min_sources
    if o_ok and b_ok:
        return CivicCall(None, f"CIViC functional evidence conflicts: {describe(o_good)}; "
                               f"{describe(b_good)}")
    if o_ok:
        return CivicCall("OS2", f"CIViC functional evidence, {len(o_src)} publications: "
                                f"{describe(o_good)}")
    if b_ok:
        return CivicCall("SBS2", f"CIViC functional evidence of unaltered function, {len(b_src)} "
                                 f"publications: {describe(b_good)}")
    seen = o_all + b_all
    if not seen:
        other = ", ".join(f"{_eid(r)} {r['significance']}" for r in records[:5])
        return CivicCall(None, f"CIViC functional evidence does not fit the {role} role: {other}")
    return CivicCall(None, f"CIViC functional evidence below the threshold (needs rating >= "
                           f"{min_rating} from >= {min_sources} publications): {describe(seen)}")


def _int(v) -> int:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0
