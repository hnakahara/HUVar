"""CIViC clinical evidence (Predictive / Diagnostic / Prognostic items and
assertions) for the AMP/ASCO/CAP 2017 tiers.

``build_civic_clinical`` keeps accepted records of single-variant molecular profiles
for gene variants and classifies how each CIViC variant is matched:

  exact   protein change (V600E, E746_A750del; HGVS numbering preferred)
  codon   any change at a residue (V600, G12) - also G12/G13 lists
  exon    exon-level events (Exon 19 Deletion, Exon 20 Insertion, Exon 11 Mutation,
          Exon 14 Skipping Mutation)
  gene    gene-level categories: 'mutation' (any protein-altering variant),
          'lof' (Loss-of-function, Deleterious / Truncating Mutation),
          'activating' (Activating Mutation, Gain-of-function)

Exon- and gene-level 'mutation'/'activating' evidence applies only to variants
classified Oncogenic / Likely oncogenic (``match_records``); 'lof' to null variants.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional

from acmg_classifier.somatic.civic import _SHORT_P, _protein_keys, _rows

CLINICAL_TYPES = ("Predictive", "Diagnostic", "Prognostic")
HEADER = ["gene", "match_type", "key", "variant_name", "record", "record_id", "evidence_type",
          "direction", "significance", "level", "rating", "doid", "disease", "therapies",
          "citation_id", "amp_category", "nccn_guideline", "regulatory_approval",
          "variant_origin"]

_CODON = re.compile(r"^([A-Z])(\d+)$")
_EXON = re.compile(r"^(?:rare )?exon (\d+)(?:-(\d+))? (deletion|insertion|mutation|frameshift|"
                   r"skipping mutation|skipping)$", re.I)
_GENE_LEVEL = {
    "mutation": "mutation", "mutations": "mutation", "alteration": "mutation",
    "loss-of-function": "lof", "loss of function": "lof", "deleterious mutation": "lof",
    "truncating mutation": "lof", "frameshift truncation": "lof", "inactivating mutation": "lof",
    "activating mutation": "activating", "gain-of-function": "activating",
    "gain of function": "activating",
}


def _normalize_short(name: str) -> str:
    return re.sub(r"(DEL|DUP|INS)", lambda m: m.group(1).lower(), name)


def classify_variant(name: str, hgvs: str) -> list[tuple[str, str]]:
    """[(match_type, key)] for a CIViC gene variant name; [] when not applicable to
    SNVs/indels (expression, amplification, fusions, ...)."""
    n = (name or "").strip()
    short = _normalize_short(n)
    if _SHORT_P.match(short):
        return [("exact", k) for k in _protein_keys(short, hgvs)]
    m = _CODON.match(n)
    if m:
        return [("codon", n)]
    if "/" in n:  # G12/G13, Q157P/R
        parts = n.split("/")
        if all(_CODON.match(p) for p in parts):
            return [("codon", p) for p in parts]
        m = re.match(r"^([A-Z]\d+)([A-Z])$", parts[0])
        if m and all(re.fullmatch(r"[A-Z]", p) for p in parts[1:]):
            return [("exact", m.group(1) + a) for a in [m.group(2)] + parts[1:]]
        return []
    m = _EXON.match(n)
    if m:
        lo, hi = int(m.group(1)), int(m.group(2) or m.group(1))
        kind = m.group(3).lower().replace("skipping mutation", "skipping")
        return [("exon", f"{e}:{kind}") for e in range(lo, hi + 1)]
    g = _GENE_LEVEL.get(n.lower())
    return [("gene", g)] if g else []


def build_civic_clinical(evidence_tsv: Path, variants_tsv: Path, assertions_tsv: Optional[Path],
                         out: Path, release: str = "") -> dict:
    mp: dict[str, list[tuple[str, str, str, str]]] = {}
    for r in _rows(variants_tsv):
        if (r.get("feature_type") or "Gene") != "Gene":
            continue
        mpid = (r.get("single_variant_molecular_profile_id") or "").strip()
        gene = (r.get("gene") or r.get("feature_name") or "").strip()
        if not mpid or not gene:
            continue
        hits = classify_variant(r.get("variant", ""), r.get("hgvs_descriptions", ""))
        if hits:
            mp[mpid] = [(gene, t, k, r.get("variant", "")) for t, k in hits]
    rows, n_ev, n_as = [], 0, 0
    for r in _rows(evidence_tsv):
        if (r.get("evidence_type") not in CLINICAL_TYPES
                or (r.get("evidence_status") or "accepted") != "accepted"
                or r.get("variant_origin") in ("Rare Germline", "Common Germline")):
            continue
        for gene, t, k, vname in mp.get((r.get("molecular_profile_id") or "").strip(), []):
            rows.append([gene, t, k, vname, "evidence", r.get("evidence_id", ""),
                         r["evidence_type"], r.get("evidence_direction", ""),
                         r.get("significance", ""), r.get("evidence_level", ""),
                         r.get("rating", ""), r.get("doid", ""), r.get("disease", ""),
                         r.get("therapies", ""), r.get("citation_id", ""), "", "", "",
                         r.get("variant_origin", "")])
            n_ev += 1
    if assertions_tsv is not None and assertions_tsv.exists():
        for r in _rows(assertions_tsv):
            if r.get("assertion_type") not in CLINICAL_TYPES:
                continue
            for gene, t, k, vname in mp.get((r.get("molecular_profile_id") or "").strip(), []):
                rows.append([gene, t, k, vname, "assertion", r.get("assertion_id", ""),
                             r["assertion_type"], r.get("assertion_direction", ""),
                             r.get("significance", ""), "", "", r.get("doid", ""),
                             r.get("disease", ""), r.get("therapies", ""), "",
                             r.get("amp_category", ""), r.get("nccn_guideline", ""),
                             r.get("regulatory_approval", ""), ""])
                n_as += 1
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# CIViC (CC0) predictive / diagnostic / prognostic evidence and assertions "
                 f"for AMP tiers; release {release or 'unknown'}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(sorted(rows, key=lambda x: (x[0], x[1], x[2], x[4], x[5])))
    return {"evidence_items": n_ev, "assertions": n_as}


@lru_cache(maxsize=4)
def load_clinical(path: Path) -> dict[str, list[dict]]:
    """gene -> records."""
    out: dict[str, list[dict]] = {}
    if path.exists():
        with path.open(encoding="utf-8") as fh:
            for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
                out.setdefault(r["gene"], []).append(r)
    return out


@dataclass(frozen=True)
class VariantFacts:
    gene: str
    protein_change: Optional[str]     # V600E / E746_A750del
    ref_aa: Optional[str]
    position: Optional[int]
    consequence: str                  # VEP consequence value
    exon: Optional[int]
    intron: Optional[int]
    is_null: bool
    oncogenic: bool                   # Oncogenic / Likely oncogenic


_EXON_KIND = {
    "deletion": {"inframe_deletion"},
    "insertion": {"inframe_insertion"},
    "frameshift": {"frameshift_variant"},
}


def match_records(db: dict, v: VariantFacts) -> list[tuple[dict, str]]:
    """CIViC records that apply to the variant: [(record, match description)]."""
    out = []
    for r in db.get(v.gene, []):
        t, k = r["match_type"], r["key"]
        if t == "exact":
            if v.protein_change and k == v.protein_change:
                out.append((r, "exact"))
        elif t == "codon":
            m = _CODON.match(k)
            if (m and v.position == int(m.group(2)) and v.ref_aa == m.group(1)
                    and v.consequence in ("missense_variant", "inframe_deletion",
                                          "inframe_insertion")):
                out.append((r, f"residue {k}"))
        elif t == "exon":
            num, kind = k.split(":")
            num = int(num)
            if kind == "skipping":
                ok = (v.consequence.startswith("splice") or v.consequence == "inframe_deletion") and (
                    v.exon == num or v.intron in (num - 1, num))
            elif kind == "mutation":
                ok = v.exon == num and v.consequence != "synonymous_variant" and v.oncogenic
            else:
                ok = v.exon == num and v.consequence in _EXON_KIND.get(kind, set())
            if ok:
                out.append((r, f"exon {num} {kind}"))
        elif t == "gene":
            if k == "lof" and v.is_null:
                out.append((r, "gene-level loss-of-function"))
            elif k in ("mutation", "activating") and v.oncogenic and v.consequence not in (
                    "synonymous_variant", "intron_variant", "5_prime_UTR_variant",
                    "3_prime_UTR_variant", "upstream_gene_variant", "downstream_gene_variant"):
                out.append((r, f"gene-level {k} (variant classified oncogenic)"))
    return out
