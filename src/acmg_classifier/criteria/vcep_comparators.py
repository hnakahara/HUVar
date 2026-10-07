"""Restrict PS1/PM5 comparators to VCEP-classified variants where the VCEP requires it.

For the genes listed in ``resources/shared/vcep_comparator_rules.tsv`` (built by
``scripts/build_vcep_comparator_rules.py`` from the CSpec texts), the previously
established pathogenic variant must have been classified by the VCEP. A ClinVar
comparator is accepted only when

* it is a ClinGen expert-panel record in ClinVar (review status >= 3 stars), or
* it is classified P/LP in the eRepo snapshot
  (``resources/<assembly>/erepo_classifications_<hg>.tsv``),

and, when the VCEP adds conditions, it was not itself classified using PS1/PM5
(``exclude_codes``, checked on the eRepo met codes) and its pathogenicity is at
the protein level, i.e. it was not classified through PVS1 (``protein_level``).
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple, Optional

from acmg_classifier.models.annotation import ClinVarRecord

# parents: [0]=criteria [1]=acmg_classifier [2]=src [3]=repo
_REPO = Path(__file__).resolve().parents[3]
_RULES = _REPO / "resources" / "shared" / "vcep_comparator_rules.tsv"


class Rule(NamedTuple):
    exclude_codes: frozenset[str]
    protein_level: bool


def _strip_chr(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


@lru_cache(maxsize=2)
def load_rules(path: Path = _RULES) -> dict[tuple[str, str], Rule]:
    out: dict[tuple[str, str], Rule] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("comparator") or "").strip() != "vcep":
                continue
            out[(r["gene"].strip(), r["criterion"].strip())] = Rule(
                frozenset(c.strip() for c in (r.get("exclude_codes") or "").split(",") if c.strip()),
                (r.get("protein_level") or "0").strip() == "1",
            )
    return out


@lru_cache(maxsize=4)
def load_erepo(assembly: str) -> dict[str, tuple[str, frozenset[str]]]:
    """``{chrom-stripped variant key: (assertion, base met codes)}`` for the snapshot."""
    hg = "hg38" if assembly == "GRCh38" else "hg19"
    path = _REPO / "resources" / assembly / f"erepo_classifications_{hg}.tsv"
    out: dict[str, tuple[str, frozenset[str]]] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8") as fh:
        rows = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for r in rows:
            chrom, pos, ref, alt = r["variant_id"].split(":")
            codes = frozenset(c.split("_")[0].strip() for c in (r.get("met_codes") or "").split(",") if c.strip())
            out[f"{_strip_chr(chrom)}:{pos}:{ref}:{alt}"] = (r["assertion"], codes)
    return out


def requires_vcep_comparator(gene: Optional[str], criterion: str) -> Optional[Rule]:
    if not gene:
        return None
    return load_rules().get((gene, criterion))


def filter_comparators(
    gene: Optional[str], criterion: str, hits: list[ClinVarRecord], assembly: str,
) -> tuple[list[ClinVarRecord], Optional[str]]:
    """Return ``(accepted_hits, note)``; ``note`` is None when no VCEP rule applies."""
    rule = requires_vcep_comparator(gene, criterion)
    if rule is None:
        return hits, None
    erepo = load_erepo(assembly)
    kept: list[ClinVarRecord] = []
    for h in hits:
        key = (f"{_strip_chr(str(h.chrom))}:{h.pos}:{h.ref}:{h.alt}"
               if h.chrom is not None and h.pos is not None else None)
        rec = erepo.get(key) if key else None
        if rec is not None:
            assertion, codes = rec
            if "Pathogenic" not in assertion:
                continue  # VCEP did not classify it P/LP
            if rule.exclude_codes & codes:
                continue  # classified using PS1/PM5 itself
            if rule.protein_level and "PVS1" in codes:
                continue  # pathogenic through a splicing / null effect
            kept.append(h)
        elif (h.star_rating or 0) >= 3:
            kept.append(h)  # ClinVar expert-panel record not in the eRepo snapshot
    note = (f"{gene} VCEP requires a VCEP-classified comparator: "
            f"{len(kept)}/{len(hits)} ClinVar comparator(s) accepted")
    return kept, note
