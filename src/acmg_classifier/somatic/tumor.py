"""Tumour types: OncoTree code -> Disease Ontology IDs, and DO hierarchy matching
(``oncotree_doid.tsv`` / ``doid_hierarchy.tsv``, scripts/build_tumor_ontology.py)."""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Optional

from acmg_classifier.somatic.resources import resolve

CANCER = "DOID:162"


def _rows(path: Path):
    with path.open(encoding="utf-8") as fh:
        yield from csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")


def norm_doid(v) -> Optional[str]:
    s = str(v or "").strip()
    if not s:
        return None
    return s if s.startswith("DOID:") else f"DOID:{s}"


@dataclass
class DiseaseOntology:
    parents: dict[str, list[str]]
    names: dict[str, str]
    by_name: dict[str, str] = field(default_factory=dict)

    def ancestors(self, doid: str) -> set[str]:
        seen: set[str] = set()
        stack = [doid]
        while stack:
            d = stack.pop()
            for p in self.parents.get(d, []):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        return seen

    def name(self, doid: Optional[str]) -> str:
        return self.names.get(doid or "", doid or "")

    def lookup_name(self, text: str) -> Optional[str]:
        return self.by_name.get((text or "").strip().lower())


@lru_cache(maxsize=2)
def load_do(path: Path) -> DiseaseOntology:
    parents, names, by_name = {}, {}, {}
    if path.exists():
        for r in _rows(path):
            parents[r["doid"]] = [p for p in (r.get("parents") or "").split(",") if p]
            names[r["doid"]] = r["name"]
            for nm in [r["name"]] + [s for s in (r.get("synonyms") or "").split("|") if s]:
                by_name.setdefault(nm.lower(), r["doid"])
    return DiseaseOntology(parents, names, by_name)


@dataclass(frozen=True)
class TumorType:
    code: str
    name: str
    doids: tuple[str, ...]
    method: str

    @property
    def direct(self) -> bool:
        """Mapped from the node itself (not from a parent / tissue root)."""
        return self.method in ("nci", "umls", "name", "curated") and CANCER not in self.doids

    def describe(self, do: DiseaseOntology) -> str:
        how = "" if self.method in ("nci", "umls", "name") else f", via {self.method}"
        return (f"{self.name} (OncoTree {self.code}; "
                f"{', '.join(f'{do.name(d)} {d}' for d in self.doids)}{how})")


@lru_cache(maxsize=2)
def load_oncotree(path: Path) -> dict[str, TumorType]:
    out = {}
    if path.exists():
        for r in _rows(path):
            ds = tuple(d for d in (r.get("doids") or "").split(",") if d)
            out[r["oncotree_code"].upper()] = TumorType(r["oncotree_code"], r["name"], ds,
                                                        r.get("method") or "")
    return out


def tumor_type(code: Optional[str], cfg=None) -> Optional[TumorType]:
    if not code:
        return None
    path = resolve(getattr(cfg, "oncotree_doid_tsv", None), "oncotree_doid.tsv")
    tt = load_oncotree(path).get(code.strip().upper())
    if tt is None:
        raise ValueError(f"unknown OncoTree code {code!r} (see resources/shared/oncotree_doid.tsv)")
    return tt


def disease_ontology(cfg=None) -> DiseaseOntology:
    return load_do(resolve(getattr(cfg, "doid_hierarchy_tsv", None), "doid_hierarchy.tsv"))


def relation(tt: TumorType, evidence_doid: Optional[str], do: DiseaseOntology) -> str:
    """'same' (the tumour type or a more general disease, e.g. 'cancer' for a
    tumour-agnostic approval), 'subtype' (a more specific disease of the tumour
    type), 'other' or 'unknown' (no DOID on the evidence)."""
    e = norm_doid(evidence_doid)
    if not e:
        return "unknown"
    for d in tt.doids:
        if e == d or e in do.ancestors(d):
            return "same"
    for d in tt.doids:
        if d != CANCER and d in do.ancestors(e):
            return "subtype"
    return "other"
