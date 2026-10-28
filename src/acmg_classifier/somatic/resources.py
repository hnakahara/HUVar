"""Bundled lookup tables for somatic oncogenicity: cancer gene roles (NCG) and
cancerhotspots.org v2 hotspots."""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Optional

_SHARED = Path(__file__).resolve().parents[3] / "resources" / "shared"

AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
    "Sec": "U",
}


def resolve(path, name: str) -> Path:
    """Configured path when it exists, else the packaged resources copy."""
    if isinstance(path, (str, Path)) and Path(path).exists():
        return Path(path)
    return _SHARED / name


def _rows(path: Path):
    with path.open(encoding="utf-8") as fh:
        yield from csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")


@dataclass(frozen=True)
class GeneRole:
    gene: str
    role: str            # TSG | oncogene | dual | unresolved
    source: str          # NCG | override
    note: str = ""

    @property
    def is_tsg(self) -> bool:
        """Bona fide tumour suppressor for OVS1: TSG role only (dual / unresolved excluded)."""
        return self.role == "TSG"

    @property
    def is_cancer_gene(self) -> bool:
        """Known oncogene or tumour suppressor for OM2 (any canonical role)."""
        return self.role in ("TSG", "oncogene", "dual")

    def describe(self) -> str:
        src = "NCG canonical driver" if self.source == "NCG" else "curated override"
        s = f"{self.gene}: {self.role} ({src})"
        return s + (f"; {self.note}" if self.note else "")


@lru_cache(maxsize=4)
def load_gene_roles(path: Path) -> dict[str, GeneRole]:
    if not path.exists():
        return {}
    return {r["gene"]: GeneRole(r["gene"], r["role"], r["source"], r.get("note") or "")
            for r in _rows(path)}


@dataclass
class HotspotSite:
    gene: str
    kind: str                     # snv | indel
    position: str                 # "600" or "229-292"
    ref_aa: str
    position_count: int
    variants: dict[str, int] = field(default_factory=dict)

    def describe(self) -> str:
        unit = "residue" if self.kind == "snv" else "indel region"
        name = f"{self.ref_aa}{self.position}" if self.kind == "snv" else self.position
        top = ", ".join(f"{name + k if self.kind == 'snv' else k} {v}"
                        for k, v in sorted(self.variants.items(), key=lambda x: -x[1])[:4])
        return (f"cancerhotspots.org v2 {self.gene} {unit} {name}: {self.position_count} "
                f"samples ({top})")


@lru_cache(maxsize=4)
def load_hotspots(path: Path) -> dict[tuple[str, str], HotspotSite]:
    """(gene, position) -> site for single residues; (gene, "indel:<range>") for
    in-frame indel regions."""
    out: dict[tuple[str, str], HotspotSite] = {}
    if not path.exists():
        return out
    for r in _rows(path):
        key = (r["gene"], r["position"] if r["type"] == "snv" else f"indel:{r['position']}")
        site = out.get(key)
        if site is None:
            site = out[key] = HotspotSite(r["gene"], r["type"], r["position"], r["ref_aa"],
                                          int(r["position_count"]))
        site.variants[r["variant"]] = site.variants.get(r["variant"], 0) + int(r["variant_count"])
    return out


def residue_site(sites: dict, gene: str, position: Optional[int]) -> Optional[HotspotSite]:
    if position is None:
        return None
    return sites.get((gene, str(position)))


def indel_site(sites: dict, gene: str, start: Optional[int],
               end: Optional[int] = None) -> Optional[HotspotSite]:
    """In-frame indel hotspot region containing residues start..end."""
    if start is None:
        return None
    end = end if end is not None else start
    for (g, key), site in sites.items():
        if g != gene or not key.startswith("indel:"):
            continue
        lo, _, hi = site.position.partition("-")
        try:
            lo_i, hi_i = int(lo), int(hi or lo)
        except ValueError:
            continue
        if lo_i <= start and end <= hi_i:
            return site
    return None


def short_protein_change(hgvs_p: Optional[str]) -> Optional[str]:
    """'NP_x:p.Glu286_Glu287del' -> 'E286_E287del' (cancerhotspots indel notation)."""
    if not hgvs_p:
        return None
    p = hgvs_p.split(":")[-1].replace("p.", "").strip("()")
    out, i = [], 0
    while i < len(p):
        tri = p[i:i + 3]
        if tri in AA3_TO_AA1:
            out.append(AA3_TO_AA1[tri])
            i += 3
        else:
            out.append(p[i])
            i += 1
    return "".join(out)
