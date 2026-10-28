"""Laboratory-curated functional evidence for OS2 / SBS2
(``resources/shared/somatic_functional_curated.tsv``; a copy in
``data/shared/`` takes precedence). Rows apply only with status=approved."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class CuratedFunctional:
    gene: str
    protein_change: str
    criterion: str          # OS2 | SBS2
    pmids: str
    evidence: str
    verification: str
    status: str             # approved | draft | ...
    curator: str = ""
    reviewed: str = ""
    note: str = ""

    @property
    def approved(self) -> bool:
        return self.status.strip().lower() == "approved"

    def describe(self) -> str:
        pm = ", ".join(f"PMID {p.strip()}" for p in self.pmids.split(",") if p.strip())
        who = f"; reviewed {self.reviewed} by {self.curator}" if self.curator or self.reviewed else ""
        note = f"; {self.note}" if self.note else ""
        return (f"lab-curated functional evidence ({self.status}): {self.evidence} [{pm}]"
                f"{who}{note}")


@lru_cache(maxsize=4)
def load_curated(path: Path) -> dict[tuple[str, str], CuratedFunctional]:
    out: dict[tuple[str, str], CuratedFunctional] = {}
    if not path.exists():
        return out
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if ln.strip() and not ln.startswith("#")),
                                delimiter="\t"):
            crit = (r.get("criterion") or "").strip().upper()
            if crit not in ("OS2", "SBS2") or not r.get("gene") or not r.get("protein_change"):
                continue
            row = CuratedFunctional(
                gene=r["gene"].strip(), protein_change=r["protein_change"].strip(),
                criterion=crit, pmids=(r.get("pmids") or "").strip(),
                evidence=(r.get("evidence") or "").strip(),
                verification=(r.get("verification") or "").strip(),
                status=(r.get("status") or "draft").strip(),
                curator=(r.get("curator") or "").strip(),
                reviewed=(r.get("reviewed") or "").strip(),
                note=(r.get("note") or "").strip(),
            )
            out[(row.gene, row.protein_change)] = row
    return out


def lookup(db: dict, gene: Optional[str], key: Optional[str]) -> Optional[CuratedFunctional]:
    if not gene or not key:
        return None
    return db.get((gene, key))
