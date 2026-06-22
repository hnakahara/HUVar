"""SCN paralogue amino-acid alignment for the PS1 analogous-residue route.

The SCN epilepsy VCEPs (GN067-070) grant PS1 from the same amino-acid change at
the *analogous* residue of a paralogue gene (SCN1A/SCN2A/SCN3A/SCN8A), per the
"Paralogous Gene Table". Unlike the RASopathy / HBA paralogues (identical residue
numbering), the SCN paralogues need a residue-correspondence map: this loads the
alignment table built by ``scripts/build_ps1_paralog_map.py`` into a
``(gene, position) -> {sibling_gene: sibling_position}`` lookup.
"""
from __future__ import annotations

import csv
from pathlib import Path


class PS1ParalogMap:
    def __init__(self, tsv_path: Path) -> None:
        self._genes: list[str] = []
        # (gene, pos) -> {sibling_gene: sibling_pos}
        self._map: dict[tuple[str, int], dict[str, int]] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        try:
            if not tsv_path.exists():
                return
            with tsv_path.open(encoding="utf-8") as fh:
                rows = list(csv.reader(fh, delimiter="\t"))
        except (OSError, TypeError, AttributeError):
            # Missing file or a non-Path (mocked) cfg attribute → no map.
            return
        if rows:
            self._genes = rows[0]
            for row in rows[1:]:
                present: dict[str, int] = {}
                for gene, cell in zip(self._genes, row):
                    cell = (cell or "").strip()
                    if cell:
                        try:
                            present[gene] = int(cell)
                        except ValueError:
                            pass
                for gene, pos in present.items():
                    sibs = {g: p for g, p in present.items() if g != gene}
                    if sibs:
                        self._map[(gene, pos)] = sibs

    def has_gene(self, gene: str | None) -> bool:
        return bool(gene) and gene in self._genes

    def analogs(self, gene: str | None, pos: int | None) -> dict[str, int]:
        """Analogous ``{sibling_gene: position}`` for *gene* residue *pos*, or {}
        when the gene/position is not in the alignment (or aligns to a gap)."""
        if not gene or pos is None:
            return {}
        return self._map.get((gene, pos), {})
