"""Per-CSpec overlay of ``disease_prevalence.tsv`` for multi-disease genes.

A few genes (RYR1, ACTA1, VWF) carry several Released, clinically distinct
ClinGen CSpecs that the batch ``disease_prevalence.tsv`` collapses into one
conservative row. ``scripts/build_disease_thresholds.py --multispec-out`` emits a
side table ``disease_prevalence_multispec.tsv`` with one row per gene×CSpec.

Every per-gene threshold loader in this package keys strictly on the
``gene_symbol`` column, so a variant can be evaluated under a specific CSpec
WITHOUT any change to the evaluators: take the conservative base table, replace
that one gene's row with the chosen CSpec's row, and point
``Config.disease_prevalence_tsv_override`` at the result. This module builds that
overlay table and lists the CSpecs available for a gene.

The batch ``disease_prevalence.tsv`` and all CLI / batch behaviour are untouched.
"""
from __future__ import annotations

import csv
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, Optional

if TYPE_CHECKING:
    from acmg_classifier.config import Config


def _read_multispec(multispec_tsv: Path) -> dict[tuple[str, str], dict]:
    """``(gene_symbol, cspec_id) -> row`` from the multispec table (empty if the
    file is absent — the app then simply offers no CSpec switch)."""
    out: dict[tuple[str, str], dict] = {}
    if not multispec_tsv.exists():
        return out
    with multispec_tsv.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = (row.get("gene_symbol") or "").strip()
            cspec = (row.get("cspec_id") or "").strip()
            if gene and cspec:
                out[(gene, cspec)] = row
    return out


def available_cspecs(gene: Optional[str], multispec_tsv: Path | str) -> list[dict]:
    """CSpecs available for *gene*, in table order.

    Returns ``[{"cspec_id", "label", "source_gn"}, ...]`` — empty when the gene
    has no multispec entry or the table is missing (the common case; only the
    curated multi-disease genes return a non-empty list)."""
    multispec_tsv = Path(multispec_tsv)
    if not gene or not multispec_tsv.exists():
        return []
    out: list[dict] = []
    with multispec_tsv.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("gene_symbol") or "").strip() == gene:
                out.append({
                    "cspec_id": (row.get("cspec_id") or "").strip(),
                    "label": (row.get("disease_label") or "").strip(),
                    "source_gn": (row.get("source_gn") or "").strip(),
                })
    return out


def overlay_tsv(
    base_tsv: Path | str,
    multispec_tsv: Path | str,
    gene: str,
    cspec_id: str,
) -> Path:
    """Write a copy of *base_tsv* with *gene*'s row replaced by the *cspec_id*
    row from *multispec_tsv*, and return the temporary file's path.

    The overlay keeps the base table's exact column set and order, so it is a
    drop-in for ``Config.disease_prevalence_tsv_override``. The caller owns the
    returned file and should delete it when done (see :func:`overlaid_config`).
    Raises ``KeyError`` if the requested gene×CSpec is not in the multispec table.
    """
    base_tsv = Path(base_tsv)
    multispec_tsv = Path(multispec_tsv)
    repl = _read_multispec(multispec_tsv).get((gene, cspec_id))
    if repl is None:
        raise KeyError(f"no multispec row for {gene}/{cspec_id} in {multispec_tsv}")

    with base_tsv.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    # Restrict the replacement to the base table's columns (the multispec table
    # carries extra leading columns — cspec_id / disease_label / source_gn —
    # that must not leak into the overlay).
    new_row = {col: repl.get(col, "") for col in fieldnames}
    for i, r in enumerate(rows):
        if (r.get("gene_symbol") or "").strip() == gene:
            rows[i] = new_row
            break
    else:  # gene absent from the base table — append so the CSpec still applies
        rows.append(new_row)

    fd, tmp = tempfile.mkstemp(suffix=".tsv", prefix=f"cspec_{gene}_{cspec_id}_")
    with os.fdopen(fd, "w", newline="\n", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return Path(tmp)


@contextmanager
def overlaid_config(
    cfg: "Config",
    multispec_tsv: Path | str,
    gene: str,
    cspec_id: str,
) -> Iterator["Config"]:
    """Yield a copy of *cfg* whose ``disease_prevalence_tsv_override`` points at a
    per-CSpec overlay table, deleting the temporary file on exit.

    Usage::

        with overlaid_config(cfg, multispec, "RYR1", "malignant_hyperthermia") as c:
            result = classify_annotated(variant, ann, c)

    The base ``cfg.disease_prevalence_tsv`` is used as the overlay base, so this
    composes with any existing override."""
    tmp = overlay_tsv(cfg.disease_prevalence_tsv, multispec_tsv, gene, cspec_id)
    try:
        yield cfg.model_copy(update={"disease_prevalence_tsv_override": tmp})
    finally:
        try:
            tmp.unlink()
        except OSError:
            pass
