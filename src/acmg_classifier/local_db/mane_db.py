"""Gene -> MANE Select transcript map.

MANE Select is GRCh38-native: VEP populates its ``mane_select`` flag only in the
GRCh38 cache. On GRCh37 no transcript carries the flag, so the MANE-first
transcript selection in ``vep_runner`` cannot fire and a non-MANE isoform may be
chosen — shifting HGVS numbering and breaking VCEP criteria keyed to the MANE
codon range (e.g. PVS1 range). This map lets the annotation layer recover the
MANE-equivalent transcript by version-stripped RefSeq/Ensembl base accession.

TSV source (built by ``scripts/build_mane_map.py`` from the MANE GFF)::

    gene_symbol<TAB>refseq<TAB>ensembl
    PTEN<TAB>NM_000314.8<TAB>ENST00000371953.8
"""
from __future__ import annotations
import csv
from functools import lru_cache
from pathlib import Path

import structlog

log = structlog.get_logger()

# Packaged fallback: <repo>/resources/shared/mane_select.tsv. Resolved relative
# to this file so an editable (`pip install -e`) checkout works without staging
# the TSV into data_dir. parents: [0]=local_db [1]=acmg_classifier [2]=src [3]=repo.
_PACKAGED_TSV = Path(__file__).resolve().parents[3] / "resources" / "shared" / "mane_select.tsv"


def _base(acc: str | None) -> str:
    """Version-stripped accession (NM_000314.8 -> NM_000314)."""
    return (acc or "").split(".")[0]


@lru_cache(maxsize=8)
def load_mane_map(tsv_path: Path) -> dict[str, tuple[str, str]]:
    """Return ``{gene_symbol: (refseq_base, ensembl_base)}``.

    Reads ``tsv_path`` (``data_dir/shared/mane_select.tsv``) when present,
    otherwise the packaged resources copy. Empty dict if neither exists.
    """
    path = tsv_path if tsv_path.exists() else _PACKAGED_TSV
    if not path.exists():
        log.warning("mane_map_missing", path=str(tsv_path))
        return {}

    result: dict[str, tuple[str, str]] = {}
    try:
        with path.open(encoding="utf-8") as fh:
            reader = csv.DictReader(
                (line for line in fh if not line.lstrip().startswith("#")),
                delimiter="\t",
            )
            for row in reader:
                gene = (row.get("gene_symbol") or "").strip()
                if gene:
                    result[gene] = (_base(row.get("refseq")), _base(row.get("ensembl")))
    except Exception as exc:  # noqa: BLE001
        log.error("mane_map_error", error=str(exc))
        return {}

    log.info("mane_map_loaded", genes=len(result), path=str(path))
    return result
