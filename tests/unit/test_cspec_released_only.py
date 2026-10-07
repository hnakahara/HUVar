"""Only Released ClinGen specifications may drive gene-specific criteria
(Genome Medicine revision, Reviewer 3, comment 24)."""
import csv
from pathlib import Path

_DP = Path(__file__).resolve().parents[2] / "resources" / "shared" / "disease_prevalence.tsv"


def test_every_gene_row_comes_from_a_released_spec():
    with _DP.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    bad = [(r["gene_symbol"], r["notes"]) for r in rows
           if not (r.get("notes") or "").split(";")[0].strip().endswith("Released")]
    assert bad == []
