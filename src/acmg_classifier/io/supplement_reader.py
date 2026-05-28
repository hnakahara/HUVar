from __future__ import annotations
import csv
from pathlib import Path

from acmg_classifier.exceptions import SupplementParseError
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.supplement import SupplementEntry

# Required TSV columns
_REQUIRED = {"variant_id", "criterion", "strength", "evidence"}


def read_supplement(tsv_path: Path) -> dict[str, list[SupplementEntry]]:
    """Parse manual evidence TSV into a dict keyed by variant_id."""
    entries: dict[str, list[SupplementEntry]] = {}
    with tsv_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or not _REQUIRED.issubset(reader.fieldnames):
            missing = _REQUIRED - set(reader.fieldnames or [])
            raise SupplementParseError(f"Supplement TSV missing columns: {missing}")
        for lineno, row in enumerate(reader, start=2):
            try:
                entry = SupplementEntry(
                    variant_id=row["variant_id"].strip(),
                    criterion=ACMGCriterion(row["criterion"].strip()),
                    strength=CriterionStrength(row["strength"].strip()),
                    evidence=row["evidence"].strip(),
                )
            except (KeyError, ValueError) as exc:
                raise SupplementParseError(
                    f"Supplement TSV line {lineno}: {exc}"
                ) from exc
            entries.setdefault(entry.variant_id, []).append(entry)
    return entries
