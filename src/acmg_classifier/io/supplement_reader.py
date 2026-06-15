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


# Friendly strength aliases for the inline `explain --evidence` form, in addition
# to the canonical CriterionStrength values (VeryStrong/Strong/Moderate/...).
_STRENGTH_ALIASES = {
    "very_strong": CriterionStrength.VERY_STRONG,
    "verystrong": CriterionStrength.VERY_STRONG,
    "vs": CriterionStrength.VERY_STRONG,
    "stand_alone": CriterionStrength.VERY_STRONG,   # BA1 stand-alone benign
    "standalone": CriterionStrength.VERY_STRONG,
    "sa": CriterionStrength.VERY_STRONG,
    "strong": CriterionStrength.STRONG,
    "three_point": CriterionStrength.THREE_POINT,
    "threepoint": CriterionStrength.THREE_POINT,
    "moderate": CriterionStrength.MODERATE,
    "mod": CriterionStrength.MODERATE,
    "supporting": CriterionStrength.SUPPORTING,
    "supp": CriterionStrength.SUPPORTING,
    "not_met": CriterionStrength.NOT_MET,
    "notmet": CriterionStrength.NOT_MET,
}


def _parse_strength(token: str) -> CriterionStrength:
    """CriterionStrength from a canonical value (case-insensitive) or a friendly
    alias (strong, very_strong, mod, supp, …)."""
    t = token.strip()
    alias = _STRENGTH_ALIASES.get(t.lower())
    if alias is not None:
        return alias
    # Canonical enum value, matched case-insensitively (e.g. "strong"->"Strong").
    for s in CriterionStrength:
        if s.value.lower() == t.lower():
            return s
    raise SupplementParseError(
        f"Unknown strength '{token}'. Use one of: very_strong, strong, "
        f"three_point, moderate, supporting (or a CriterionStrength value)."
    )


def parse_inline_evidence(specs: list[str], variant_id: str) -> list[SupplementEntry]:
    """Parse ``explain --evidence`` strings into SupplementEntry rows for one variant.

    Each spec is ``CRITERION:STRENGTH`` or ``CRITERION:STRENGTH:free text``, e.g.
    ``PS3:strong``, ``PM1:moderate:hotspot``, ``BA1:stand_alone:gnomAD AF 5%``.
    Criterion and strength are case-insensitive; the optional third field is the
    evidence note (may itself contain ':')."""
    out: list[SupplementEntry] = []
    for spec in specs:
        parts = spec.split(":", 2)
        if len(parts) < 2:
            raise SupplementParseError(
                f"Invalid --evidence '{spec}'. Expected CRITERION:STRENGTH[:note]."
            )
        crit_tok, strength_tok = parts[0].strip(), parts[1].strip()
        note = parts[2].strip() if len(parts) == 3 else "manual evidence (explain --evidence)"
        try:
            criterion = ACMGCriterion(crit_tok.upper())
        except ValueError as exc:
            raise SupplementParseError(
                f"Unknown criterion '{crit_tok}' in --evidence '{spec}'."
            ) from exc
        out.append(SupplementEntry(
            variant_id=variant_id,
            criterion=criterion,
            strength=_parse_strength(strength_tok),
            evidence=note,
            source="explain --evidence",
        ))
    return out
