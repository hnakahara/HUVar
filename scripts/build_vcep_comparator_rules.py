"""Build ``resources/shared/vcep_comparator_rules.tsv`` — VCEPs that require the
PS1/PM5 comparator variant to have been classified by the VCEP itself.

Several CSpec specifications state that the previously established pathogenic
variant used for PS1 or PM5 must be classified by the panel (or with the panel's
rules) — e.g. the SCID, InSiGHT, TP53, VHL, VWD, Peroxisomal, Monogenic Diabetes,
LGMD and Antibody Deficiencies VCEPs. Some add conditions: the comparator must
have been classified pathogenic *without using PS1/PM5* (CTLA4, PIK3R1), or be
pathogenic at the protein level and not through a splicing defect (InSiGHT MMR).

The table is mined from the PS1/PM5 code texts in ``cspec_summary.json`` with the
patterns below and should be reviewed manually after each CSpec update.

Columns: gene, criterion (PS1/PM5), comparator (vcep), exclude_codes
(comparators whose eRepo met codes include these are not used), protein_level
(1 = exclude comparators classified via PVS1, i.e. splicing), source_text.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

_VCEP_RE = re.compile(
    r"(classified|curated|assessed|determined|established)[^.]{0,60}(by|using rules from|per|"
    r"according to)[^.]{0,20}(the )?([A-Za-z0-9/\- ]{0,25})(VCEP|Expert Panel|EP\b|panel)|"
    r"(VCEP|expert panel)[- ](classified|approved|curated|assertion)|"
    r"(classified|curated|established|determined)\s+(by|according to|using)\s+(the\s+)?"
    r"([A-Za-z0-9/\- ]{0,25})(VCEP|MDEP|EP)\b|"
    r"by this VCEP|by (ClinGen )?[A-Z]+ VCEP specifications|"
    r"classified according to [A-Za-z ]+ VCEP specifications",
    re.IGNORECASE,
)
_WITHOUT_RE = re.compile(r"without (?:using )?(PS1|PM5)\b", re.IGNORECASE)
_PROTEIN_RE = re.compile(r"protein level|not due to aberrant splicing|not a predicted or confirmed splice",
                         re.IGNORECASE)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cspec", type=Path,
                    default=Path("resources/clingen/cspec_json/cspec_summary.json"))
    ap.add_argument("--out", type=Path, default=Path("resources/shared/vcep_comparator_rules.tsv"))
    a = ap.parse_args()
    data = json.loads(a.cspec.read_text(encoding="utf-8"))["data"]
    rows: dict[tuple[str, str], dict] = {}
    for spec in data:
        genes = [g["label"] for g in spec.get("genes", [])]
        for c in spec.get("codes", []):
            crit = c.get("label")
            if crit not in ("PS1", "PM5"):
                continue
            t = (c.get("text") or "").replace("\n", " ")
            if not _VCEP_RE.search(t):
                continue
            excl = sorted({m.upper() for m in _WITHOUT_RE.findall(t)})
            prot = "1" if _PROTEIN_RE.search(t) else "0"
            for g in genes:
                r = rows.setdefault((g, crit), dict(gene=g, criterion=crit, comparator="vcep",
                                                    exclude_codes=set(), protein_level="0",
                                                    source_text=t[:160]))
                r["exclude_codes"] |= set(excl)
                if prot == "1":
                    r["protein_level"] = "1"
    with a.out.open("w", encoding="utf-8") as fh:
        fh.write("gene\tcriterion\tcomparator\texclude_codes\tprotein_level\tsource_text\n")
        for (g, crit) in sorted(rows):
            r = rows[(g, crit)]
            fh.write("\t".join([g, crit, "vcep", ",".join(sorted(r["exclude_codes"])),
                                r["protein_level"], r["source_text"].replace("\t", " ")]) + "\n")
    print(f"wrote {len(rows)} gene/criterion rules -> {a.out}")


if __name__ == "__main__":
    main()
