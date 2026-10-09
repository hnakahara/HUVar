"""ClinGen TP53 VCEP functional rules (PS3 / BS3) from systematic assay data.

Data: resources/shared/tp53_functional.tsv (built by scripts/build_tp53_functional.py
from the NCI TP53 Database R21): Kato 2003 transactivation class and the "other
eligible assays" Giacomelli 2018 (LOF = etoposide Z <= -0.21), Kotler 2018 (DNA-
binding domain only) and Kawaguchi 2005 (oligomerisation; abnormal counted as LOF).

VCEP flowchart for missense variants ("majority" = more than half of the assays
with data, counting a non-functional / functional Kato result as LOF / no LOF;
validated against the TP53 VCEP calls in the eRepo):

  Kato non-functional  + LOF by the majority of assays            -> PS3 (Strong)
  Kato non-functional  + abnormal oligomerisation (Kawaguchi)      -> PS3_Supporting
  Kato partially func. + LOF by the majority of other assays      -> PS3_Moderate
  Kato partially func. + no LOF by all other assays               -> BS3_Supporting
  Kato functional      + no LOF by the majority of assays         -> BS3 (Strong)
  no Kato data         + LOF / no LOF by the majority of assays   -> PS3_Supporting / BS3_Supporting

Caveats (applied in the registry): not with a SpliceAI-based PP3; PS3 not applied
with PVS1 at full strength and downgraded to Moderate with PVS1_Strong. Small
deletions (Kotler) and Funk 2025 are not implemented.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength

SRC_TP53_FUNCTIONAL = "[TP53 VCEP functional data]"

_AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
        "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
        "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}


@dataclass(frozen=True)
class TP53Assays:
    kato: str = ""        # nonfunctional / partial / functional / ""
    giacomelli: str = ""  # LOF / noLOF / ""
    kotler: str = ""
    kawaguchi: str = ""   # LOF (= abnormal oligomerisation) / noLOF / ""

    def others(self) -> dict[str, str]:
        return {k: v for k, v in (("Giacomelli", self.giacomelli), ("Kotler", self.kotler),
                                  ("Kawaguchi", self.kawaguchi)) if v}

    def describe(self) -> str:
        parts = [f"Kato={self.kato or 'NA'}"] + [f"{k}={v}" for k, v in self.others().items()]
        return ", ".join(parts)


class TP53Functional:
    """Loader + decision for the TP53 VCEP functional flowchart. A missing file
    degrades to "no data" (the caller then falls back to its generic path)."""

    def __init__(self, tsv_path: Optional[Path]) -> None:
        self._d: dict[str, TP53Assays] = {}
        try:
            if tsv_path is None:
                return
            p = Path(tsv_path)
            if not p.exists():
                p = Path(__file__).resolve().parents[3] / "resources" / "shared" / p.name
            if not p.exists():
                return
            with p.open(encoding="utf-8") as fh:
                rows = (ln for ln in fh if not ln.startswith("#"))
                for r in csv.DictReader(rows, delimiter="\t"):
                    self._d[r["aa_change"]] = TP53Assays(
                        r.get("kato") or "", r.get("giacomelli") or "",
                        r.get("kotler") or "", r.get("kawaguchi") or "")
        except Exception:
            self._d.clear()

    def __bool__(self) -> bool:
        return bool(self._d)

    @staticmethod
    def aa_change(amino_acids: Optional[str], position: Optional[int],
                  hgvs_p: Optional[str] = None) -> Optional[str]:
        """One-letter missense change (e.g. R175H) from VEP amino_acids/position,
        or from a three-letter HGVS p. string."""
        if amino_acids and position and "/" in amino_acids:
            ref, alt = amino_acids.split("/", 1)
            if len(ref) == 1 and len(alt) == 1 and ref != alt and "*" not in (ref, alt):
                return f"{ref}{position}{alt}"
        if isinstance(hgvs_p, str) and hgvs_p:
            import re
            m = re.search(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)?$", hgvs_p.split(":")[-1])
            if m and m.group(1) in _AA3 and m.group(3) in _AA3 and m.group(1) != m.group(3):
                return f"{_AA3[m.group(1)]}{m.group(2)}{_AA3[m.group(3)]}"
        return None

    def lookup(self, aa: Optional[str]) -> Optional[TP53Assays]:
        return self._d.get(aa) if aa else None

    @staticmethod
    def decide(a: TP53Assays) -> Optional[tuple[ACMGCriterion, CriterionStrength, str]]:
        o = a.others()
        n = len(o)
        n_lof = sum(v == "LOF" for v in o.values())
        n_no = sum(v == "noLOF" for v in o.values())
        all_no = n > 0 and n_no == n
        # "Majority of available assays": the flowchart counts Kato among the
        # eligible assays (non-functional = LOF, functional = no LOF); a partially
        # functional Kato result is not counted either way.
        k_lof = 1 if a.kato == "nonfunctional" else 0
        k_no = 1 if a.kato == "functional" else 0
        tot = n + k_lof + k_no
        maj_lof = tot > 0 and 2 * (n_lof + k_lof) > tot
        maj_no = tot > 0 and 2 * (n_no + k_no) > tot
        S, PS3, BS3 = CriterionStrength, ACMGCriterion.PS3, ACMGCriterion.BS3
        if a.kato == "nonfunctional":
            if maj_lof:
                return PS3, S.STRONG, "Kato non-functional and LOF by the majority of assays"
            if a.kawaguchi == "LOF":
                return PS3, S.SUPPORTING, "Kato non-functional and abnormal oligomerisation (Kawaguchi)"
        elif a.kato == "partial":
            if maj_lof:
                return PS3, S.MODERATE, "Kato partially functional and LOF by the majority of other assays"
            if all_no:
                return BS3, S.SUPPORTING, "Kato partially functional and no LOF by all other assays"
        elif a.kato == "functional":
            if maj_no:
                return BS3, S.STRONG, "Kato functional and no LOF by the majority of assays"
        else:
            if maj_lof:
                return PS3, S.SUPPORTING, "no Kato data; LOF by the majority of available assays"
            if maj_no:
                return BS3, S.SUPPORTING, "no Kato data; no LOF by the majority of available assays"
        return None
