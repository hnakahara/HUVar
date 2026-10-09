#!/usr/bin/env python3
"""Build resources/shared/tp53_functional.tsv for the TP53 VCEP PS3/BS3 rules.

Source: NCI TP53 Database R21 (https://tp53.cancer.gov/get_tp53data)
  MutationView_r21.csv            TransactivationClass (Kato 2003, yeast, median of
                                  8 promoters: non-functional <=20%, partially
                                  functional 20-75%, functional >75%, supertrans
                                  treated as functional) and DNE_LOFclass
                                  (Giacomelli 2018, A549 growth suppression; LOF =
                                  etoposide Z <= -0.21)
  FunctionIshiokaDownload_r21.csv Oligomerisation_yeast (Kawaguchi 2005; MON/DIM =
                                  abnormal, TETR = normal)
  FunctionDownload_r21.csv        Kotler 2018 (H1299 NGS growth assay, DNA-binding
                                  domain): "GS (in vitro)" in Loss_of_Function =
                                  LOF, in Conserved_WT_Function = no LOF
Funk 2025 (CRISPR, a few exons) is not in the database and is not used.

Output: one row per missense protein change (one-letter, e.g. R175H):
  aa_change  kato  giacomelli  kotler  kawaguchi
with kato in {nonfunctional, partial, functional}, the others in {LOF, noLOF}
(kawaguchi: LOF = abnormal oligomerisation), empty = no data.

  python scripts/build_tp53_functional.py --dir tp53_functional \
      --out resources/shared/tp53_functional.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd

KATO = {"non-functional": "nonfunctional", "partially functional": "partial",
        "functional": "functional", "supertrans": "functional"}
GIAC = {"DNE_LOF": "LOF", "notDNE_LOF": "LOF", "notDNE_notLOF": "noLOF"}
KAWA = {"MON": "LOF", "DIM": "LOF", "TETR": "noLOF"}
_MISSENSE = re.compile(r"^p\.([A-Z])(\d+)([A-Z])$")


def _one(series: pd.Series, name: str) -> dict[str, str]:
    """aa_change -> value; drop protein changes with conflicting values."""
    out: dict[str, str] = {}
    bad: set[str] = set()
    for aa, v in series.dropna().items():
        if aa in out and out[aa] != v:
            bad.add(aa)
        out[aa] = v
    for aa in bad:
        out.pop(aa)
    if bad:
        print(f"{name}: {len(bad)} protein changes with conflicting classes dropped")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("resources/shared/tp53_functional.tsv"))
    a = ap.parse_args()

    mv = pd.read_csv(a.dir / "MutationView_r21.csv", dtype=str, low_memory=False)
    mv = mv[mv["ProtDescription"].fillna("").str.match(_MISSENSE)]
    mv = mv.assign(aa=mv["ProtDescription"].str[2:]).set_index("aa")
    kato = _one(mv["TransactivationClass"].map(KATO), "kato")
    giac = _one(mv["DNE_LOFclass"].map(GIAC), "giacomelli")

    ish = pd.read_csv(a.dir / "FunctionIshiokaDownload_r21.csv", dtype=str)
    ish = ish[ish["ProtDescription"].fillna("").str.match(_MISSENSE)]
    ish = ish.assign(aa=ish["ProtDescription"].str[2:]).set_index("aa")
    kawa = _one(ish["Oligomerisation_yeast"].str.strip().map(KAWA), "kawaguchi")

    fn = pd.read_csv(a.dir / "FunctionDownload_r21.csv", dtype=str, encoding_errors="replace")
    ko = fn[fn["Authors"].str.contains("Kotler", na=False)
            & fn["Method"].eq("NGS-based growth assay")
            & fn["ProtDescription"].fillna("").str.match(_MISSENSE)]
    lof = ko["Loss_of_Function"].fillna("").str.contains("GS (in vitro)", regex=False)
    keep = ko["Conserved_WT_Function"].fillna("").str.contains("GS (in vitro)", regex=False)
    kv = pd.Series(None, index=ko.index, dtype=object)
    kv[lof & ~keep] = "LOF"
    kv[keep & ~lof] = "noLOF"
    kotler = _one(pd.Series(kv.values, index=ko["ProtDescription"].str[2:]), "kotler")

    keys = sorted(set(kato) | set(giac) | set(kawa) | set(kotler),
                  key=lambda s: (int(re.sub(r"\D", "", s)), s))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="", encoding="utf-8") as fh:
        fh.write("# TP53 functional assay classes for the ClinGen TP53 VCEP PS3/BS3 rules. "
                 "Derived from The TP53 Database (R21, Jan 2025), https://tp53.cancer.gov "
                 "(de Andrade et al., Cell Death Differ 2022): Kato 2003, Giacomelli 2018, "
                 "Kotler 2018, Kawaguchi 2005. Free use with acknowledgement; use resulting "
                 "in monetization is prohibited (https://tp53.cancer.gov/about). "
                 "Built by scripts/build_tp53_functional.py\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["aa_change", "kato", "giacomelli", "kotler", "kawaguchi"])
        for k in keys:
            w.writerow([k, kato.get(k, ""), giac.get(k, ""), kotler.get(k, ""), kawa.get(k, "")])
    print(f"{len(keys)} protein changes: kato {len(kato)}, giacomelli {len(giac)}, "
          f"kotler {len(kotler)}, kawaguchi {len(kawa)} -> {a.out}")


if __name__ == "__main__":
    main()
