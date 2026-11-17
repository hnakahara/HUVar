#!/usr/bin/env python3
"""Build resources/shared/vcep_functional_codes.tsv: variant-level PS3/BS3 codes
taken from functional-evidence files that ClinGen VCEPs attach to their CSpec
specifications (https://cspec.genome.network/cspec/ui/svi/doc/GNxxx).

Sources (download the files from each specification's document page):

* TP53 (GN009) ``Functional-worksheet.xlsx`` — Supplementary Table S3: the VCEP's
  preliminary PS3 / BS3 code for every possible p53 missense variant (TP53
  specifications v2: Kato, Funk 2025, Giacomelli, Kotler and Kawaguchi assays).
  When the code depends on the nucleotide change the row carries it in
  ``seq_change`` (a codon change such as ``CGG>CTA``, a base change within the
  codon such as ``G>T``, or a c. change).
* OTC (GN156) ``OTC-Functional-Supplement.xlsx`` — Table S7 of Lo et al. 2023
  (PMID 37146589): yeast complementation growth of 1,570 OTC missense variants,
  wild type = 1. Coded with the VCEP thresholds (spec text and
  ``OTC_Functional-assay-SVI-documentation``): growth <= 50% PS3_Moderate,
  > 50% and <= 65% PS3_Supporting, > 65% and < 80% BS3_Supporting, >= 80%
  BS3_Moderate. Not applied to the SMG motif (aa 264-276) or to variants the study
  flags for a possible splicing impact.
* GALT (GN158) ``PS3_Functional-assay_GALT_final_for-Clin-Gen_updated-July-2024.xlsx``
  — sheet "All variants (final)": residual activity (% of wild type) and the
  VCEP code (PS3_Supporting <= 20%, BS3_Supporting > 50%).
* PAH (GN006) ``PS3_Functional-assay_PAH.xlsx`` — sheet "PAH variants (Table S8)",
  the VCEP pilot assignments of Kanavy et al. 2019 (Genome Med 11:77). Variants the
  VCEP has since classified (``--erepo``: the eRepo snapshot
  resources/GRCh38/erepo_classifications_hg38.tsv) take the PS3/BS3 code and
  strength of that current classification, including "no PS3/BS3". For the other
  variants a pilot PS3 (Strong) is stored as PS3_Moderate: the current
  specification makes PS3_Strong not applicable.
* IDUA (GN091) ``Appendix-3-functional-assays_IDUA.xlsx`` — sheet "Variants":
  residual activity (% of wild type; numeric cells are percent-formatted) per
  study, coded with the study thresholds of the "IDUA functional assays" sheet:
  PS3_Supporting below 1% (Yu 2020) or 2% (Beesley 2001, Yogalingam 2004,
  Matte 2003); PS3_Moderate for the Matte 2003 variants with Western blot
  evidence (p.Arg48del, p.Gly208Asp, p.Ser260Phe, p.Asn350Ile, p.Pro533Arg,
  p.Arg628Pro); BS3_Supporting above 10% (Yu: "~10%", the range that holds the
  four pseudodeficiency variants). A variant studied in several papers is coded
  only when all of them agree; values between the thresholds give no code.

Usage:
    python scripts/build_vcep_functional_codes.py \\
        --tp53 Functional-worksheet.xlsx --otc OTC-Functional-Supplement.xlsx \\
        --galt PS3_Functional-assay_GALT_final_for-Clin-Gen_updated-July-2024.xlsx \\
        --pah PS3_Functional-assay_PAH.xlsx \\
        --erepo resources/GRCh38/erepo_classifications_hg38.tsv \\
        --idua Appendix-3-functional-assays_IDUA.xlsx \\
        --out resources/shared/vcep_functional_codes.tsv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

AA = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
      "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
      "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "X": "Ter",
      "*": "Ter"}
COLUMNS = ("gene", "transcript", "aa_change", "seq_change", "hgvs_c", "code", "strength",
           "value", "source", "note")


def _sheet_rows(path: Path, sheet_match: str):
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = next(ws for ws in wb.worksheets if sheet_match.lower() in ws.title.lower())
    for row in ws.iter_rows(values_only=True):
        yield ["" if v is None else str(v).strip() for v in row]


def _one_to_three(change: str) -> str | None:
    """'R175H' -> 'Arg175His' (also 'R204X' -> 'Arg204Ter'), else None."""
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z*])", change.strip())
    if not m or m.group(1) not in AA or m.group(3) not in AA:
        return None
    return f"{AA[m.group(1)]}{m.group(2)}{AA[m.group(3)]}"


def _code(text: str) -> tuple[str, str]:
    """'PS3_Supporting' / 'PS3 supporting' / 'BS3' / 'No evidence' -> (code, strength)."""
    t = text.strip().lower().replace(" ", "_")
    if not t or t in ("no_evidence", "not_met", "no", "none"):
        return "none", ""
    m = re.fullmatch(r"(ps3|bs3)(?:_(supporting|moderate|strong|very_strong))?", t)
    if not m:
        raise ValueError(f"unrecognised functional code {text!r}")
    return m.group(1).upper(), (m.group(2) or "strong")


def tp53(path: Path) -> list[dict]:
    out, header = [], None
    for r in _sheet_rows(path, "S3"):
        if r and r[0] == "Protein change":
            header = r
            continue
        if header is None or not r or not r[0] or not _one_to_three(r[0]):
            continue
        code, strength = _code(r[7])
        seq = r[1]
        note = (f"Kato={r[2]}; Funk={r[3]}; Giacomelli={r[4]}; Kotler={r[5]}; "
                f"Kawaguchi={r[6]}")
        if seq.startswith("AASub"):
            note += "; multi-nucleotide substitution (no Kato data)"
            seq = ""
        elif seq == "NA":
            seq = ""
        out.append(dict(gene="TP53", transcript="NM_000546", aa_change=_one_to_three(r[0]),
                        seq_change=seq, hgvs_c="", code=code, strength=strength, value="",
                        source="TP53 VCEP functional worksheet (Table S3)", note=note))
    return out


def otc(path: Path) -> list[dict]:
    out, started = [], False
    for r in _sheet_rows(path, "Table S7"):
        if r and r[0] == "Human_AA_num":
            started = True
            continue
        if not started or not r or not r[1]:
            continue
        aa = r[1]
        pos = int(r[0])
        growth = float(r[3])
        splice = r[5] == "True"
        if 264 <= pos <= 276:
            code, strength, note = "none", "", "SMG motif (aa 264-276): yeast assay not applied"
        elif splice:
            code, strength, note = "none", "", "possible splicing impact (Lo 2023): not applied"
        elif growth <= 0.50:
            code, strength, note = "PS3", "moderate", "growth <= 50% of wild type"
        elif growth <= 0.65:
            code, strength, note = "PS3", "supporting", "growth > 50% and <= 65% of wild type"
        elif growth < 0.80:
            code, strength, note = "BS3", "supporting", "growth > 65% and < 80% of wild type"
        else:
            code, strength, note = "BS3", "moderate", "growth >= 80% of wild type"
        out.append(dict(gene="OTC", transcript="NM_000531", aa_change=aa, seq_change="",
                        hgvs_c="", code=code, strength=strength, value=f"{growth:.3f}",
                        source="Lo et al. 2023 yeast complementation (PMID 37146589)",
                        note=f"{note}; class {r[7]}"))
    return out


def galt(path: Path) -> list[dict]:
    by_aa: dict[str, list[list[str]]] = {}
    for r in _sheet_rows(path, "All variants"):
        if not r or not r[0] or r[0] == "Variant":
            continue
        aa = _one_to_three(r[0])
        if aa is not None:
            by_aa.setdefault(aa, []).append(r)
    out = []
    for aa, rs in by_aa.items():
        codes = {_code(r[6]) for r in rs}
        desc = "; ".join(
            f"activity {r[1]}% of wild type ({', '.join(x for x in (r[2], r[3], r[4]) if x)})"
            + (f", group {r[5]}" if r[5] else "") + f" -> {r[6] or 'no code'}" for r in rs)
        if len(codes) == 1:
            code, strength = codes.pop()
        else:  # the table lists discordant results from different assays
            code, strength = "none", ""
            desc = "discordant results, no code: " + desc
        out.append(dict(gene="GALT", transcript="NM_000155", aa_change=aa, seq_change="",
                        hgvs_c="", code=code, strength=strength,
                        value="/".join(r[1] for r in rs),
                        source="GALT VCEP functional assay table (July 2024)", note=desc))
    return out


def _erepo_functional(path: Path | None, gene: str) -> dict[str, tuple[str, str, str]]:
    """{c.: (code, strength, label)} of the gene's eRepo VCEP classifications."""
    import csv
    out: dict[str, tuple[str, str, str]] = {}
    if path is None:
        return out
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            if r["gene"] != gene:
                continue
            codes = [c for c in r["met_codes"].split(",") if c.startswith(("PS3", "BS3"))]
            code, strength = _code(codes[0]) if codes else ("none", "")
            out[r["hgvs_c"]] = (code, strength,
                                f"{r['assertion']}, approved {r['approval_date']}")
    return out


def pah(path: Path, erepo: Path | None = None) -> list[dict]:
    current = _erepo_functional(erepo, "PAH")
    out, started = [], False
    for r in _sheet_rows(path, "S8)"):  # sheet "PAH variants (Tablle S8)"
        if r and len(r) > 1 and r[1] == "cDNA":
            started = True
            continue
        if not started or not r or not r[1].startswith("c.") or not r[3].strip():
            continue
        c = r[1].strip()
        code, strength = _code(r[3])
        pilot = "none" if code == "none" else f"{code}_{strength}"
        note = re.sub(r"\s+", " ", r[5])[:300]
        if c in current:
            code, strength, label = current[c]
            source = f"PAH VCEP classification (eRepo: {label})"
            note = f"Table S8 pilot code {pilot}; {note}"
        else:
            source = "PAH VCEP pilot (Kanavy et al. 2019, Table S8)"
            if code == "PS3" and strength == "strong":
                strength = "moderate"
                note = ("pilot PS3 stored as PS3_Moderate (current specification: "
                        f"PS3_Strong not applicable); {note}")
        out.append(dict(gene="PAH", transcript="NM_000277", aa_change="", seq_change="",
                        hgvs_c=c, code=code, strength=strength, value=r[6],
                        source=source, note=note))
    return out


# IDUA study thresholds (Appendix 3, "IDUA functional assays"): abnormal below, normal above.
_IDUA_STUDY = {"Yu": (1.0, 10.0, True), "Beesley": (2.0, 10.0, False),
               "Yogalingam": (2.0, 10.0, False), "Matte": (2.0, 10.0, False)}
_IDUA_WESTERN = {"Arg48del", "Gly208Asp", "Ser260Phe", "Asn350Ile", "Pro533Arg", "Arg628Pro"}


def _pct_range(cell) -> tuple[float, float, bool] | None:
    """Residual activity as (low, high, high_is_exclusive) in % of wild type."""
    v, fmt = cell
    if isinstance(v, (int, float)):
        x = v * 100 if "%" in fmt else float(v)
        return x, x, False
    t = str(v).lower()
    if "not detected" in t and not re.search(r"\d", t) or "no activity" in t:
        return 0.0, 0.0, False
    m = re.search(r"<\s*([\d.]+)\s*%", t)
    if m:
        return 0.0, float(m.group(1)), True
    nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", t.split(" of normal")[0])]
    if "not detected" in t:
        nums = [0.0] + nums
    if not nums:
        return None
    return min(nums[:2]), max(nums[:2]), False


def _idua_call(study: str, rng: tuple[float, float, bool]) -> str:
    abn, normal, approx = _IDUA_STUDY[study]
    lo, hi, excl = rng
    if hi < abn or (excl and hi <= abn):
        return "abnormal"
    mid = (lo + hi) / 2
    if mid > normal or (approx and mid >= normal):
        return "normal"
    return "indeterminate"


def idua(path: Path) -> list[dict]:
    import openpyxl
    ws = openpyxl.load_workbook(path, data_only=True)["Variants"]
    out = []
    for row in ws.iter_rows(min_row=2):
        c_raw, p_raw, refs_raw = (str(x.value or "").strip() for x in row[:3])
        if not c_raw.startswith("c."):
            continue
        c = c_raw.split()[0]
        aa = p_raw.replace("p.", "").replace("*", "Ter").strip()
        cell = row[3]
        refs = [x.strip() for x in refs_raw.split(";") if x.strip()]
        vals = ([(cell.value, cell.number_format)] if not isinstance(cell.value, str)
                else [(x.strip(), "") for x in cell.value.split(";") if x.strip()])
        calls, desc = [], []
        for i, ref in enumerate(refs):
            study = next((k for k in _IDUA_STUDY if ref.startswith(k)), None)
            val = vals[i] if len(vals) == len(refs) else vals[0]
            rng = _pct_range(val)
            if study is None or rng is None:
                continue  # study outside the approved list (e.g. Brooks) / no value
            call = _idua_call(study, rng)
            calls.append((study, call))
            shown = (f"{rng[1]:.3g}%" if rng[0] == rng[1] else
                     (f"<{rng[1]:.3g}%" if rng[2] else f"{rng[0]:.3g}-{rng[1]:.3g}%"))
            desc.append(f"{study} {shown} -> {call}")
        kinds = {k for _, k in calls}
        if kinds == {"abnormal"}:
            moderate = aa in _IDUA_WESTERN and any(s == "Matte" for s, _ in calls)
            code, strength = "PS3", ("moderate" if moderate else "supporting")
            if moderate:
                desc.append("Matte 2003 Western blot: abnormal synthesis / processing")
        elif kinds == {"normal"}:
            code, strength = "BS3", "supporting"
        else:
            code, strength = "none", ""
            if len(kinds) > 1:
                desc.append("studies disagree")
        out.append(dict(gene="IDUA", transcript="NM_000203", aa_change=aa, seq_change="",
                        hgvs_c=c, code=code, strength=strength,
                        value="/".join(d.split(" ")[1] for d in desc if " -> " in d),
                        source="IDUA VCEP Appendix 3 (functional assays)",
                        note="; ".join(desc)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tp53", type=Path)
    ap.add_argument("--otc", type=Path)
    ap.add_argument("--galt", type=Path)
    ap.add_argument("--pah", type=Path)
    ap.add_argument("--erepo", type=Path, help="eRepo classifications TSV (PAH current codes)")
    ap.add_argument("--idua", type=Path)
    ap.add_argument("--out", type=Path, default=Path("resources/shared/vcep_functional_codes.tsv"))
    a = ap.parse_args()
    rows: list[dict] = []
    for fn, path in ((tp53, a.tp53), (otc, a.otc), (galt, a.galt), (idua, a.idua)):
        if path:
            rows += fn(path)
    if a.pah:
        rows += pah(a.pah, a.erepo)
    with a.out.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("# Variant-level PS3/BS3 codes from functional-evidence files attached to ClinGen "
                 "VCEP specifications (CSpec). Built by scripts/build_vcep_functional_codes.py.\n")
        fh.write("# Matching: hgvs_c when given, otherwise aa_change on the MANE transcript; "
                 "seq_change (codon / base / c. change) separates codon-specific TP53 codes.\n")
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in COLUMNS) + "\n")
    print(f"{len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
