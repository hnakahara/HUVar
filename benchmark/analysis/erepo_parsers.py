#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
erepo_parsers.py

Shared parsers / normalisation for the eRepo benchmark.
Sources:
  - original           : erepo-original.csv          (reference / ground truth)
  - intervar           : erepo_{build}_{build}_multianno_txt.intervar
  - bias               : erepo_{build}_bias-2015.tsv
  - bias_nopp5bp6      : BIAS-2015 with PP5/BP6 disabled (hg19=..._ignore_PP5_BP6.tsv,
                         hg38=..._PP5_BP6.tsv; note that the file names differ between builds)
  - huhvar_default     : erepo_{build}_huhvar_default_results.tsv
  - huhvar_revel       : erepo_{build}_huhvar_revel_results.tsv
  - huhvar_spliceai    : erepo_{build}_huhvar_spliceai_results.tsv

HUVar provides two classifications (classification_2015 / classification_bayesian).

Each record is returned as a dict with the keys:
  source, build, system, gene, transcript (without version), hgvsc_raw, hgvsc_norm,
  hgvsp, coord (chr:pos:ref:alt; InterVar indels in ANNOVAR style), cls_raw, cls
"""

import csv
import re

csv.field_size_limit(10 ** 7)

# ---- Canonical labels of the 5 categories ----
#   P=Pathogenic, LP=Likely Pathogenic, VUS=Uncertain, LB=Likely Benign, B=Benign
CANON_ORDER = ["P", "LP", "VUS", "LB", "B"]

_LABEL_MAP = {
    "pathogenic": "P",
    "likely pathogenic": "LP",
    "likelypathogenic": "LP",
    "likely_pathogenic": "LP",
    "uncertain significance": "VUS",
    "uncertain_significance": "VUS",
    "uncertain": "VUS",
    "vus": "VUS",
    "likely benign": "LB",
    "likelybenign": "LB",
    "likely_benign": "LB",
    "benign": "B",
}


def canon_label(s):
    """Normalise a raw label to one of the 5 categories (P/LP/VUS/LB/B); None if unknown."""
    if s is None:
        return None
    return _LABEL_MAP.get(s.strip().lower())


def strip_ver(nm):
    """Remove the version of a transcript accession (NM_005026.3 -> NM_005026)."""
    if not nm:
        return nm
    return re.sub(r"\.\d+$", "", nm.strip())


def norm_c(c):
    """Normalise HGVSc for comparison.
    - remove a leading 'transcript:'
    - remove whitespace
    - remove trailing bases of del / dup (c.47_48delCT -> c.47_48del, c.40dupA -> c.40dup)
  (inserted bases of ins / delins are kept because they are needed for uniqueness)
    """
    if not c:
        return c
    c = c.strip()
    if ":" in c:
        c = c.split(":", 1)[1]
    c = c.replace(" ", "")
    c = re.sub(r"(del)[ACGTN]+$", r"\1", c, flags=re.IGNORECASE)
    c = re.sub(r"(dup)[ACGTN]+$", r"\1", c, flags=re.IGNORECASE)
    return c


# ---- ACMG codes ----
import json

ACMG_CODES = [
    "PVS1",
    "PS1", "PS2", "PS3", "PS4",
    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
    "PP1", "PP2", "PP3", "PP4", "PP5",
    "BA1",
    "BS1", "BS2", "BS3", "BS4",
    "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7",
]


def codes_from_erepo(s):
    """Extract the base codes from the eRepo 'Applied Evidence Codes (Met)'.
    e.g. 'PM6,PM2,PS4_Supporting' -> [PM6,PM2,PS4]"""
    out = []
    for tok in (s or "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        m = re.match(r"([A-Za-z]+\d+)", tok)
        if m:
            out.append(m.group(1).upper())
    return out


def codes_from_bias_rationale(rat):
    """Extract the codes with score > 0 from the BIAS rationale (JSON)."""
    try:
        d = json.loads(rat)
    except Exception:
        return []
    out = []
    for _grp, members in d.items():
        if isinstance(members, dict):
            for code, val in members.items():
                if isinstance(val, list) and val and \
                        isinstance(val[0], (int, float)) and val[0] > 0:
                    out.append(code.upper())
    return out


# InterVar arrays are 0-based (PS[0]=PS1 ...). Trailing padding is ignored by using only maxn codes.
_IV_ARRAYS = [("PS", 4), ("PM", 6), ("PP", 5), ("BS", 4), ("BP", 7)]


def codes_from_intervar(val):
    out = []
    if re.search(r"PVS1=([1-9])", val or ""):
        out.append("PVS1")
    if re.search(r"BA1=([1-9])", val or ""):
        out.append("BA1")
    for prefix, maxn in _IV_ARRAYS:
        m = re.search(prefix + r"=\[([0-9,\s]+)\]", val or "")
        if m:
            arr = [int(x) for x in m.group(1).split(",")]
            for i, v in enumerate(arr):
                if i < maxn and v > 0:
                    out.append(f"{prefix}{i + 1}")
    return out


def codes_from_huhvar(row):
    out = []
    for c in ACMG_CODES:
        v = str(row.get(c, "0")).strip()
        if v and v not in ("0", "NotMet", "False", "."):
            out.append(c)
    return out


def _rec(source, build, system, gene, transcript, hgvsc_raw, hgvsp, coord, cls_raw):
    return {
        "source": source,
        "build": build,
        "system": system,
        "gene": (gene or "").strip(),
        "transcript": strip_ver(transcript),
        "hgvsc_raw": (hgvsc_raw or "").strip(),
        "hgvsc_norm": norm_c(hgvsc_raw),
        "hgvsp": (hgvsp or "").strip(),
        "coord": coord,
        "cls_raw": (cls_raw or "").strip(),
        "cls": canon_label(cls_raw),
    }


# ----------------------------------------------------------------------
# Reference (original)
# ----------------------------------------------------------------------
def load_original(path):
    """Read the reference and return two lookup dictionaries and the list of records.

    Returns: dict(by_nm=..., by_gene=..., records=...)
      by_nm  : (transcript_noversion, hgvsc_norm) -> cls(canonical)
      by_gene: (gene, hgvsc_norm)                 -> cls(canonical)
    """
    by_nm = {}
    by_gene = {}
    records = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            gene = row["Gene"].strip()
            nm = strip_ver(row["RefSeq (NM_)"])
            craw = row["HGVSc"].strip()
            cn = norm_c(craw)
            cls = canon_label(row["Assertion"])
            rec = _rec("original", None, None, gene, nm, craw,
                       row.get("HGVSp", ""), None, row["Assertion"])
            rec["evidence_codes"] = row.get("Applied Evidence Codes (Met)", "")
            rec["codes"] = codes_from_erepo(rec["evidence_codes"])
            records.append(rec)
            if not cn or cls is None:
                continue
            # several assertions may share a key; the first one is kept (rare)
            by_nm.setdefault((nm, cn), cls)
            by_gene.setdefault((gene, cn), cls)
    return {"by_nm": by_nm, "by_gene": by_gene, "records": records}


# ----------------------------------------------------------------------
# BIAS-2015
# ----------------------------------------------------------------------
def load_bias(path, build, system="bias"):
    """Read a BIAS-2015 TSV.
    system: system tag of the output records; "bias" for the published version and e.g.
        "bias_nopp5bp6" for the version with PP5/BP6 disabled (same columns, same processing).
    """
    out = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            hc = row["hgvsc"].strip()  # e.g. "NM_005026.3:c.4C>A"
            if ":" in hc:
                tr, craw = hc.split(":", 1)
            else:
                tr, craw = row.get("transcript", ""), hc
            coord = f"{row['chromosome']}:{row['position']}:{row['refAllele']}:{row['altAllele']}"
            rec = _rec("bias", build, system, row["geneName"], tr,
                       craw, row.get("hgvsp", ""), coord,
                       row["acmgClassification"])
            rec["codes"] = codes_from_bias_rationale(row.get("rationale", ""))
            out.append(rec)
    return out


# ----------------------------------------------------------------------
# HUHVar (default / spliceai), system = 2015 | bayesian
# ----------------------------------------------------------------------
def load_huhvar(path, build, source, system):
    assert system in ("2015", "bayesian")
    col = "classification_2015" if system == "2015" else "classification_bayesian"
    out = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            coord = row["variant_id"].strip()  # already chr:pos:ref:alt (VCF style)
            rec = _rec(source, build, system, row["gene"], row["transcript"],
                       row["hgvs_c"], row.get("hgvs_p", ""), coord, row[col])
            rec["codes"] = codes_from_huhvar(row)
            out.append(rec)
    return out


# ----------------------------------------------------------------------
# InterVar (ANNOVAR multianno + InterVar)
# ----------------------------------------------------------------------
_IV_RE = re.compile(r"InterVar:\s*(.*?)\s+PVS1=")


def _intervar_class(val):
    m = _IV_RE.search(val or "")
    return m.group(1).strip() if m else (val or "").strip()


def load_intervar(path, build):
    out = []
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        iv_col = [c for c in r.fieldnames if "InterVar" in c][0]
        oi_col = r.fieldnames[-1]  # Otherinfo: "GENE:HGVSc"
        for row in r:
            oi = row[oi_col].strip()
            if ":" in oi:
                g_oi, craw = oi.split(":", 1)
            else:
                g_oi, craw = "", oi
            gene = (row.get("Ref.Gene") or g_oi).strip()
            cls_raw = _intervar_class(row[iv_col])
            # InterVar coordinates are in ANNOVAR style; they are resolved to VCF coordinates with the truth annovar2coord map.
            rec = _rec("intervar", build, "intervar", gene, "",
                       craw, "", None, cls_raw)
            rec["annovar"] = (row["#Chr"], int(row["Start"]), int(row["End"]),
                              row["Ref"], row["Alt"])
            rec["codes"] = codes_from_intervar(row[iv_col])
            out.append(rec)
    return out


# ----------------------------------------------------------------------
# Matching helpers
# ----------------------------------------------------------------------
def match_assertion(rec, orig):
    """Match one record to the reference and return (assertion, key_used); (None, reason) if not matched.
    key_used: 'nm+c' | 'gene+c' | None
    reason  : 'no_hgvsc' | 'transcript_unmatched' | 'hgvsc_unmatched'
    """
    cn = rec["hgvsc_norm"]
    if not cn:
        return None, "no_hgvsc"
    nm = rec["transcript"]
    gene = rec["gene"]
    if nm and (nm, cn) in orig["by_nm"]:
        return orig["by_nm"][(nm, cn)], "nm+c"
    if gene and (gene, cn) in orig["by_gene"]:
        return orig["by_gene"][(gene, cn)], "gene+c"
    # infer the reason for a failed match
    nm_in_orig = nm and any(k[0] == nm for k in orig["by_nm"])
    gene_in_orig = gene and any(k[0] == gene for k in orig["by_gene"])
    if not (nm_in_orig or gene_in_orig):
        return None, "transcript_unmatched"
    return None, "hgvsc_unmatched"


# ----------------------------------------------------------------------
# Truth coordinate map based on the VCF (TransVar conversion)
# ----------------------------------------------------------------------
def load_vcf(path):
    """Read the TransVar-derived VCF and return (chrom, pos, ref, alt, name) for each row.
    name is INFO NAME= (or the ID column) = 'GENE:HGVSc'.
    """
    out = []
    with open(path, newline="", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 8:
                continue
            m = re.search(r"NAME=([^;]+)", p[7])
            name = m.group(1) if m else p[2]
            out.append((p[0], int(p[1]), p[3], p[4], name))
    return out


def vcf_to_annovar(chrom, pos, ref, alt):
    """Convert a VCF record (1-based, with anchor base) deterministically to ANNOVAR (chr,start,end,ref,alt).
    No reference genome is needed (derived from the REF/ALT bases).
    Used to match InterVar (ANNOVAR multianno) coordinates.
    """
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if len(ref) == 1 and len(alt) == 1:                 # SNV
        return (c, pos, pos, ref, alt)
    if ref and alt and ref[0] == alt[0]:                # anchored (shared first base)
        r, a = ref[1:], alt[1:]
        if r == "" and a != "":                         # insertion
            return (c, pos, pos, "-", a)
        if a == "" and r != "":                         # deletion
            return (c, pos + 1, pos + len(ref) - 1, r, "-")
        if r != "" and a != "":                         # delins
            return (c, pos + 1, pos + len(ref) - 1, r, a)
    return (c, pos, pos + len(ref) - 1, ref, alt)       # block substitution


def build_truth(vcf_path, original):
    """Build the build-specific truth map from the VCF and the original.

    Returns a dict:
  coord2assert : {"chrom:pos:ref:alt" -> cls(canonical)}   (VCF-style coordinates)
  annovar2coord: {(c,start,end,ref,alt) -> "chrom:pos:ref:alt"}  (for matching InterVar)
  conflicts    : set of coordinate strings excluded because of conflicting assertions
      stats        : {n_vcf, n_in_original, n_unique_coord, n_conflict}
    """
    name2assert = {}
    # GENE:HGVSc of the original -> cls
    for r in original["records"]:
        if not r["hgvsc_raw"] or r["cls"] is None:
            continue
        name2assert.setdefault(f"{r['gene']}:{r['hgvsc_raw']}", r["cls"])

    vcf = load_vcf(vcf_path)
    by_coord = {}
    annovar2coord = {}
    n_in_orig = 0
    for chrom, pos, ref, alt, name in vcf:
        coord = f"{chrom}:{pos}:{ref}:{alt}"
        annovar2coord[vcf_to_annovar(chrom, pos, ref, alt)] = coord
        a = name2assert.get(name)
        if a is None:
            continue
        n_in_orig += 1
        by_coord.setdefault(coord, set()).add(a)

    coord2assert = {}
    conflicts = set()
    for coord, asserts in by_coord.items():
        if len(asserts) == 1:
            coord2assert[coord] = next(iter(asserts))
        else:
            conflicts.add(coord)            # several truths for the same coordinate -> excluded from evaluation
    # conflicting coordinates may still be reached through annovar2coord but are absent from coord2assert, so they are excluded
    return {
        "coord2assert": coord2assert,
        "annovar2coord": annovar2coord,
        "conflicts": conflicts,
        "stats": {
            "n_vcf": len(vcf),
            "n_in_original": n_in_orig,
            "n_unique_coord": len(by_coord),
            "n_conflict": len(conflicts),
        },
    }


def resolve_truth(rec, truth):
    """Resolve one record to the VCF truth coordinate and return (assertion, coord); (None, None) if unresolved.
    - intervar: ANNOVAR 5-tuple -> VCF coordinate -> assertion
    - others: rec['coord'] (chr:pos:ref:alt) is matched directly
    """
    if rec["source"] == "intervar":
        coord = truth["annovar2coord"].get(rec.get("annovar"))
    else:
        coord = rec.get("coord")
    if coord is None:
        return None, None
    return truth["coord2assert"].get(coord), coord


# ----------------------------------------------------------------------
# File name registry
# ----------------------------------------------------------------------
def file_registry(base):
    base = base.rstrip("/")
    reg = {"original": f"{base}/erepo-original.csv", "builds": {}}
    # Note: the file names of BIAS-2015 with PP5/BP6 disabled differ between builds.
    #   hg19: erepo_hg19_bias-2015_ignore_PP5_BP6.tsv
    #   hg38: erepo_hg38_bias-2015_PP5_BP6.tsv
    # Both are treated as the same method, "BIAS-2015 ignoring PP5/BP6".
    # If the actual file names differ, only this dictionary needs to be changed.
    bias_nopp5bp6_name = {
        "hg19": "erepo_hg19_bias-2015_ignore_PP5_BP6.tsv",
        "hg38": "erepo_hg38_bias-2015_PP5_BP6.tsv",
    }
    for b in ("hg19", "hg38"):
        reg["builds"][b] = {
            "vcf": f"{base}/erepo_{b}.vcf",
            "intervar": f"{base}/erepo_{b}_{b}_multianno_txt.intervar",
            "bias": f"{base}/erepo_{b}_bias-2015.tsv",
            "bias_nopp5bp6": f"{base}/{bias_nopp5bp6_name[b]}",
            "huhvar_default": f"{base}/erepo_{b}_huhvar_default_results.tsv",
            "huhvar_revel": f"{base}/erepo_{b}_huhvar_revel_results.tsv",
            "huhvar_spliceai": f"{base}/erepo_{b}_huhvar_spliceai_results.tsv",
        }
    return reg
