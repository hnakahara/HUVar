#!/usr/bin/env python3
"""Build a per-gene PM1 hotspot table from the ClinGen cspec_summary.json export.

PM1 ("located in a mutational hotspot / critical functional domain") is defined
per gene in free text by each VCEP. This script mines those PM1 descriptions for
machine-readable hotspot evidence — residue RANGES (e.g. "codons 167-931",
"amino acids 271-292") and explicit RESIDUE positions (e.g. "Arg158", "R107",
"residues 175, 245, 248") — and writes one row per (gene, strength) to
``pm1_hotspots.tsv`` (columns: gene_symbol, strength, regions, residues).

Coverage is best-effort: rules expressed only as exon numbers, external tables
("Supp. Table 4"), or cancerhotspots.org occurrence counts cannot be resolved to
residue positions and are skipped. Multi-gene PM1 entries are skipped because a
region quoted in the prose usually applies to one gene, not all of them.

Usage:
    python scripts/build_pm1_hotspots.py \
        --summary resources/clingen/cspec_json/cspec_summary.json \
        --out resources/clingen/pm1_hotspots.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import re

_AA3 = "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val"
_AA1 = "ACDEFGHIKLMNPQRSTVWY"

# Tokens that carry stray numbers (transcripts, PMIDs, citations, HTML, coding
# HGVS) — stripped before parsing so they are not mistaken for residues.
# NOTE: the bracket rule strips ONLY numeric-citation brackets ("[12]",
# "[1-3]", "[10, 11]"). VCEPs (RASopathy etc.) write the actual hotspot ranges
# inside square brackets as "[AA 10-17]" — those contain letters and MUST be
# preserved, or every bracketed-domain PM1 (KRAS/HRAS/PTPN11/RAF1/…) loses its
# regions and then trips the not-applicable fallback. So the bracket class is
# restricted to digits/whitespace/separators only.
_NOISE = re.compile(
    r"ENS[TPG]\d+(?:\.\d+)?|N[MPR]_?\d+(?:\.\d+)?|c\.-?\d+|PMID:?\s*\d+|pmid_\d+"
    r"|\[[\s\d,;.–-]+\]|<[^>]+>|&[a-z]+;",
    re.IGNORECASE,
)
# An "AA<n> - AA<m>" span (e.g. "Ser151 - Pro153") → residue range.
_AA_RANGE = re.compile(
    rf"(?:{_AA3}|[{_AA1}])\s?(\d{{1,4}})\s*[-–]\s*(?:{_AA3}|[{_AA1}])\s?(\d{{1,4}})"
)
# A bare numeric range "167-931". Occurrence-count ranges (cancerhotspots) are
# excluded by a trailing-context check below.
_RANGE = re.compile(r"(?<![\d.])(\d{1,4})\s*[-–]\s*(\d{1,4})(?!\d)")
_AA3_RES = re.compile(rf"\b(?:{_AA3})\s?(\d{{1,4}})\b")
_AA1_RES = re.compile(rf"\b([{_AA1}])(\d{{2,4}})\b")
# "AA <n>"-prefixed single residues (e.g. PTPN11 "AA 247, AA 251, AA 256").
# Ranges in this notation ("AA 10-17") are already captured by _RANGE; this
# picks up the isolated residues a VCEP lists alongside them.
_AA_PREFIX_RES = re.compile(r"\bAA\s?(\d{1,4})\b", re.IGNORECASE)
# A comma-separated bare-number list introduced by codon/residue wording.
_CODON_LIST = re.compile(
    r"(?:codons?|residues?|amino acids?)[^:.\n]*?[:\s]\s*((?:\d{1,4}\s*,\s*)+\d{1,4})",
    re.IGNORECASE,
)
# Occurrence/instance wording that follows a count range, not a residue range.
_OCCURRENCE = re.compile(r"\b(occurrence|instances?|somatic occurrence)", re.IGNORECASE)
_NOT_APPLICABLE = re.compile(
    r"does not apply|not applicable|highly polymorphic", re.IGNORECASE
)
# A positive applicability clause ("Applicable only to ... domains/residues").
# When present, a trailing "Not applicable to specific amino acid residues
# (see PM5)" is a PM5 redirect for the *rest* of the gene — NOT a gene-level
# negation — so the gene must NOT be marked not_applicable.
_POSITIVE_APPLICABLE = re.compile(r"\bapplicable\b\s+(?:only\s+)?to\b", re.IGNORECASE)

_VALID_STRENGTH = {"Supporting", "Moderate", "Strong"}


def _normalise(text: str) -> str:
    # Remove markdown escape backslashes ("NM\_00546.4", "PM1\_strong") so the
    # noise regex can match transcript/identifier tokens; drop thousands
    # separators ("2,101" -> "2101"); then strip noise tokens carrying unrelated
    # numbers (transcripts, PMIDs, citations).
    text = text.replace("\\", "")
    text = re.sub(r"(\d),(\d{3})\b", r"\1\2", text)
    return _NOISE.sub(" ", text)


def parse_regions(text: str) -> tuple[list[tuple[int, int]], list[int]]:
    t = _normalise(text)
    ranges: set[tuple[int, int]] = set()
    residues: set[int] = set()

    # b > 4 rejects enumeration / domain-numbering noise ("Exons 1-3",
    # "Cys2-Cys3", "1. ... 4.") — no real protein hotspot ends within the first
    # four residues.
    for m in _AA_RANGE.finditer(t):
        a, b = int(m.group(1)), int(m.group(2))
        if 0 < a < b <= 9999 and b > 4 and b - a < 2000:
            ranges.add((a, b))

    for m in _RANGE.finditer(t):
        a, b = int(m.group(1)), int(m.group(2))
        if not (0 < a < b <= 9999 and b > 4 and b - a < 2000):
            continue
        # Skip "2-9 occurrences" style count ranges (cancerhotspots wording).
        if _OCCURRENCE.search(t[m.end():m.end() + 30]):
            continue
        ranges.add((a, b))

    for m in _AA3_RES.finditer(t):
        residues.add(int(m.group(1)))
    for m in _AA1_RES.finditer(t):
        residues.add(int(m.group(2)))
    for m in _AA_PREFIX_RES.finditer(t):
        residues.add(int(m.group(1)))
    for m in _CODON_LIST.finditer(t):
        for n in re.findall(r"\d{1,4}", m.group(1)):
            residues.add(int(n))

    # Drop residues <5 (disulfide-bond numbering like "Cys2-Cys3", start codon)
    # and those already covered by a range, to keep the table compact and clean.
    residues = {
        r for r in residues
        if r >= 5 and not any(a <= r <= b for a, b in ranges)
    }
    return sorted(ranges), sorted(residues)


def build(summary_path: str) -> dict[tuple[str, str], tuple[set, set]]:
    with open(summary_path, encoding="utf-8") as fh:
        data = json.load(fh)["data"]

    # (gene, strength) -> (ranges set, residues set); plus not_applicable genes.
    table: dict[tuple[str, str], tuple[set, set]] = {}
    not_applicable: set[str] = set()

    for item in data:
        genes = [
            g.get("label") if isinstance(g, dict) else g
            for g in (item.get("genes") or [])
        ]
        genes = [g for g in genes if g]
        if len(genes) != 1:  # single-gene PM1 only (avoid mis-assigning regions)
            continue
        gene = genes[0]
        for code in item.get("codes", []):
            if code.get("label") != "PM1":
                continue
            strength = code.get("strengthDescriptor")
            if strength not in _VALID_STRENGTH:
                continue
            text = (code.get("text") or "").strip()
            if not text:
                continue
            ranges, residues = parse_regions(text)
            if not ranges and not residues:
                # No resolvable region. Record an explicit not-applicable only
                # when the VCEP negates the rule for the WHOLE gene — not when
                # "not applicable" is merely a PM5 redirect sub-clause sitting
                # next to a positive "Applicable only to … domains" statement.
                if _NOT_APPLICABLE.search(text) and not _POSITIVE_APPLICABLE.search(text):
                    not_applicable.add(gene)
                continue
            key = (gene, strength)
            r, res = table.setdefault(key, (set(), set()))
            r.update(ranges)
            res.update(residues)

    # Materialise not_applicable as its own strength row (only if the gene has no
    # positive hotspot rows from any spec).
    positive_genes = {g for g, _ in table}
    for gene in not_applicable - positive_genes:
        table[(gene, "not_applicable")] = (set(), set())
    return table


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="resources/clingen/cspec_json/cspec_summary.json")
    ap.add_argument("--out", default="resources/clingen/pm1_hotspots.tsv")
    args = ap.parse_args()

    table = build(args.summary)
    rows = []
    for (gene, strength), (ranges, residues) in sorted(table.items()):
        rows.append({
            "gene_symbol": gene,
            "strength": strength,
            "regions": ";".join(f"{a}-{b}" for a, b in sorted(ranges)),
            "residues": ",".join(str(r) for r in sorted(residues)),
        })
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["gene_symbol", "strength", "regions", "residues"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    n_genes = len({r["gene_symbol"] for r in rows})
    print(f"PM1 hotspot rows: {len(rows)} | genes: {n_genes} | written → {args.out}")


if __name__ == "__main__":
    main()
