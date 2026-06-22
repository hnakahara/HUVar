#!/usr/bin/env python3
"""Build a coding-exon coordinate table for the VCEP PVS1 genes from a MANE
RefSeq GFF3.

The exon-aware PVS1 splice refinement (``vcep_pvs1_exons.py``) needs, for each
clinical transcript, the *coding* length of every exon so it can decide whether
skipping that exon is in-frame (length divisible by 3) and how much of the
protein it removes. This script extracts those numbers straight from the
authoritative MANE annotation so the exon numbering is guaranteed to match the
RefSeq transcript VEP annotates against (genomic 5'→3' order, i.e. VEP's
``exon = n/N`` and ``intron = n/N`` fields).

Why this matters: a hand-entered exon→strength table is dangerous because many
PVS1 genes have a NON-CODING exon 1 (HNF4A, CDKL5, MECP2, DICER1, …); VEP counts
it as exon 1 while a VCEP decision tree may count the first CODING exon as
"exon 1", giving an off-by-one error that silently mis-scores splice variants.
Generating the table from the same coordinate space VEP uses removes that risk.

Output TSV columns (one row per CODING exon of each target transcript):
    gene, transcript, n_exons, exon, cds_len_nt, cum_cds_before_nt,
    inframe_skip, pct_protein, is_last_coding, is_penultimate_coding

  * ``exon``               1-based exon index in transcript (VEP order)
  * ``cds_len_nt``         coding bases contributed by this exon
  * ``cum_cds_before_nt``  coding bases in all earlier exons (for codon math)
  * ``inframe_skip``       "1" if skipping the exon preserves the reading frame
  * ``pct_protein``        cds_len_nt / total_cds, rounded to 4 dp
  * ``is_last_coding`` / ``is_penultimate_coding``  NMD-escape geometry helpers

Usage:
    python scripts/build_vcep_pvs1_exons.py \
        --gff data/MANE.GRCh38.v1.4.refseq_genomic.gff.gz \
        --out data/GRCh38/vcep_pvs1_exons.tsv

The MANE GFF3 is downloadable from
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/ .
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

# Clinical transcripts the VCEP PVS1 module reasons about. Keep in sync with
# ``vcep_pvs1.py`` (versionless match — the GFF may carry a different minor
# version than the spec).
TARGET_TRANSCRIPTS: dict[str, str] = {
    "RPE65": "NM_000329", "CYP1B1": "NM_000104", "VHL": "NM_000551",
    "GCK": "NM_000162", "RAG1": "NM_000448", "ATM": "NM_000051",
    "GP9": "NM_000174", "IDUA": "NM_000203", "ACVRL1": "NM_000020",
    "PAH": "NM_000277", "HNF1A": "NM_000545", "GJB2": "NM_004004",
    "FOXG1": "NM_005249", "DICER1": "NM_177438", "PALB2": "NM_024675",
    "FBN1": "NM_000138", "GP1BA": "NM_000173", "CDH1": "NM_004360",
    "AIPL1": "NM_014336", "ACADVL": "NM_000018", "TP53": "NM_000546",
    "GAA": "NM_000152", "GAMT": "NM_000156", "HNF4A": "NM_175914",
    "RUNX1": "NM_001754", "CDKL5": "NM_001323289", "RPGR": "NM_001034853",
    "IL2RG": "NM_000206", "MECP2": "NM_001110792", "F9": "NM_000133",
    "ABCD1": "NM_000033",
    # cspec re-examination batch
    "ADA": "NM_000022", "DCLRE1C": "NM_001033855", "JAK3": "NM_000215",
    "IL7R": "NM_002185", "FOXN1": "NM_001369369", "RAG2": "NM_000536",
    "CTLA4": "NM_005214", "KCNQ1": "NM_000218", "MLH1": "NM_000249",
    "MSH2": "NM_000251", "MSH6": "NM_000179", "PMS2": "NM_000535",
    "OTC": "NM_000531", "SLC9A6": "NM_001379110", "TCF4": "NM_001083962",
    "UBE3A": "NM_130839", "GUCY2D": "NM_000180", "RS1": "NM_000330",
    # decision-tree batch (files supplied)
    "ENG": "NM_001114753", "GP1BB": "NM_000407", "SCN1B": "NM_001037",
    "SCN2A": "NM_001371246", "SCN3A": "NM_006922", "SCN8A": "NM_014191",
    "NEB": "NM_001164508", "F8": "NM_000132", "PTEN": "NM_000314",
    "MYBPC3": "NM_000256", "HBB": "NM_000518", "HBA2": "NM_000517",
    # GALT decision-tree (GN158): in-frame exon 6/7/9/10 skip strengths.
    "GALT": "NM_000155",
}


def _attr(col9: str, key: str) -> str | None:
    m = re.search(rf"{key}=([^;]+)", col9)
    return m.group(1) if m else None


def _open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def _is_primary(seqid: str) -> bool:
    """Primary-assembly chromosome? GRCh38 MANE uses UCSC 'chr*'; GRCh37 RefSeq
    uses 'NC_*'. Alt loci / patches (NW_*, NT_*) and unplaced scaffolds are
    rejected so a transcript annotated on several scaffolds is not double-counted."""
    return seqid.startswith("chr") or seqid.startswith("NC_")


def build(gff: Path) -> list[dict]:
    want = {v: g for g, v in TARGET_TRANSCRIPTS.items()}
    # (transcript_versionless, seqid) -> {strand, exons, cds}. The GRCh37 file is
    # NOT MANE-limited: a clinical transcript may appear on the primary
    # chromosome AND on alt scaffolds / PAR_Y, so loci are kept separate and one
    # is chosen per transcript below.
    loci: dict[tuple[str, str], dict] = defaultdict(
        lambda: {"strand": "+", "exons": set(), "cds": set()}
    )

    with _open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            parent = _attr(f[8], "transcript_id") or _attr(f[8], "Parent") or ""
            base = parent.split(".")[0].replace("rna-", "")
            if base not in want:
                continue
            rec = loci[(base, f[0])]
            rec["strand"] = f[6]
            rec[("exons" if f[2] == "exon" else "cds")].add((int(f[3]), int(f[4])))

    # Choose one locus per transcript: prefer a primary chromosome, then the one
    # with the most coding bases (guards against truncated alt-scaffold copies).
    best: dict[str, dict] = {}
    for (base, seqid), rec in loci.items():
        if not rec["cds"]:
            continue
        cds_bases = sum(e - s + 1 for (s, e) in rec["cds"])
        key = (_is_primary(seqid), cds_bases)
        if base not in best or key > best[base]["_key"]:
            best[base] = {"_key": key, "seqid": seqid, **rec}

    rows: list[dict] = []
    for base, rec in best.items():
        gene = want[base]
        exons = sorted(rec["exons"])
        cds = sorted(rec["cds"])
        if not exons or not cds:
            continue
        # Order exons 5'→3' (VEP numbering): ascending on +, descending on −.
        if rec["strand"] == "-":
            exons = exons[::-1]
        cds_ranges = cds
        # Coding length contributed by each exon = overlap with any CDS segment.
        per_exon_cds: list[int] = []
        for (es, ee) in exons:
            c = sum(max(0, min(ee, ce) - max(es, cs) + 1) for (cs, ce) in cds_ranges)
            per_exon_cds.append(c)
        total_cds = sum(per_exon_cds)
        coding_idx = [i for i, c in enumerate(per_exon_cds) if c > 0]
        last_coding = coding_idx[-1]
        penult_coding = coding_idx[-2] if len(coding_idx) >= 2 else -1

        cum = 0
        n_exons = len(exons)
        for i, (es, ee) in enumerate(exons):
            c = per_exon_cds[i]
            if c == 0:
                continue
            rows.append({
                "gene": gene,
                "transcript": base,
                "n_exons": n_exons,
                "exon": i + 1,
                "cds_len_nt": c,
                "cum_cds_before_nt": cum,
                "inframe_skip": int(c % 3 == 0),
                "pct_protein": round(c / total_cds, 4),
                "is_last_coding": int(i == last_coding),
                "is_penultimate_coding": int(i == penult_coding),
            })
            cum += c
    rows.sort(key=lambda r: (r["gene"], r["exon"]))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff", required=True, type=Path, help="MANE RefSeq GFF3 (.gff/.gff.gz)")
    ap.add_argument("--out", required=True, type=Path, help="output TSV")
    args = ap.parse_args()

    rows = build(args.gff)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["gene", "transcript", "n_exons", "exon", "cds_len_nt",
            "cum_cds_before_nt", "inframe_skip", "pct_protein",
            "is_last_coding", "is_penultimate_coding"]
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    genes = sorted({r["gene"] for r in rows})
    print(f"Wrote {len(rows)} coding-exon rows for {len(genes)} genes to {args.out}")
    missing = sorted(set(TARGET_TRANSCRIPTS) - set(genes))
    if missing:
        print(f"WARNING: no rows for {len(missing)} gene(s): {', '.join(missing)}")


if __name__ == "__main__":
    main()
