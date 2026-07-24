"""Build ``resources/shared/mane_select.tsv`` — a gene -> MANE Select transcript map.

MANE Select is a GRCh38-native concept: VEP's ``mane_select`` flag is populated
only in the GRCh38 cache. On GRCh37 no transcript carries the flag, so the
pipeline's MANE-first transcript selection cannot fire and a non-MANE transcript
(e.g. a longer alternative isoform) may be chosen — changing HGVS numbering and,
worse, VCEP criteria that are keyed to the MANE codon range (PVS1 range, PM1
hotspots, ...). This map lets the annotation layer recover the MANE-equivalent
transcript by RefSeq/Ensembl base accession even when VEP provides no flag.

Source: the committed MANE GFF (``resources/gff/mane_phase16_rs.ucsc_seqids.gff``,
NCBI MANE release). Only ``mRNA`` features tagged ``MANE Select`` are used.

Output format (consumed by ``local_db.mane_db.load_mane_map``)::

    gene_symbol<TAB>refseq<TAB>ensembl
    PTEN<TAB>NM_000314.8<TAB>ENST00000371953.8

Accessions are stored versioned; the loader matches by version-stripped base.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

# GFF attribute helpers -------------------------------------------------------
_GENE_RE = re.compile(r"(?:^|;)gene=([^;]+)")
_TXID_RE = re.compile(r"(?:^|;)transcript_id=([^;]+)")
_ENST_RE = re.compile(r"Ensembl:(ENST[0-9.]+)")


def _parse_attrs(attr: str) -> tuple[str, str, str] | None:
    """Extract (gene, refseq_nuc, ensembl_nuc) from a MANE Select mRNA line."""
    if "tag=MANE Select" not in attr:
        return None
    gene_m = _GENE_RE.search(attr)
    tx_m = _TXID_RE.search(attr)
    enst_m = _ENST_RE.search(attr)
    if not (gene_m and tx_m):
        return None
    gene = gene_m.group(1).strip()
    refseq = tx_m.group(1).strip()
    ensembl = enst_m.group(1).strip() if enst_m else ""
    # MANE Select is RefSeq (NM_/NR_); ignore anything else defensively.
    if not gene or not refseq.startswith(("NM_", "NR_")):
        return None
    return gene, refseq, ensembl


def build(gff_path: Path, out_path: Path) -> int:
    """Parse the MANE GFF and write the gene->transcript TSV. Returns row count."""
    seen: dict[str, tuple[str, str]] = {}
    with gff_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "mRNA":
                continue
            parsed = _parse_attrs(cols[8])
            if parsed is None:
                continue
            gene, refseq, ensembl = parsed
            # One MANE Select per gene; first wins (GFF lists each mRNA once).
            seen.setdefault(gene, (refseq, ensembl))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_symbol", "refseq", "ensembl"])
        for gene in sorted(seen):
            refseq, ensembl = seen[gene]
            w.writerow([gene, refseq, ensembl])
    return len(seen)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--gff",
        type=Path,
        default=Path("./resources/gff/mane_phase16_rs.ucsc_seqids.gff"),
        help="MANE GFF (default: resources/gff/mane_phase16_rs.ucsc_seqids.gff)",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=Path("./resources/shared/mane_select.tsv"),
        help="Output TSV (default: resources/shared/mane_select.tsv)",
    )
    args = ap.parse_args()
    n = build(args.gff, args.output)
    print(f"wrote {n} genes -> {args.output}")


if __name__ == "__main__":
    main()
