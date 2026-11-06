"""Translate MANE Select CDS from the genome FASTA + MANE GFF (helper for
scripts/build_kcnq1_paralog_conservation.py).

usage: extract_mane_proteins.py GENOME.fa MANE.gff NM_x.y [NM_...] > proteins.fa
"""
from __future__ import annotations

import sys

import pysam

_CODE = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_B):
    for j, b in enumerate(_B):
        for k, c in enumerate(_B):
            _CODE[a + b + c] = _AA[16 * i + 4 * j + k]
_COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def cds_segments(gff: str, tx: str):
    segs, chrom, strand, gene = [], None, None, None
    with open(gff, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS" or f"Parent=rna-{tx};" not in f[8]:
                continue
            chrom, strand = f[0], f[6]
            segs.append((int(f[3]), int(f[4])))
            gene = f[8].split("gene=")[1].split(";")[0]
    return chrom, strand, sorted(segs), gene


def main() -> None:
    fa = pysam.FastaFile(sys.argv[1])
    for tx in sys.argv[3:]:
        chrom, strand, segs, gene = cds_segments(sys.argv[2], tx)
        if chrom not in fa.references and chrom.startswith("chr"):
            chrom = chrom[3:]                       # Ensembl-style FASTA (1, 2, ..., X)
        seq = "".join(fa.fetch(chrom, s - 1, e) for s, e in segs).upper()
        if strand == "-":
            seq = seq.translate(_COMP)[::-1]
        prot = "".join(_CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
        print(f">{gene} {tx}\n{prot.rstrip('*')}")


if __name__ == "__main__":
    main()
