"""Build ``resources/shared/mane_exons.tsv`` — transcript structures for NMD prediction.

The PVS1 decision tree needs to know where the premature termination codon (PTC)
lies relative to the last exon-exon junction (NMD is predicted when the PTC is
more than 50 nt upstream of the last junction; Abou Tayoun et al. 2018). VEP
reports only the exon number of the *variant*, which is not enough: a frameshift
can create its stop codon several exons downstream, and the 50-nt rule needs
transcript coordinates. This script extracts, for every MANE Select / MANE Plus
Clinical transcript, the exon lengths (5'->3', in transcript orientation) and the
cDNA position of the first coding nucleotide.

Source: the MANE GFF (``resources/gff/mane_phase16_rs.ucsc_seqids.gff``, NCBI MANE
Release v1.5; the GFF itself is not version-controlled). Transcript structure is
assembly-independent, so the same table serves GRCh37 and GRCh38.

Output columns::

    refseq  ensembl  gene  mane_tag  strand  cds_start  cds_end  exon_lengths
    NM_000314.8  ENST00000371953.8  PTEN  MANE Select  +  1358  2569  1361,85,45,44,239,142,167,223,2648

``cds_start`` / ``cds_end`` are 1-based cDNA positions of the first / last coding
nucleotide (stop codon included). ``exon_lengths`` is comma-separated.

Usage::

    python scripts/build_mane_exons.py \
        --gff resources/gff/mane_phase16_rs.ucsc_seqids.gff \
        --out resources/shared/mane_exons.tsv
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

_ID = re.compile(r"(?:^|;)ID=([^;]+)")
_PARENT = re.compile(r"(?:^|;)Parent=([^;]+)")
_GENE = re.compile(r"(?:^|;)gene=([^;]+)")
_ENST = re.compile(r"Ensembl:(ENST[0-9.]+)")
_TAG = re.compile(r"(?:^|;)tag=([^;]+)")


def build(gff: Path) -> list[dict]:
    mrna: dict[str, dict] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    cds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    annotation = ""
    with gff.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#!annotation-source"):
                annotation = line.split(" ", 1)[1].strip()
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype, start, end, strand, attrs = f[2], int(f[3]), int(f[4]), f[6], f[8]
            if ftype == "mRNA":
                tag = _TAG.search(attrs)
                if not tag or "MANE" not in tag.group(1):
                    continue
                rid = _ID.search(attrs).group(1)
                enst = _ENST.search(attrs)
                gene = _GENE.search(attrs)
                mrna[rid] = dict(
                    refseq=rid.removeprefix("rna-"),
                    ensembl=enst.group(1) if enst else "",
                    gene=gene.group(1) if gene else "",
                    mane_tag=tag.group(1),
                    strand=strand,
                )
            elif ftype in ("exon", "CDS"):
                par = _PARENT.search(attrs)
                if par:
                    (exons if ftype == "exon" else cds)[par.group(1)].append((start, end))
    rows = []
    for rid, info in mrna.items():
        ex = sorted(exons.get(rid, []))
        cd = sorted(cds.get(rid, []))
        if not ex or not cd:
            continue
        rev = info["strand"] == "-"
        if rev:
            ex = ex[::-1]
        lengths = [e - s + 1 for s, e in ex]
        # cDNA position of a genomic coordinate (transcript orientation).
        def cdna(pos: int) -> int | None:
            acc = 0
            for (s, e), ln in zip(ex, lengths):
                if s <= pos <= e:
                    return acc + ((e - pos) if rev else (pos - s)) + 1
                acc += ln
            return None
        first = cd[-1][1] if rev else cd[0][0]
        last = cd[0][0] if rev else cd[-1][1]
        cs, ce = cdna(first), cdna(last)
        if cs is None or ce is None:
            continue
        rows.append(dict(info, cds_start=cs, cds_end=ce,
                         exon_lengths=",".join(map(str, lengths))))
    return rows, annotation


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--gff", type=Path, default=Path("resources/gff/mane_phase16_rs.ucsc_seqids.gff"))
    ap.add_argument("--out", type=Path, default=Path("resources/shared/mane_exons.tsv"))
    a = ap.parse_args()
    rows, annotation = build(a.gff)
    cols = ["refseq", "ensembl", "gene", "mane_tag", "strand", "cds_start", "cds_end", "exon_lengths"]
    with a.out.open("w", encoding="utf-8") as fh:
        fh.write(f"# source: {a.gff.name} ({annotation})\n")
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["gene"], r["refseq"])):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {len(rows)} transcripts -> {a.out}")


if __name__ == "__main__":
    main()
