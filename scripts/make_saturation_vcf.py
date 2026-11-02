#!/usr/bin/env python3
"""Write a VCF of every possible SNV in the coding region and splice regions of
the MANE Select transcript of the given genes (saturation test set).

Regions (per transcript, GRCh38, from the MANE GFF):

* ``cds``          — every coding base (CDS features, stop codon included);
* ``exon_splice``  — the ``--exon-splice`` exonic bases next to each splice site
                     (default 3, VEP's splice_region) where they are not already
                     coding, i.e. at UTR-exon boundaries;
* ``intron_flank`` — ``--flank`` intronic bases on each side of every intron
                     (default 20: donor +1..+20, acceptor -20..-1; covers the
                     canonical +/-1,2 sites, the donor region and the
                     polypyrimidine tract).

At each position all three alternative alleles are written, so the file contains
3 x (number of positions) SNVs. Indels are not enumerated.

The reference allele comes from an indexed FASTA (``--fasta``; HUVar's
``<data>/GRCh38/genome/GRCh38.p14.fa``). Contig names with or without the ``chr``
prefix are both accepted; the VCF uses ``chr``-prefixed names like the eRepo VCF.

Because the MANE GFF is large (and not part of the repository), the regions can be
written once with ``--write-regions`` and reused elsewhere with ``--regions``:

    # on a machine with the GFF: extract the regions only
    python scripts/make_saturation_vcf.py --gff resources/gff/mane_phase16_rs.ucsc_seqids.gff \\
        --genes BRCA1 BRCA2 --write-regions inputs/brca12_saturation_regions_hg38.tsv

    # on the server with the reference FASTA: write the VCF
    python scripts/make_saturation_vcf.py --regions inputs/brca12_saturation_regions_hg38.tsv \\
        --fasta /data/GRCh38/genome/GRCh38.p14.fa -o brca12_saturation_hg38.vcf
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

_BASES = "ACGT"


def _attrs(col9: str) -> dict:
    out = {}
    for item in col9.strip().split(";"):
        k, sep, v = item.partition("=")
        if sep:
            out[k] = v
    return out


def read_mane_gff(gff: Path, genes: set[str]) -> dict:
    """{gene: {transcript, chrom, strand, exons:[(s,e)], cds:[(s,e)]}} (1-based, inclusive)
    for the MANE Select transcript of each requested gene."""
    tx: dict[str, dict] = {}
    with gff.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            a = _attrs(f[8])
            gene = a.get("gene")
            if gene not in genes or "MANE Select" not in a.get("tag", ""):
                continue
            tid = a.get("transcript_id") or a.get("Parent", "").replace("rna-", "")
            t = tx.setdefault(gene, {"transcript": tid, "chrom": f[0], "strand": f[6],
                                     "exons": [], "cds": []})
            if t["transcript"] != tid:
                raise SystemExit(f"{gene}: more than one MANE Select transcript ({t['transcript']}, {tid})")
            t["exons" if f[2] == "exon" else "cds"].append((int(f[3]), int(f[4])))
    missing = genes - set(tx)
    if missing:
        raise SystemExit(f"no MANE Select transcript in the GFF for: {', '.join(sorted(missing))}")
    for t in tx.values():
        t["exons"].sort()
        t["cds"].sort()
    return tx


def regions_from_transcripts(tx: dict, flank: int, exon_splice: int) -> list[tuple]:
    """[(chrom, start, end, gene, transcript, kind, strand)] 1-based inclusive, merged per kind."""
    out = []
    for gene, t in sorted(tx.items()):
        chrom, tid, strand = t["chrom"], t["transcript"], t["strand"]
        coding = set()
        for s, e in t["cds"]:
            out.append((chrom, s, e, gene, tid, "cds", strand))
            coding.update(range(s, e + 1))
        ex = t["exons"]
        for (s1, e1), (s2, e2) in zip(ex, ex[1:]):          # each intron e1+1 .. s2-1
            out.append((chrom, e1 + 1, min(e1 + flank, s2 - 1), gene, tid, "intron_flank", strand))
            out.append((chrom, max(s2 - flank, e1 + 1), s2 - 1, gene, tid, "intron_flank", strand))
            for s, e in ((max(s1, e1 - exon_splice + 1), e1), (s2, min(e2, s2 + exon_splice - 1))):
                extra = [p for p in range(s, e + 1) if p not in coding]
                if extra:
                    out.append((chrom, extra[0], extra[-1], gene, tid, "exon_splice", strand))
    # merge overlapping intervals of the same kind (short introns)
    merged: list[tuple] = []
    for r in sorted(out, key=lambda r: (r[0], r[3], r[5], r[1])):
        if merged and merged[-1][0] == r[0] and merged[-1][3] == r[3] and merged[-1][5] == r[5] \
                and r[1] <= merged[-1][2] + 1:
            m = merged[-1]
            merged[-1] = (m[0], m[1], max(m[2], r[2]), m[3], m[4], m[5], m[6])
        else:
            merged.append(r)
    return sorted(merged, key=lambda r: (r[0], r[1]))


def write_regions(regions: list[tuple], path: Path) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("#chrom\tstart0\tend\tgene\ttranscript\tregion\tstrand\n")   # BED: 0-based start
        for c, s, e, g, t, k, st in regions:
            fh.write(f"{c}\t{s - 1}\t{e}\t{g}\t{t}\t{k}\t{st}\n")


def read_regions(path: Path) -> list[tuple]:
    out = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            c, s0, e, g, t, k = f[:6]
            out.append((c, int(s0) + 1, int(e), g, t, k, f[6] if len(f) > 6 else "."))
    return out


def _contig(fasta, chrom: str) -> str:
    names = set(fasta.references)
    for c in (chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom):
        if c in names:
            return c
    raise SystemExit(f"contig {chrom} not in the FASTA")


_COMP = str.maketrans("ACGT", "TGCA")


def _check_codons(regions: list[tuple], fasta) -> None:
    """Sanity check of the FASTA against the transcript model: each CDS must start
    with ATG and end with a stop codon. A wrong assembly or contig naming shows up
    here instead of as a VCF of wrong REF alleles."""
    cds: dict = defaultdict(list)
    for c, s, e, g, t, k, st in regions:
        if k == "cds":
            cds[(g, c, st)].append((s, e))
    bad = []
    for (g, c, st), iv in cds.items():
        if st not in "+-":
            continue
        ctg = _contig(fasta, c)
        iv.sort()
        seq = "".join(fasta.fetch(ctg, s - 1, e).upper() for s, e in iv)
        if st == "-":
            seq = seq.translate(_COMP)[::-1]
        start, stop = seq[:3], seq[-3:]
        status = "ok" if start == "ATG" and stop in ("TAA", "TAG", "TGA") and len(seq) % 3 == 0 else "MISMATCH"
        print(f"  check {g}: CDS {len(seq)} bp, start {start}, stop {stop} -> {status}")
        if status != "ok":
            bad.append(g)
    if bad:
        raise SystemExit(f"reference FASTA does not match the GRCh38 transcript model for "
                         f"{', '.join(bad)} (wrong assembly or FASTA?)")


def write_vcf(regions: list[tuple], fasta_path: Path, out: Path) -> dict:
    import pysam  # HUVar dependency
    fasta = pysam.FastaFile(str(fasta_path))
    # one region label per position (cds > exon_splice > intron_flank)
    rank = {"cds": 0, "exon_splice": 1, "intron_flank": 2}
    _check_codons(regions, fasta)
    pos_info: dict[tuple, tuple] = {}
    for c, s, e, g, t, k, _st in regions:
        for p in range(s, e + 1):
            cur = pos_info.get((c, p))
            if cur is None or rank[k] < rank[cur[2]]:
                pos_info[(c, p)] = (g, t, k)
    counts: dict = defaultdict(int)
    chrom_order = sorted({c for c, _ in pos_info},
                         key=lambda c: (len(c.replace("chr", "")), c.replace("chr", "")))
    with out.open("w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##source=HUVar make_saturation_vcf.py; reference={fasta_path.name}\n")
        fh.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene (MANE Select)">\n')
        fh.write('##INFO=<ID=TX,Number=1,Type=String,Description="MANE Select transcript">\n')
        fh.write('##INFO=<ID=REGION,Number=1,Type=String,Description="cds / exon_splice / intron_flank">\n')
        for c in chrom_order:
            fh.write(f"##contig=<ID={c if c.startswith('chr') else 'chr' + c}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c in chrom_order:
            ctg = _contig(fasta, c)
            vcf_chrom = c if c.startswith("chr") else "chr" + c
            for p in sorted(p for cc, p in pos_info if cc == c):
                ref = fasta.fetch(ctg, p - 1, p).upper()
                if ref not in _BASES:
                    counts["skipped_N"] += 1
                    continue
                g, t, k = pos_info[(c, p)]
                for alt in _BASES:
                    if alt == ref:
                        continue
                    fh.write(f"{vcf_chrom}\t{p}\t{g}_{vcf_chrom}_{p}_{ref}_{alt}\t{ref}\t{alt}\t.\tPASS\t"
                             f"GENE={g};TX={t};REGION={k}\n")
                    counts[(g, k)] += 1
    return counts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--gff", type=Path, help="MANE GFF (e.g. resources/gff/mane_phase16_rs.ucsc_seqids.gff)")
    src.add_argument("--regions", type=Path, help="regions BED written earlier with --write-regions")
    ap.add_argument("--genes", nargs="+", default=["BRCA1", "BRCA2"])
    ap.add_argument("--flank", type=int, default=20, help="intronic bases per splice site (default 20)")
    ap.add_argument("--exon-splice", type=int, default=3,
                    help="exonic bases per splice site added at non-coding exon boundaries (default 3)")
    ap.add_argument("--write-regions", type=Path, help="write the regions BED here")
    ap.add_argument("--fasta", type=Path, help="indexed reference FASTA (GRCh38) for the REF alleles")
    ap.add_argument("-o", "--output", type=Path, help="output VCF")
    a = ap.parse_args()

    if a.gff:
        regions = regions_from_transcripts(read_mane_gff(a.gff, set(a.genes)), a.flank, a.exon_splice)
    else:
        regions = [r for r in read_regions(a.regions) if r[3] in set(a.genes)]
        if not regions:
            ap.error(f"no regions for {', '.join(a.genes)} in {a.regions}")
    if a.write_regions:
        write_regions(regions, a.write_regions)
        n = defaultdict(int)
        for c, s, e, g, t, k, _st in regions:
            n[(g, k)] += e - s + 1
        for (g, k), v in sorted(n.items()):
            print(f"  {g:8s} {k:13s} {v:6d} bp")
        print(f"regions -> {a.write_regions}")
    if a.output:
        if not a.fasta:
            ap.error("--fasta is required to write the VCF (REF alleles)")
        counts = write_vcf(regions, a.fasta, a.output)
        total = sum(v for k, v in counts.items() if isinstance(k, tuple))
        for (g, k), v in sorted((k, v) for k, v in counts.items() if isinstance(k, tuple)):
            print(f"  {g:8s} {k:13s} {v:6d} SNVs")
        print(f"{total} SNVs -> {a.output}" + (f" (skipped {counts['skipped_N']} N positions)"
                                              if counts.get("skipped_N") else ""))
    if not a.write_regions and not a.output:
        ap.error("nothing to do: give --write-regions and/or -o/--fasta")
    return 0


if __name__ == "__main__":
    sys.exit(main())
