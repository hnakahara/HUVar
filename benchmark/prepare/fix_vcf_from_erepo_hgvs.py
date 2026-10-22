#!/usr/bin/env python3
"""Fix eRepo benchmark VCF records that TransVar converted incorrectly.

TransVar (HGVS c. -> VCF) mis-converts some eRepo variants: about 20% of the
``delins`` records (e.g. HNF4A c.281_282delinsC became a G>C SNV, PALB2
c.2524_2535delinsTCAGA an in-frame 12->6 change), and a few variants placed with
another transcript's numbering (MECP2 c.916C>T, GCK c.772G>T, UTR/intronic
variants). Each eRepo record carries ClinGen Allele Registry genomic HGVS
(``HGVS Expressions`` in erepo-tabbed.tsv: NC_ accessions for GRCh38 / GRCh37 /
NCBI36). This script compares every VCF record with that genomic HGVS and
rewrites the mismatching records from it, using the reference FASTA for the
REF / anchor bases.

The variant_id (chrom:pos:ref:alt) of a rewritten record changes, so the same
old->new mapping is applied to the HUVar eRepo resource tables keyed by
variant_id (erepo_classifications_*.tsv, erepo_manual_criteria_*.tsv) given with
--remap. Every file that is modified is first copied to <file>.bak_vcffix.
Line endings (CRLF / LF) are preserved.

Usage (GRCh38):
  python fix_vcf_from_erepo_hgvs.py --assembly 38 \
      --vcf inputs/erepo_hg38_may_2026.vcf --fasta /path/GRCh38.fa \
      --tabbed /path/erepo-tabbed.tsv --original inputs/erepo-original.csv \
      --remap <HUVar>/resources/GRCh38/erepo_classifications_hg38.tsv \
              <HUVar>/resources/GRCh38/erepo_manual_criteria_hg38.tsv \
      --log vcffix_hg38.tsv
(GRCh37: --assembly 37 with the hg19 VCF / FASTA / resource files.)
Add --dry-run to only report.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)


# --------------------------------------------------------------------------- FASTA
class Fasta:
    """Minimal indexed FASTA reader (.fai required; samtools faidx <fa>)."""

    def __init__(self, path: str) -> None:
        self.path = path
        self.idx: dict[str, tuple[int, int, int, int]] = {}
        with open(path + ".fai") as fh:
            for ln in fh:
                f = ln.split("\t")
                self.idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        self.fh = open(path, "rb")

    def _name(self, chrom: str) -> str:
        for c in (chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom,
                  "MT" if chrom in ("chrM", "M") else chrom):
            if c in self.idx:
                return c
        raise KeyError(chrom)

    def seq(self, chrom: str, start: int, end: int) -> str:
        """1-based inclusive."""
        length, offset, lb, lw = self.idx[self._name(chrom)]
        out = []
        for p in range(start - 1, end):
            self.fh.seek(offset + (p // lb) * lw + p % lb)
            out.append(self.fh.read(1).decode())
        return "".join(out).upper()


# --------------------------------------------------------------------------- HGVS
_CHR = {f"NC_0000{n:02d}": (f"chr{n}" if n <= 22 else {23: "chrX", 24: "chrY"}[n])
        for n in range(1, 25)}


def genomic_hgvs(tabbed: str, original: str, assembly: str) -> dict:
    """{(gene, c.): (chrom, g.)} for the eRepo records, using the eRepo record's
    own RefSeq transcript (erepo-original) to avoid c. collisions between
    transcripts. GRCh38 = highest NC_ version, GRCh37 = second highest."""
    nm: dict = {}
    with open(original, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            nm.setdefault((r["Gene"], r["HGVSc"]), set()).add(r["RefSeq (NM_)"])
    gmap: dict = {}
    with open(tabbed, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ex = [x.strip() for x in (r.get("HGVS Expressions") or "").split(",")]
            nc = [x for x in ex if re.match(r"NC_0000\d\d\.\d+:g\.", x)]
            if not nc:
                continue
            vers = sorted({int(x.split(":")[0].split(".")[1]) for x in nc}, reverse=True)
            k = 0 if assembly == "38" else 1
            if len(vers) <= k:
                continue
            g = [x for x in nc if int(x.split(":")[0].split(".")[1]) == vers[k]][0]
            acc, gh = g.split(":", 1)
            chrom = _CHR.get(acc.split(".")[0])
            for x in ex:
                m = re.match(r"(NM_\d+)\.\d+:(c\.\S+)", x)
                if m and chrom:
                    gmap.setdefault((r["HGNC Gene Symbol"], m.group(1), m.group(2)), set()).add((chrom, gh))
    out = {}
    for (gene, c), nms in nm.items():
        gs = set()
        for n in nms:
            gs |= gmap.get((gene, n, c), set())
        if len(gs) == 1:
            out[(gene, c)] = next(iter(gs))
    return out


def g_to_vcf(chrom: str, g: str, fa: Fasta) -> tuple[int, str, str] | None:
    m = re.fullmatch(r"g\.(\d+)([ACGT])>([ACGT])", g)
    if m:
        pos = int(m.group(1))
        ref = fa.seq(chrom, pos, pos)
        if ref != m.group(2):
            raise ValueError(f"reference mismatch {g}: genome has {ref}")
        return pos, ref, m.group(3)
    m = re.fullmatch(r"g\.(\d+)(?:_(\d+))?(del|dup|ins|delins)([ACGT]*)", g)
    if not m:
        return None
    a = int(m.group(1)); b = int(m.group(2) or a); typ = m.group(3); s = m.group(4)
    if typ == "del":
        ref = fa.seq(chrom, a - 1, b)
        return a - 1, ref, ref[0]
    if typ == "dup":
        anc = fa.seq(chrom, b, b)
        return b, anc, anc + fa.seq(chrom, a, b)
    if typ == "ins":
        anc = fa.seq(chrom, a, a)
        return a, anc, anc + s
    ref = fa.seq(chrom, a, b)                       # delins
    if len(ref) == 1 and len(s) == 1:
        return a, ref, s
    if ref[0] == s[:1]:
        return a, ref, s
    anc = fa.seq(chrom, a - 1, a - 1)               # anchor base keeps the record VCF-style
    return a - 1, anc + ref, anc + s


def matches(pos: int, ref: str, alt: str, g: str) -> bool:
    m = re.fullmatch(r"g\.(\d+)([ACGT])>([ACGT])", g)
    if m:
        return len(ref) == 1 and len(alt) == 1 and pos == int(m.group(1)) and ref == m.group(2) and alt == m.group(3)
    m = re.fullmatch(r"g\.(\d+)(?:_(\d+))?(del|dup|ins|delins)([ACGT]*)", g)
    if not m:
        return True                                 # unparsed -> leave untouched
    a = int(m.group(1)); b = int(m.group(2) or a); typ = m.group(3); s = m.group(4)
    net = {"del": -(b - a + 1), "dup": b - a + 1, "ins": len(s), "delins": len(s) - (b - a + 1)}[typ]
    return (len(alt) - len(ref)) == net and abs(pos - a) <= 60


# --------------------------------------------------------------------------- IO
def read_lines(path: Path) -> tuple[list[str], str]:
    raw = path.read_bytes().decode("utf-8")
    eol = "\r\n" if "\r\n" in raw else "\n"
    return raw.split(eol), eol


def backup(path: Path) -> None:
    bak = path.with_name(path.name + ".bak_vcffix")
    if not bak.exists():
        shutil.copy2(path, bak)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly", choices=["38", "37"], required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--tabbed", required=True)
    ap.add_argument("--original", required=True)
    ap.add_argument("--remap", nargs="*", default=[])
    ap.add_argument("--log", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    fa = Fasta(args.fasta)
    gmap = genomic_hgvs(args.tabbed, args.original, args.assembly)
    vcf = Path(args.vcf)
    lines, eol = read_lines(vcf)
    mapping: dict[str, str] = {}
    log = []
    new_lines = []
    for ln in lines:
        if not ln or ln.startswith("#"):
            new_lines.append(ln)
            continue
        f = ln.split("\t")
        chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
        gene, c = vid.split(":", 1)
        hit = gmap.get((gene, c))
        if not hit:
            new_lines.append(ln)
            continue
        gchrom, g = hit
        if gchrom != chrom or matches(pos, ref, alt, g):
            new_lines.append(ln)
            continue
        try:
            fixed = g_to_vcf(chrom, g, fa)
        except Exception as e:                       # noqa: BLE001
            log.append((vid, f"{chrom}:{pos}:{ref}:{alt}", "", g, f"ERROR {e}"))
            new_lines.append(ln)
            continue
        if fixed is None:
            log.append((vid, f"{chrom}:{pos}:{ref}:{alt}", "", g, "unparsed g."))
            new_lines.append(ln)
            continue
        npos, nref, nalt = fixed
        old_id = f"{chrom}:{pos}:{ref}:{alt}"
        new_id = f"{chrom}:{npos}:{nref}:{nalt}"
        mapping[old_id] = new_id
        f[1], f[3], f[4] = str(npos), nref, nalt
        new_lines.append("\t".join(f))
        log.append((vid, old_id, new_id, g, "fixed"))

    print(f"GRCh{args.assembly}: {sum(1 for x in log if x[4] == 'fixed')} records fixed, "
          f"{sum(1 for x in log if x[4] != 'fixed')} not fixed")
    for x in log:
        print("\t".join(x))
    if args.log:
        with open(args.log, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["erepo_variant", "old_variant_id", "new_variant_id", "erepo_g_hgvs", "status"])
            w.writerows(log)
    if args.dry_run:
        return

    # VCF: keep the original record order (callers sort if needed).
    backup(vcf)
    vcf.write_bytes(eol.join(new_lines).encode("utf-8"))
    for p in args.remap:
        path = Path(p)
        rl, reol = read_lines(path)
        n = 0
        out = []
        for ln in rl:
            if ln and not ln.startswith("#"):
                head, sep, rest = ln.partition("\t")
                if head in mapping:
                    ln = mapping[head] + sep + rest
                    n += 1
            out.append(ln)
        backup(path)
        path.write_bytes(reol.join(out).encode("utf-8"))
        print(f"{path}: {n} rows remapped")


if __name__ == "__main__":
    main()
