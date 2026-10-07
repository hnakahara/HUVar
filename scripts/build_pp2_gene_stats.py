"""Build ``resources/<assembly>/pp2_gene_stats.tsv`` — gnomAD common missense per gene.

The PP2 fallback heuristic (genes without a ClinGen VCEP PP2 decision) asks
whether a gene has "a low rate of benign missense variation". ClinVar alone
under-counts benign missense variation in genes that are rarely tested, so the
benign side is supplemented with gnomAD missense variants that are common enough
to meet the gene's own BS1 threshold (Genome Medicine revision, Reviewer 3).

For every MANE Select CDS this script

1. selects PASS gnomAD SNVs whose filtering allele frequency (FAF95 grpmax; the
   grpmax point AF for genes whose VCEP uses it, ``af_basis=popmax``) is at least
   the global BS1 floor,
2. annotates them with the local VEP install (offline cache, MANE Select pick),
3. keeps missense variants in the gene whose frequency meets that gene's BS1
   threshold (``disease_prevalence.tsv``; flat default 0.5%),

and writes ``gene  common_missense  bs1_threshold``. Simplifications relative to
the per-variant BS1 evaluator: male-only and non-cancer AF subsets and the
``bs1_exclude`` variant lists are not applied.

The table is bundled in ``resources/GRCh38/`` and ``resources/GRCh37/`` (gene-level
counts do not depend on the assembly, so the GRCh38 build is copied to GRCh37). A
copy placed in ``data/<assembly>/`` overrides the bundled one. When neither is
present, PP2 falls back to ClinVar B/LB counts only.

Usage (on the server holding the gnomAD DuckDB and the VEP cache)::

    python scripts/build_pp2_gene_stats.py --data-dir $DATA --assembly GRCh38 \
        --mane-gff resources/gff/mane_phase16_rs.ucsc_seqids.gff --workers 16
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from acmg_classifier.config import Config
from acmg_classifier.criteria.allele_frequency import _BS1_FLOOR, DiseaseThresholds
from acmg_classifier.local_db.gnomad_db import _pass_filter
from acmg_classifier.models.enums import Assembly

_PARENT = re.compile(r"(?:^|;)Parent=([^;]+)")
_ID = re.compile(r"(?:^|;)ID=([^;]+)")
_GENE = re.compile(r"(?:^|;)gene=([^;]+)")
_TAG = re.compile(r"(?:^|;)tag=([^;]+)")


def _strip_chr(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def mane_cds_intervals(gff: Path) -> list[tuple[str, int, int, str]]:
    """(chrom, start, end, gene) of MANE Select CDS segments (1-based, inclusive)."""
    select: dict[str, str] = {}
    cds: list[tuple[str, int, int, str]] = []
    with gff.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                tag = _TAG.search(f[8])
                if tag and "MANE Select" in tag.group(1):
                    g = _GENE.search(f[8])
                    select[_ID.search(f[8]).group(1)] = g.group(1) if g else ""
            elif f[2] == "CDS":
                par = _PARENT.search(f[8])
                if par:
                    cds.append((f[0], int(f[3]), int(f[4]), par.group(1)))
    return [(_strip_chr(c), s, e, select[p]) for c, s, e, p in cds if p in select]


def candidate_variants(duckdb_path: Path, intervals, floor: float):
    """{(chrom, pos, ref, alt): (faf95, popmax_af)} — PASS SNVs in the CDS
    intervals with FAF95 (or grpmax AF when FAF is missing) >= floor."""
    import duckdb

    con = duckdb.connect(str(duckdb_path), read_only=True)
    sample = con.execute("SELECT chrom FROM variants LIMIT 1").fetchone()
    with_chr = bool(sample) and str(sample[0]).lower().startswith("chr")
    con.execute("CREATE TEMP TABLE cds (chrom TEXT, s INTEGER, e INTEGER)")
    con.executemany("INSERT INTO cds VALUES (?,?,?)",
                    [((f"chr{c}" if with_chr else c), s, e) for c, s, e, _ in intervals])
    rows = con.execute(
        """
        SELECT DISTINCT v.chrom, v.pos, v.ref, v.alt, v.faf95_popmax, v.popmax_af, v.filters
        FROM variants v JOIN cds ON v.chrom = cds.chrom AND v.pos BETWEEN cds.s AND cds.e
        WHERE length(v.ref) = 1 AND length(v.alt) = 1
          AND coalesce(v.faf95_popmax, v.popmax_af, 0) >= ?
        """, [floor]).fetchall()
    con.close()
    out: dict[tuple, list] = {}
    for chrom, pos, ref, alt, faf, pmax, filters in rows:
        if not _pass_filter(filters):
            continue
        key = (_strip_chr(str(chrom)), int(pos), ref, alt)
        cur = out.setdefault(key, [None, None])
        if faf is not None:
            cur[0] = faf if cur[0] is None else max(cur[0], faf)
        if pmax is not None:
            cur[1] = pmax if cur[1] is None else max(cur[1], pmax)
    return out


def run_vep(cfg: Config, variants, workers: int, vep_cmd: str, tmp: Path) -> dict[str, tuple[str, str]]:
    """{variant id: (gene symbol, consequence)} using the MANE Select transcript."""
    vcf = tmp / "pp2_candidates.vcf"
    with vcf.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in sorted(variants, key=lambda k: (k[0], k[1])):
            fh.write(f"{chrom}\t{pos}\t{chrom}:{pos}:{ref}:{alt}\t{ref}\t{alt}\t.\t.\t.\n")
    out = tmp / "pp2_candidates.vep.tsv"
    subprocess.run([
        vep_cmd, "--input_file", str(vcf), "--output_file", str(out),
        "--format", "vcf", "--tab", "--cache", "--offline", "--cache_version", "111",
        "--dir_cache", str(cfg.vep_cache_dir), "--assembly", cfg.assembly.value,
        "--merged", "--fasta", str(cfg.genome_fasta), "--symbol", "--mane_select",
        "--pick", "--pick_order", "mane_select,canonical,rank",
        "--fields", "Uploaded_variation,SYMBOL,Consequence,MANE_SELECT",
        "--no_stats", "--no_progress", "--fork", str(workers), "--force_overwrite",
    ], check=True)
    res: dict[str, tuple[str, str]] = {}
    with out.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            vid, sym, csq, mane = (line.rstrip("\n").split("\t") + ["", "", "", ""])[:4]
            if mane and mane != "-":
                res[vid] = (sym, csq)
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--assembly", default="GRCh38", choices=["GRCh38", "GRCh37"])
    ap.add_argument("--mane-gff", type=Path,
                    default=Path("resources/gff/mane_phase16_rs.ucsc_seqids.gff"))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--vep-cmd", default=None, help="VEP executable (default: autodetect)")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()

    cfg = Config(data_dir=a.data_dir, assembly=Assembly(a.assembly))
    if a.assembly == "GRCh37":
        # MANE coordinates are GRCh38; GRCh37 runs reuse the GRCh38 table
        # (gene-level counts are assembly-independent).
        raise SystemExit("Build on GRCh38 and copy the TSV to resources/GRCh37/.")
    if a.vep_cmd is None:
        from acmg_classifier.setup.vep_installer import find_vep_cmd
        a.vep_cmd = find_vep_cmd()
    out_path = a.out or (Path(__file__).resolve().parents[1] / "resources"
                         / a.assembly / "pp2_gene_stats.tsv")

    thresholds = DiseaseThresholds(cfg.disease_prevalence_tsv)
    intervals = mane_cds_intervals(a.mane_gff)
    print(f"MANE Select CDS segments: {len(intervals)}")
    cands = candidate_variants(cfg.gnomad_duckdb, intervals, _BS1_FLOOR)
    print(f"gnomAD PASS SNV candidates (AF >= {_BS1_FLOOR}): {len(cands)}")

    with tempfile.TemporaryDirectory() as td:
        vep = run_vep(cfg, cands, a.workers, a.vep_cmd, Path(td))

    counts: dict[str, int] = defaultdict(int)
    genes = {g for *_, g in intervals if g}
    for (chrom, pos, ref, alt), (faf, pmax) in cands.items():
        hit = vep.get(f"{chrom}:{pos}:{ref}:{alt}")
        if not hit or "missense_variant" not in hit[1].split("&"):
            continue
        gene = hit[0]
        gt = thresholds.get(gene)
        freq = faf if faf is not None else pmax
        if gt.af_basis == "popmax" and cfg.popmax_af_basis and pmax is not None:
            freq = pmax
        if freq is not None and freq >= gt.bs1:
            counts[gene] += 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# built {_dt.date.today().isoformat()} from {cfg.gnomad_duckdb.name}, "
                 f"{a.mane_gff.name}, {cfg.disease_prevalence_tsv.name}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "common_missense", "bs1_threshold"])
        for g in sorted(genes):
            w.writerow([g, counts.get(g, 0), f"{thresholds.get(g).bs1:g}"])
    n_pos = sum(1 for g in genes if counts.get(g))
    print(f"wrote {len(genes)} genes ({n_pos} with >=1 common missense) -> {out_path}")


if __name__ == "__main__":
    main()
