#!/usr/bin/env python3
"""
scripts/setup_data.py

Automated setup of data files for acmg-classifier.
Existing files are skipped; only missing ones are downloaded/built.

Usage:
  python scripts/setup_data.py --data-dir ./data

  # Point at existing files to skip downloads
  python scripts/setup_data.py --data-dir ./data \\
      --genome-fasta /db/reference/GRCh38/hg38.fa \\
      --gnomad-vcf-dir /db/gnomad/v4.1/exomes/vcf

  # Skip gnomAD (~300 GB) and set up everything else
  python scripts/setup_data.py --data-dir ./data --skip-gnomad

  # Download gnomAD for specific chromosomes only
  python scripts/setup_data.py --data-dir ./data --gnomad-chromosomes chr1 chr2 chrX
"""
from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ENSEMBL_RELEASE = 111

# MMSplice gene-annotation GTF output names (must match config.mmsplice_gtf and
# downloader._MMSPLICE_GTF). Filtered to protein-coding genes during setup.
GTF_NAMES = {
    "GRCh38": "Homo_sapiens.GRCh38.111.protein_coding.gtf",
    "GRCh37": "Homo_sapiens.GRCh37.87.protein_coding.gtf",
}

URLS: dict[str, dict[str, str]] = {
    "GRCh38": {
        "genome": (
            f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        ),
        "vep_cache": (
            f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/variation/indexed_vep_cache/"
            f"homo_sapiens_merged_vep_{ENSEMBL_RELEASE}_GRCh38.tar.gz"
        ),
        "clinvar_vcf": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
        "clinvar_xml": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_xml_old_format/ClinVarFullRelease_00-latest.xml.gz",
        "alphamissense": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz",
        "esm1b_zip": (
            "https://huggingface.co/spaces/ntranoslab/esm_variants/resolve/main/"
            "ALL_hum_isoforms_ESM1b_LLR.zip"
        ),
        "gtf": (
            f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/gtf/homo_sapiens/"
            f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"
        ),
        "gnomad_constraint": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/"
            "gnomad.v4.1.constraint_metrics.tsv"
        ),
        "gnomad_vcf": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/"
            "gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz"
        ),
        "gnomad_vcf_tbi": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/"
            "gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz.tbi"
        ),
        "repeatmasker": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
    },
    "GRCh37": {
        "genome": (
            f"https://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
        ),
        "vep_cache": (
            f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/variation/indexed_vep_cache/"
            f"homo_sapiens_merged_vep_{ENSEMBL_RELEASE}_GRCh37.tar.gz"
        ),
        "clinvar_vcf": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
        "clinvar_xml": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_xml_old_format/ClinVarFullRelease_00-latest.xml.gz",
        "alphamissense": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz",
        "esm1b_zip": (
            "https://huggingface.co/spaces/ntranoslab/esm_variants/resolve/main/"
            "ALL_hum_isoforms_ESM1b_LLR.zip"
        ),
        "gtf": (
            "https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/"
            "Homo_sapiens.GRCh37.87.gtf.gz"
        ),
        "gnomad_constraint": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/"
            "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
        ),
        "gnomad_vcf": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/"
            "gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz"
        ),
        "gnomad_vcf_tbi": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/"
            "gnomad.exomes.r2.1.1.sites.{chrom}.vcf.bgz.tbi"
        ),
        "repeatmasker": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz",
    },
}

# GRCh38 (v4.1): with chr prefix / GRCh37 (v2.1.1): without
GNOMAD_CHROMS: dict[str, list[str]] = {
    "GRCh38": [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]],
    "GRCh37": [str(c) for c in list(range(1, 23)) + ["X", "Y"]],
}

# Candidate directories for existing gnomAD / genome files on the server
_GENOME_SEARCH = [
    "/db/reference",
    "/mnt/department/db/genome",
    "/opt/resources/genome",
    "/data/reference",
    "/shared/genome",
]
_GNOMAD_SEARCH = [
    "/db/gnomad",
    "/mnt/department/db/gnomad",
    "/opt/resources/gnomad",
    "/data/gnomad",
    "/shared/databases/gnomad",
]

_GENOME_NAMES = {
    "GRCh38": "GRCh38.p14.fa",
    "GRCh37": "GRCh37.p13.fa",
}
_GNOMAD_VER = {"GRCh38": "4.1", "GRCh37": "2.1.1"}

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _has(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _require(*cmds: str) -> None:
    missing = [c for c in cmds if not _has(c)]
    if missing:
        sys.exit(
            f"[ERROR] Required tools not found: {', '.join(missing)}\n"
            f"  Install example: conda install -c bioconda {' '.join(missing)}"
        )


def _run(cmd: list, **kw) -> None:
    print("  $", " ".join(shlex.quote(str(x)) for x in cmd))
    subprocess.run([str(x) for x in cmd], check=True, **kw)


def _run_shell(cmd: str) -> None:
    print(f"  $ {cmd}")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def _download(url: str, dest: Path, label: str = "") -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  ↓  {label or dest.name}")
    if _has("wget"):
        _run(["wget", "-c", "--no-verbose", "--show-progress", "-O", str(dest), url])
    elif _has("curl"):
        _run(["curl", "-C", "-", "-L", "--progress-bar", "-o", str(dest), url])
    else:
        sys.exit("[ERROR] wget or curl is required")


def _verify_size(path: Path, min_bytes: int, label: str = "") -> bool:
    """Treat a file far smaller than expected as corrupt and delete it."""
    if not path.exists():
        return False
    sz = path.stat().st_size
    if sz < min_bytes:
        print(f"  [WARN] {label or path.name} is too small ({sz:,} bytes < {min_bytes:,})")
        print(f"         Corrupt or incomplete — deleting and re-downloading")
        path.unlink(missing_ok=True)
        return False
    return True


def _verify_sqlite_has_rows(path: Path, table: str = "variants") -> bool:
    """True if the SQLite file has the table with rows; delete it and return False if empty."""
    if not path.exists():
        return False
    try:
        import sqlite3
        con = sqlite3.connect(str(path))
        try:
            cur = con.execute(f"SELECT COUNT(*) FROM {table}")
            count = cur.fetchone()[0]
        finally:
            con.close()
        if count > 0:
            return True
        print(f"  [WARN] {path.name} is empty (0 rows in {table}) — deleting and rebuilding")
    except Exception as e:
        print(f"  [WARN] {path.name} is corrupt ({e}) — deleting and rebuilding")
    path.unlink(missing_ok=True)
    return False


def _verify_duckdb_has_rows(path: Path, table: str = "variants") -> bool:
    """True if the DuckDB file has the table with rows; delete it and return False if empty."""
    if not path.exists():
        return False
    try:
        import duckdb
        con = duckdb.connect(str(path), read_only=True)
        try:
            count = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        finally:
            con.close()
        if count > 0:
            return True
        print(f"  [WARN] {path.name} is empty (0 rows in {table}) — deleting and rebuilding")
    except Exception as e:
        print(f"  [WARN] {path.name} is corrupt ({e}) — deleting and rebuilding")
    path.unlink(missing_ok=True)
    return False


def _find_existing(search_dirs: list[str], patterns: list[str]) -> Path | None:
    for d in search_dirs:
        p = Path(d)
        if not p.is_dir():
            continue
        for pat in patterns:
            for f in p.rglob(pat):
                print(f"  Found existing file: {f}")
                return f
    return None


def _add_src_to_path() -> None:
    src = Path(__file__).parent.parent / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


# ---------------------------------------------------------------------------
# Setup steps
# ---------------------------------------------------------------------------


def step_genome(asm_dir: Path, assembly: str, urls: dict, genome_fasta: Path | None, skip: bool) -> bool:
    fa_name = _GENOME_NAMES[assembly]
    dest = asm_dir / "genome" / fa_name
    fai = Path(str(dest) + ".fai")
    if dest.exists() and fai.exists():
        print(f"  [SKIP] {dest.name}")
        return True
    if skip:
        print("  [SKIP] --skip-genome specified")
        return True

    _require("samtools")
    dest.parent.mkdir(parents=True, exist_ok=True)

    # The .fai is derived locally from the FASTA, so a missing index never
    # implies a version mismatch — only re-acquire the (large) FASTA when the
    # FASTA itself is absent, then always (re)build the index below.
    if not dest.exists():
        existing = genome_fasta
        if existing is None:
            patterns = [f"*{assembly}*primary*assembly*.fa", f"*{assembly}*primary*assembly*.fasta",
                        "hg38.fa" if assembly == "GRCh38" else "hg19.fa"]
            existing = _find_existing(_GENOME_SEARCH, patterns)

        if existing and existing.exists():
            if str(existing).endswith(".gz"):
                print(f"  Decompressing: {existing.name} → {fa_name} (~3 GB)...")
                _run_shell(f"zcat {shlex.quote(str(existing))} > {shlex.quote(str(dest))}")
            else:
                print(f"  Creating symlink: {existing} → {dest}")
                dest.symlink_to(existing.resolve())
        else:
            gz = dest.with_suffix(".fa.gz")
            _download(urls["genome"], gz, f"Ensembl {assembly} primary assembly FASTA (~880 MB)")
            print("  Decompressing (~3 GB)...")
            _run_shell(f"zcat {shlex.quote(str(gz))} > {shlex.quote(str(dest))}")
            gz.unlink(missing_ok=True)

    print("  Indexing with samtools faidx...")
    _run(["samtools", "faidx", str(dest)])
    return True


def step_vep_cache(data_dir: Path, assembly: str, urls: dict, skip: bool) -> bool:
    cache_dir = data_dir / "vep_cache"
    asm_cache = cache_dir / "homo_sapiens_merged" / f"{ENSEMBL_RELEASE}_{assembly}"
    if asm_cache.exists() and any(asm_cache.iterdir()):
        print(f"  [SKIP] VEP cache exists ({asm_cache})")
        return True
    if skip:
        print("  [SKIP] --skip-vep-cache specified")
        return True

    tar = cache_dir / f"homo_sapiens_merged_vep_{ENSEMBL_RELEASE}_{assembly}.tar.gz"
    _download(urls["vep_cache"], tar, f"Ensembl VEP merged cache {assembly} (~14 GB)")
    print("  Extracting VEP cache...")
    cache_dir.mkdir(parents=True, exist_ok=True)
    _run(["tar", "-xzf", str(tar), "-C", str(cache_dir)])
    tar.unlink(missing_ok=True)
    return True


def step_clinvar_vcf(asm_dir: Path, assembly: str, urls: dict) -> bool:
    dest = asm_dir / "clinvar" / f"clinvar_{assembly}.vcf.gz"
    tbi = Path(str(dest) + ".tbi")
    if dest.exists() and tbi.exists():
        print(f"  [SKIP] {dest.name}")
        return True

    _require("tabix")
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Coupled set: the .vcf.gz and its .tbi MUST come from the same ClinVar
    # release. ClinVar is a rolling weekly release at a fixed URL, so a
    # separately-downloaded .tbi can index a newer release than the local
    # .vcf.gz, producing "Invalid BGZF header" errors on tabix random access.
    # If the set is incomplete, re-acquire the whole set: drop any orphaned
    # remnant, re-download the .vcf.gz, and (re)build the .tbi locally from the
    # on-disk bytes so the pair is guaranteed consistent.
    tbi.unlink(missing_ok=True)
    dest.unlink(missing_ok=True)
    _download(urls["clinvar_vcf"], dest, "ClinVar VCF (~120 MB)")
    print("  Building index with tabix...")
    _run(["tabix", "-p", "vcf", str(dest)])
    return True


def step_clinvar_sqlite(asm_dir: Path, assembly: str, urls: dict, workers: int | None = None) -> bool:
    dest = asm_dir / "clinvar" / f"clinvar_ps1_pm5_{assembly}.sqlite"
    if _verify_sqlite_has_rows(dest):
        print(f"  [SKIP] {dest.name}")
        return True

    xml_gz = asm_dir / "clinvar" / "ClinVarFullRelease.xml.gz"
    # ClinVarFullRelease is ~5 GB; treat anything under 1 GB as corrupt
    _verify_size(xml_gz, min_bytes=1_000_000_000, label="ClinVar XML")
    if not xml_gz.exists():
        _download(urls["clinvar_xml"], xml_gz, "ClinVar XML (~5 GB)")
        if not _verify_size(xml_gz, min_bytes=1_000_000_000, label="ClinVar XML"):
            raise RuntimeError("Failed to download ClinVar XML")

    _add_src_to_path()
    from acmg_classifier.setup.clinvar_builder import build_clinvar_sqlite  # type: ignore
    build_clinvar_sqlite(xml_gz, dest, assembly=assembly, workers=workers)
    return True


def step_alphamissense(asm_dir: Path, assembly: str, urls: dict) -> bool:
    suffix = "hg38" if assembly == "GRCh38" else "hg19"
    dest = asm_dir / "alphamissense" / f"AlphaMissense_{suffix}.tsv.gz"
    tbi = Path(str(dest) + ".tbi")
    if dest.exists() and tbi.exists():
        print(f"  [SKIP] {dest.name}")
        return True

    _require("tabix")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        _download(urls["alphamissense"], dest, "AlphaMissense (~600 MB)")
    if not tbi.exists():
        # AlphaMissense: col1=CHROM, col2=POS, header lines start with #
        _run(["tabix", "-s", "1", "-b", "2", "-e", "2", "-f", str(dest)])
    return True


def step_esm1b(data_dir: Path, urls: dict, skip: bool) -> bool:
    """Build ESM1b LLR SQLite from Brandes 2023 archive.

    Assembly-independent (protein-coordinate), so the SQLite lives at
    `data_dir/esm1b/esm1b_llr.sqlite` and is reused across GRCh37/GRCh38.
    """
    out_dir = data_dir / "esm1b"
    dest = out_dir / "esm1b_llr.sqlite"
    if _verify_sqlite_has_rows(dest, table="scores"):
        print(f"  [SKIP] {dest.name}")
        return True
    if skip:
        print("  [SKIP] --skip-esm1b specified")
        return True

    out_dir.mkdir(parents=True, exist_ok=True)
    zip_path = out_dir / "ALL_hum_isoforms_ESM1b_LLR.zip"
    _verify_size(zip_path, min_bytes=500_000_000, label="ESM1b zip")
    if not zip_path.exists():
        _download(urls["esm1b_zip"], zip_path, "ESM1b LLR archive (~1.34 GB)")

    print("  Building ESM1b SQLite...")
    _add_src_to_path()
    from acmg_classifier.setup.esm1b_builder import build_esm1b_sqlite  # type: ignore
    build_esm1b_sqlite(zip_path, dest)
    return True


def step_gnomad_constraint(asm_dir: Path, assembly: str, urls: dict) -> bool:
    ver = _GNOMAD_VER[assembly]
    dest = asm_dir / "gnomad" / f"gnomad_v{ver}_constraint.tsv"
    if dest.exists():
        print(f"  [SKIP] {dest.name}")
        return True

    dest.parent.mkdir(parents=True, exist_ok=True)
    url = urls["gnomad_constraint"]
    if url.endswith(".bgz"):
        # v2.1.1 is bgzip-compressed
        gz = dest.with_suffix(".tsv.bgz")
        _download(url, gz, f"gnomAD v{ver} constraint TSV (bgzipped)")
        _run_shell(f"zcat {shlex.quote(str(gz))} > {shlex.quote(str(dest))}")
        gz.unlink(missing_ok=True)
    else:
        _download(url, dest, f"gnomAD v{ver} constraint TSV (~3 MB)")
    return True


def step_gnomad_duckdb(
    asm_dir: Path,
    assembly: str,
    urls: dict,
    vcf_dir: Path | None,
    chromosomes: list[str],
    skip: bool,
    workers: int | None = None,
) -> bool:
    ver = _GNOMAD_VER[assembly]
    dest = asm_dir / "gnomad" / f"gnomad_v{ver}_exomes.duckdb"
    if _verify_duckdb_has_rows(dest):
        print(f"  [SKIP] {dest.name}")
        return True

    # Look for an existing VCF directory
    if vcf_dir is None:
        found = _find_existing(_GNOMAD_SEARCH, ["gnomad.exomes*.vcf.bgz"])
        if found:
            vcf_dir = found.parent

    staging = asm_dir / "gnomad" / "vcf"

    if vcf_dir is None or not list(vcf_dir.glob("*.vcf.bgz")):
        if skip:
            print(f"  [SKIP] Skipping gnomAD VCF download (--skip-gnomad)")
            print(f"  Once the VCFs are ready, re-run with --gnomad-vcf-dir <path>")
            return False
        total = len(chromosomes)
        print(f"  Downloading gnomAD v{ver} exomes ({total} chromosomes, ~300 GB total)")
        print("  This takes a while. Ctrl+C to interrupt; re-run to resume.")
        staging.mkdir(parents=True, exist_ok=True)
        for chrom in chromosomes:
            vcf_url = urls["gnomad_vcf"].format(chrom=chrom)
            tbi_url = urls["gnomad_vcf_tbi"].format(chrom=chrom)
            vcf_f = staging / f"gnomad.exomes.v{ver}.sites.{chrom}.vcf.bgz"
            tbi_f = Path(str(vcf_f) + ".tbi")
            # Re-download existing files that are suspiciously small (corrupt)
            _verify_size(vcf_f, min_bytes=10_000_000, label=f"gnomAD {chrom} VCF")
            _verify_size(tbi_f, min_bytes=1_000, label=f"gnomAD {chrom} .tbi")
            if not vcf_f.exists():
                _download(vcf_url, vcf_f, f"gnomAD {chrom}")
            if not tbi_f.exists():
                _download(tbi_url, tbi_f, f"gnomAD {chrom} .tbi")
        vcf_dir = staging

    n_workers = workers if workers is not None else max(1, (os.cpu_count() or 2) - 1)
    print(f"  Building gnomAD DuckDB (~25 GB, {n_workers} workers)...")
    _add_src_to_path()
    from acmg_classifier.setup.gnomad_builder import build_gnomad_duckdb  # type: ignore
    dest.parent.mkdir(parents=True, exist_ok=True)
    build_gnomad_duckdb(vcf_dir, dest, assembly=assembly, gnomad_version=ver,
                        max_workers=n_workers)
    return True


def step_repeatmasker(asm_dir: Path, assembly: str, urls: dict) -> bool:
    suffix = "hg38" if assembly == "GRCh38" else "hg19"
    dest = asm_dir / "repeats" / f"repeatmasker_dfam_{suffix}.bed.gz"
    tbi = Path(str(dest) + ".tbi")
    if dest.exists() and tbi.exists():
        print(f"  [SKIP] {dest.name}")
        return True

    _require("bgzip", "tabix", "awk", "sort")
    dest.parent.mkdir(parents=True, exist_ok=True)

    rmsk_raw = dest.parent / "rmsk.txt.gz"
    if not rmsk_raw.exists():
        _download(urls["repeatmasker"], rmsk_raw, f"UCSC RepeatMasker {suffix} (~120 MB)")

    if not dest.exists():
        print("  Converting to sorted BED + bgzip...")
        # rmsk.txt: col6=chrom, col7=start(0-based), col8=end, col12=class, col13=family
        cmd = (
            f"zcat {shlex.quote(str(rmsk_raw))} | "
            r"awk 'BEGIN{OFS=\"\t\"}{print $6, $7, $8, $12\"/\"$13}' | "
            f"sort -k1,1 -k2,2n | bgzip > {shlex.quote(str(dest))}"
        )
        subprocess.run(["bash", "-c", cmd], check=True)

    if not tbi.exists():
        _run(["tabix", "-p", "bed", str(dest)])
    return True


def step_mmsplice_gtf(asm_dir: Path, assembly: str, urls: dict, skip: bool) -> bool:
    """DISABLED — not registered in the steps list (MMSplice integration is off).
    Retained for re-enabling later. See the SpliceTool enum for context.

    Download the Ensembl GTF and filter to protein-coding genes for MMSplice.

    MMSplice (the open-source runtime splice predictor) needs a gene-annotation
    GTF. The upstream docs recommend filtering to protein-coding genes, which we
    do with a grep so the dataloader emits fewer, cleaner predictions.

    The `mmsplice` Python package itself (TensorFlow/Keras) is an OPTIONAL pip
    dependency and is intentionally NOT installed here — install it separately
    with `pip install -e .[mmsplice]` when you want to use --splice-tool mmsplice.
    """
    dest = asm_dir / "mmsplice" / GTF_NAMES[assembly]
    if dest.exists():
        print(f"  [SKIP] {dest.name}")
        return True
    if skip:
        print("  [SKIP] --skip-mmsplice-gtf specified")
        return True

    dest.parent.mkdir(parents=True, exist_ok=True)
    gz = dest.parent / "ensembl.gtf.gz"
    _verify_size(gz, min_bytes=10_000_000, label="Ensembl GTF")
    if not gz.exists():
        _download(urls["gtf"], gz, f"Ensembl GTF {assembly} (~50 MB)")

    print("  Filtering to protein-coding genes...")
    # Keep header lines (#) and any feature line whose attributes declare
    # gene_biotype "protein_coding".
    _run_shell(
        f"zcat {shlex.quote(str(gz))} | "
        f"grep -E '^#|gene_biotype \"protein_coding\"' > {shlex.quote(str(dest))}"
    )
    gz.unlink(missing_ok=True)
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Automated setup of data files for acmg-classifier",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--data-dir", type=Path, default=Path("./data"), metavar="PATH",
                        help="Root data directory (default: ./data)")
    parser.add_argument("--assembly", choices=["GRCh38", "GRCh37"], default="GRCh38")
    parser.add_argument("--genome-fasta", type=Path, default=None, metavar="PATH",
                        help="Path to an existing genome FASTA (skips download)")
    parser.add_argument("--gnomad-vcf-dir", type=Path, default=None, metavar="PATH",
                        help="Directory of existing gnomAD *.vcf.bgz files")
    parser.add_argument("--gnomad-chromosomes", nargs="+", default=None, metavar="CHR",
                        help="Chromosomes to download (default: all 24)")
    parser.add_argument("--gnomad-workers", type=int, default=None, metavar="N",
                        help="DuckDB build parallelism (default: CPU cores - 1)")
    parser.add_argument("--clinvar-workers", type=int, default=None, metavar="N",
                        help="ClinVar XML parse parallelism (default: 4, max: 24)")
    parser.add_argument("--skip-gnomad", action="store_true",
                        help="Skip gnomAD download (~300 GB)")
    parser.add_argument("--skip-genome", action="store_true",
                        help="Skip genome FASTA download (~880 MB)")
    parser.add_argument("--skip-vep-cache", action="store_true",
                        help="Skip VEP cache download (~14 GB)")
    parser.add_argument("--skip-esm1b", action="store_true",
                        help="Skip ESM1b download/build (~1.34 GB)")
    # MMSplice GTF DISABLED (MMSplice integration is off). Re-enable with:
    # parser.add_argument("--skip-mmsplice-gtf", action="store_true",
    #                     help="Skip MMSplice GTF download/filter (~50 MB)")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    assembly = args.assembly
    asm_dir = data_dir / assembly
    urls = URLS[assembly]
    chroms = args.gnomad_chromosomes or GNOMAD_CHROMS[assembly]

    sep = "=" * 62
    print(f"\n{sep}")
    print(f"  ACMG Classifier — data setup")
    print(f"  Data directory : {data_dir}")
    print(f"  Assembly       : {assembly}")
    print(f"{sep}\n")

    steps: list[tuple[str, object]] = [
        ("Reference genome",  lambda: step_genome(asm_dir, assembly, urls, args.genome_fasta, args.skip_genome)),
        ("VEP cache",         lambda: step_vep_cache(data_dir, assembly, urls, args.skip_vep_cache)),
        ("ClinVar VCF",       lambda: step_clinvar_vcf(asm_dir, assembly, urls)),
        ("ClinVar SQLite",    lambda: step_clinvar_sqlite(asm_dir, assembly, urls, args.clinvar_workers)),
        ("AlphaMissense",     lambda: step_alphamissense(asm_dir, assembly, urls)),
        ("ESM1b",             lambda: step_esm1b(data_dir, urls, args.skip_esm1b)),
        # MMSplice GTF DISABLED (MMSplice integration is off). Re-enable with:
        # ("MMSplice GTF",      lambda: step_mmsplice_gtf(asm_dir, assembly, urls, args.skip_mmsplice_gtf)),
        ("gnomAD constraint", lambda: step_gnomad_constraint(asm_dir, assembly, urls)),
        ("gnomAD DuckDB",     lambda: step_gnomad_duckdb(asm_dir, assembly, urls, args.gnomad_vcf_dir, chroms, args.skip_gnomad, args.gnomad_workers)),
        ("RepeatMasker",      lambda: step_repeatmasker(asm_dir, assembly, urls)),
    ]

    ok_steps, failed_steps = [], []
    for name, fn in steps:
        print(f"── {name} ──")
        try:
            result = fn()
            (ok_steps if result else failed_steps).append(name)
        except KeyboardInterrupt:
            print("\n[Interrupted]")
            sys.exit(130)
        except Exception as e:
            print(f"  [ERROR] {e}")
            failed_steps.append(name)
        print()

    print(sep)
    print(f"  Done: {len(ok_steps)} / {len(steps)} steps")
    if failed_steps:
        print(f"  Failed: {', '.join(failed_steps)}")
        sys.exit(1)
    else:
        print("  All files ready")
    print(sep + "\n")


if __name__ == "__main__":
    main()
