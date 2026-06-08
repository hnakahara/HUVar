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
from concurrent.futures import ThreadPoolExecutor, as_completed
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
        "clinvar_xml": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_00-latest.xml.gz",
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
        # JOINT (combined exome+genome) sites VCF — de-duplicated frequencies
        # (AF_joint / fafmax_faf95_max_joint) used for BA1/BS1/PM2/BS2.
        "gnomad_vcf": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/"
            "gnomad.joint.v4.1.sites.{chrom}.vcf.bgz"
        ),
        "gnomad_vcf_tbi": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/"
            "gnomad.joint.v4.1.sites.{chrom}.vcf.bgz.tbi"
        ),
        "repeatmasker": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz",
        "phylop": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw",
        # REVEL ships a single zip carrying BOTH hg19 and GRCh38 coordinates.
        "revel": "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip",
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
        "clinvar_xml": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_release/ClinVarRCVRelease_00-latest.xml.gz",
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
        # v2.1.1 has no joint release, so the genomes callset is loaded alongside
        # exomes and merged at query time (per-field MAX). Genomes have no chrY.
        "gnomad_vcf_genomes": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/"
            "gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz"
        ),
        "gnomad_vcf_genomes_tbi": (
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/"
            "gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz.tbi"
        ),
        "repeatmasker": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz",
        # hg19's bigWig is named hg19.100way.phyloP100way.bw (unlike hg38's
        # hg38.phyloP100way.bw); we still save it locally as hg19.phyloP100way.bw.
        "phylop": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw",
        # Same single REVEL zip as GRCh38; the hg19 position column is selected
        # when building the per-assembly TSV.
        "revel": "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip",
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
# gnomAD release "kind" per assembly (drives the DuckDB filename): v4.1 (GRCh38)
# uses the JOINT (combined exome+genome) sites VCF; v2.1.1 (GRCh37) has no joint
# release, so exomes AND genomes are both loaded and merged at query time.
_GNOMAD_KIND = {"GRCh38": "joint", "GRCh37": "exome_genome"}

# gnomAD is mirrored byte-identically on Google Cloud and AWS: the object key
# after the bucket prefix is the same on both, so we can swap only the prefix.
# Downloading from both mirrors in parallel roughly doubles throughput.
# Verified 2026-06: both return HTTP 200 for v4.1 (joint) and v2.1.1
# (exomes/genomes). The URLS dict above stores the Google form.
_GNOMAD_GCS_PREFIX = "https://storage.googleapis.com/gcp-public-data--gnomad/"
_GNOMAD_AWS_PREFIX = "https://gnomad-public-us-east-1.s3.amazonaws.com/"
# (label, base URL) pairs, in round-robin assignment order.
_GNOMAD_MIRRORS: list[tuple[str, str]] = [
    ("google", _GNOMAD_GCS_PREFIX),
    ("amazon", _GNOMAD_AWS_PREFIX),
]
# Concurrent downloads per mirror (1 google + 1 amazon = 2 total in flight).
_GNOMAD_PER_MIRROR = 1

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


def _remote_size(url: str) -> int | None:
    """Return the remote file's Content-Length in bytes, or None if unknown.

    Used to tell a fully-downloaded file from a partial one left by an
    interrupted run: a partial VCF can be larger than the corruption
    threshold yet still be incomplete, so existence alone is not enough.
    """
    try:
        if _has("curl"):
            out = subprocess.run(
                ["curl", "-sIL", "-m", "30", url],
                capture_output=True, text=True, check=True,
            ).stdout
        elif _has("wget"):
            # --spider does a HEAD; server response goes to stderr.
            out = subprocess.run(
                ["wget", "--spider", "--server-response", "-T", "30", "-t", "1", url],
                capture_output=True, text=True,
            ).stderr
        else:
            return None
    except subprocess.SubprocessError:
        return None
    # With redirects there can be several Content-Length lines; the final one
    # corresponds to the actual object (200 response).
    size: int | None = None
    for line in out.splitlines():
        if line.lower().strip().startswith("content-length:"):
            try:
                size = int(line.split(":", 1)[1].strip())
            except ValueError:
                pass
    return size


def _download_quiet(url: str, dest: Path, label: str = "") -> None:
    """Download without an interactive progress bar.

    Used by the parallel gnomAD downloader where concurrent wget/curl progress
    bars would otherwise interleave into unreadable output. Keeps the resume
    flag (wget -c / curl -C -) so an interrupted transfer continues from its
    current offset instead of restarting.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if _has("wget"):
        cmd = ["wget", "-c", "--no-verbose", "-O", str(dest), url]
    elif _has("curl"):
        cmd = ["curl", "-C", "-", "-L", "-s", "-S", "-o", str(dest), url]
    else:
        sys.exit("[ERROR] wget or curl is required")
    subprocess.run([str(x) for x in cmd], check=True)


def _download_gnomad_chrom(
    label: str,
    base: str,
    name: str,
    chrom: str,
    vcf_suffix: str,
    tbi_suffix: str,
    vcf_f: Path,
    tbi_f: Path,
) -> str:
    """Download one chromosome's VCF + .tbi from the assigned mirror.

    On failure, falls back to the other mirror (the gnomAD GCS/AWS mirrors are
    byte-identical) before giving up. Runs in a worker thread; returns a short
    status string for logging. Already-present, non-corrupt files are skipped.
    """
    def _fetch(suffix: str, dest: Path, what: str) -> None:
        # Try the assigned mirror first, then the other mirror as a fallback.
        attempts = [(label, base)] + [m for m in _GNOMAD_MIRRORS if m[0] != label]
        last_exc: Exception | None = None
        for mlabel, mbase in attempts:
            url = mbase + suffix
            remote = _remote_size(url)
            local = dest.stat().st_size if dest.exists() else 0
            # Skip ONLY when the local file matches the remote size exactly. A
            # partial file from an interrupted run (local < remote) falls
            # through to wget -c, which resumes from the current offset rather
            # than being mistaken for complete and skipped.
            if dest.exists() and remote is not None and local == remote:
                return
            try:
                action = "resuming" if local > 0 else "downloading"
                print(f"  ↓  gnomAD {name} {chrom} {what}  [{mlabel}] ({action})")
                _download_quiet(url, dest, f"gnomAD {name} {chrom} {what}")
                # Confirm completeness against the remote size when known.
                if remote is not None and dest.exists() and dest.stat().st_size != remote:
                    raise RuntimeError(
                        f"size mismatch after download "
                        f"({dest.stat().st_size} != {remote})"
                    )
                return
            except (subprocess.CalledProcessError, RuntimeError) as exc:
                last_exc = exc
                print(f"  [WARN] {mlabel} failed for {name} {chrom} {what}; trying next mirror")
        raise RuntimeError(f"All mirrors failed for gnomAD {name} {chrom} {what}") from last_exc

    _fetch(vcf_suffix, vcf_f, "VCF")
    _fetch(tbi_suffix, tbi_f, ".tbi")
    return f"{name} {chrom} [{label}]"


def _download_gnomad_parallel(jobs: list[tuple], per_mirror: int = _GNOMAD_PER_MIRROR) -> None:
    """Download all gnomAD chromosome jobs across both mirrors in parallel.

    Each job is one chromosome (VCF + .tbi). Jobs are assigned round-robin to
    the mirrors, and each mirror gets its own thread pool of `per_mirror`
    workers — so at steady state there are `per_mirror` downloads in flight
    from each mirror simultaneously.
    """
    buckets: dict[str, list[tuple]] = {label: [] for label, _ in _GNOMAD_MIRRORS}
    for i, job in enumerate(jobs):
        label = _GNOMAD_MIRRORS[i % len(_GNOMAD_MIRRORS)][0]
        buckets[label].append(job)

    base_by_label = dict(_GNOMAD_MIRRORS)
    pools: list[ThreadPoolExecutor] = []
    futures = []
    try:
        for label, _ in _GNOMAD_MIRRORS:
            if not buckets[label]:
                continue
            pool = ThreadPoolExecutor(max_workers=per_mirror, thread_name_prefix=f"dl-{label}")
            pools.append(pool)
            base = base_by_label[label]
            for job in buckets[label]:
                futures.append(pool.submit(_download_gnomad_chrom, label, base, *job))

        errors: list[Exception] = []
        for fut in as_completed(futures):
            try:
                done = fut.result()
                print(f"  ✓  {done}")
            except Exception as exc:  # noqa: BLE001 — aggregate and report below
                errors.append(exc)
                print(f"  [ERROR] {exc}")
        if errors:
            raise errors[0]
    finally:
        for pool in pools:
            pool.shutdown()


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

    # ClinVarRCVRelease (RCV_release) — the current, weekly/monthly-updated RCV
    # XML. The legacy RCV_xml_old_format/ClinVarFullRelease_00-latest is frozen
    # (its content stopped at 2025-07), so PM5/PS1 comparators built from it lag
    # ClinVar by ~a year. The distinct local name ensures a stale, previously
    # downloaded ClinVarFullRelease.xml.gz is not silently reused.
    xml_gz = asm_dir / "clinvar" / "ClinVarRCVRelease.xml.gz"
    # ClinVarRCVRelease is ~6 GB; treat anything under 1 GB as corrupt
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


def step_revel(asm_dir: Path, assembly: str, urls: dict, enabled: bool) -> bool:
    """Build the per-assembly REVEL TSV for PP3/BP4 (--insilico-tool revel).

    REVEL distributes one zip (`revel_with_transcript_ids`, a CSV carrying both
    hg19 and GRCh38 coordinates). We extract a clean, tab-separated, position-
    sorted, bgzipped 5-column file (chrom, pos, ref, alt, REVEL) indexed on the
    assembly's coordinate column, so the query side stays assembly-agnostic.

    OPT-IN (--with-revel): the source zip is ~600 MB and ESM1b is the default
    in-silico tool, so REVEL is only fetched when explicitly requested."""
    suffix = "grch38" if assembly == "GRCh38" else "grch37"
    dest = asm_dir / "revel" / f"revel_{suffix}.tsv.gz"
    tbi = Path(str(dest) + ".tbi")
    if dest.exists() and tbi.exists():
        print(f"  [SKIP] {dest.name}")
        return True
    if not enabled:
        print("  [SKIP] REVEL not requested (pass --with-revel to enable "
              "--insilico-tool revel; ~600 MB zip)")
        return True

    _require("unzip", "awk", "sort", "bgzip", "tabix")
    dest.parent.mkdir(parents=True, exist_ok=True)

    zip_path = dest.parent / "revel-v1.3_all_chromosomes.zip"
    if not zip_path.exists():
        _download(urls["revel"], zip_path, "REVEL v1.3 (~600 MB)")

    if not dest.exists():
        # CSV columns: chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,txids.
        # GRCh38 reads the grch38_pos column (3) and drops un-lifted rows ('.');
        # GRCh37 reads hg19_pos (2). REVEL is column 8.
        pos_col = "$3" if assembly == "GRCh38" else "$2"
        guard = '$3!="" && $3!="."' if assembly == "GRCh38" else '$2!="" && $2!="."'
        print("  Converting REVEL CSV to sorted TSV + bgzip...")
        awk = (
            f"BEGIN{{OFS=\"\\t\"}} NR>1 && {guard} "
            f"{{print $1, {pos_col}, $4, $5, $8}}"
        )
        cmd = (
            f"unzip -p {shlex.quote(str(zip_path))} | "
            f"awk -F, {shlex.quote(awk)} | "
            f"sort -k1,1 -k2,2n | bgzip > {shlex.quote(str(dest))}"
        )
        subprocess.run(["bash", "-c", cmd], check=True)

    if not tbi.exists():
        _run(["tabix", "-s", "1", "-b", "2", "-e", "2", "-f", str(dest)])
    return True


# OpenSpliceAI OSAI_MANE pretrained models. Each flanking size has a 5-model
# ensemble (random seeds rs10–rs14), mirroring SpliceAI's 5-model averaging.
# Hosted on the JHU CCB FTP; the openspliceai CLI takes the per-flank directory
# via -m and averages every .pt it finds.
_OSAI_FLANKS = (80, 400, 2000, 10000)
_OSAI_SEEDS = (10, 11, 12, 13, 14)
_OSAI_FTP = (
    "ftp://ftp.ccb.jhu.edu/pub/data/OpenSpliceAI/OSAI-MANE/"
    "{flank}nt/model_{flank}nt_rs{seed}.pt"
)
# Gene annotation tables (#NAME/CHROM/STRAND/TX_START/TX_END/EXON_START/EXON_END),
# passed to `openspliceai variant -A`. Distributed in the OpenSpliceAI repo.
_OSAI_ANNOTATION = (
    "https://raw.githubusercontent.com/Kuanhao-Chao/OpenSpliceAI/main/data/{name}"
)


def step_openspliceai(asm_dir: Path, assembly: str, skip: bool) -> bool:
    """Download the OSAI_MANE model ensembles AND the gene annotation table.

    The default --splice-tool is openspliceai, so this runs by default (opt out
    with --skip-openspliceai). It downloads the 5-model ensemble for ALL four
    flanking sizes (80/400/2000/10000 nt) into data/<asm>/openspliceai/<flank>nt/,
    exactly the layout OpenSpliceAIPredictor / Config.openspliceai_model_path
    expect. The model weights are sequence-based (assembly-independent).

    It ALSO downloads the gene annotation table (grch38.txt / grch37.txt). The
    `openspliceai variant -A grch38` keyword is NOT usable from an installed
    package — openspliceai maps it to the relative path "./data/vcf/<asm>.txt"
    which only resolves inside the openspliceai source tree — so we stage the
    table explicitly and pass its path via Config.openspliceai_annotation.

    The `openspliceai` CLI itself is a regular project dependency (pyproject),
    installed with the package; this step stages the runtime data files."""
    if skip:
        print("  [SKIP] --skip-openspliceai specified")
        return True

    base = asm_dir / "openspliceai"
    n_have = 0
    for flank in _OSAI_FLANKS:
        d = base / f"{flank}nt"
        d.mkdir(parents=True, exist_ok=True)
        for seed in _OSAI_SEEDS:
            dest = d / f"model_{flank}nt_rs{seed}.pt"
            if dest.exists():
                n_have += 1
                continue
            _download(_OSAI_FTP.format(flank=flank, seed=seed), dest,
                      f"OSAI_MANE {flank}nt rs{seed}")
            n_have += 1
    print(f"  OSAI_MANE models ready: {n_have}/{len(_OSAI_FLANKS) * len(_OSAI_SEEDS)} "
          f"across {len(_OSAI_FLANKS)} flanking sizes")

    # Gene annotation table (-A): explicit file, required because the keyword is
    # broken from an installed package (see docstring).
    ann_name = "grch38.txt" if assembly == "GRCh38" else "grch37.txt"
    ann_dest = base / ann_name
    if ann_dest.exists():
        print(f"  [SKIP] {ann_name}")
    else:
        base.mkdir(parents=True, exist_ok=True)
        _download(_OSAI_ANNOTATION.format(name=ann_name), ann_dest,
                  f"OpenSpliceAI annotation table {ann_name}")
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
    kind = _GNOMAD_KIND[assembly]
    dest = asm_dir / "gnomad" / f"gnomad_v{ver}_{kind}.duckdb"
    if _verify_duckdb_has_rows(dest):
        print(f"  [SKIP] {dest.name}")
        return True

    # gnomAD callsets to load. GRCh38: the single joint VCF. GRCh37: exomes AND
    # genomes (no joint release; merged at query time). genomes have no chrY.
    primary = "joint" if assembly == "GRCh38" else "exomes"
    sources = [(primary, urls["gnomad_vcf"], urls["gnomad_vcf_tbi"], chromosomes)]
    if "gnomad_vcf_genomes" in urls:
        genome_chroms = [c for c in chromosomes if c not in ("chrY", "Y", "chrM", "MT")]
        sources.append((
            "genomes", urls["gnomad_vcf_genomes"], urls["gnomad_vcf_genomes_tbi"],
            genome_chroms,
        ))

    # Look for an existing VCF directory (any gnomAD callset VCFs present).
    if vcf_dir is None:
        found = _find_existing(_GNOMAD_SEARCH, [f"gnomad.{name}*.vcf.bgz" for name, *_ in sources])
        if found:
            vcf_dir = found.parent

    staging = asm_dir / "gnomad" / "vcf"

    if vcf_dir is None or not list(vcf_dir.glob("*.vcf.bgz")):
        if skip:
            print(f"  [SKIP] Skipping gnomAD VCF download (--skip-gnomad)")
            print(f"  Once the VCFs are ready, re-run with --gnomad-vcf-dir <path>")
            return False
        print(f"  Downloading gnomAD v{ver} ({', '.join(n for n, *_ in sources)}; "
              f"~300+ GB total) from {len(_GNOMAD_MIRRORS)} mirrors, "
              f"{_GNOMAD_PER_MIRROR} concurrent each. Ctrl+C to interrupt; re-run to resume.")
        staging.mkdir(parents=True, exist_ok=True)

        # Build one job per chromosome (VCF + .tbi). URL suffixes are derived by
        # stripping the Google prefix so each job can be fetched from either
        # mirror; _download_gnomad_parallel splits them across both.
        jobs: list[tuple] = []
        for name, vcf_tmpl, tbi_tmpl, src_chroms in sources:
            for chrom in src_chroms:
                vcf_suffix = vcf_tmpl.format(chrom=chrom).removeprefix(_GNOMAD_GCS_PREFIX)
                tbi_suffix = tbi_tmpl.format(chrom=chrom).removeprefix(_GNOMAD_GCS_PREFIX)
                vcf_f = staging / f"gnomad.{name}.v{ver}.sites.{chrom}.vcf.bgz"
                tbi_f = Path(str(vcf_f) + ".tbi")
                jobs.append((name, chrom, vcf_suffix, tbi_suffix, vcf_f, tbi_f))

        _download_gnomad_parallel(jobs)
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


def step_phylop(asm_dir: Path, assembly: str, urls: dict, enabled: bool) -> bool:
    """Download the phyloP100way conservation bigWig for the BP7 conservation
    gate. OPT-IN (--with-phylop) because the track is ~9.2 GB. Skipped by
    default; BP7 then falls back to its splice-only logic."""
    suffix = "hg38" if assembly == "GRCh38" else "hg19"
    dest = asm_dir / "conservation" / f"{suffix}.phyloP100way.bw"
    if not enabled and not dest.exists():
        print("  [SKIP] phyloP not requested (pass --with-phylop to enable the "
              "BP7 conservation gate; ~9.2 GB)")
        return True
    # Skip ONLY when the local file matches the remote size exactly. A partial
    # file from an interrupted run (local < remote) must NOT be treated as
    # complete — it falls through to wget -c, which resumes from the current
    # offset rather than being mistakenly skipped.
    if dest.exists():
        remote = _remote_size(urls["phylop"])
        local = dest.stat().st_size
        if remote is not None and local == remote:
            print(f"  [SKIP] {dest.name}")
            return True
        if remote is None:
            # Cannot verify size; assume an existing file is complete to avoid
            # re-downloading ~9 GB on every run when the server omits a HEAD size.
            print(f"  [SKIP] {dest.name} (remote size unknown; assuming complete)")
            return True
        print(f"  ↻  {dest.name} is partial ({local}/{remote} bytes) — resuming")
    dest.parent.mkdir(parents=True, exist_ok=True)
    _download(urls["phylop"], dest, f"UCSC phyloP100way {suffix} bigWig (~9.2 GB)")
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
    parser.add_argument("--workers", type=int, default=None, metavar="N",
                        help="Build parallelism for the ClinVar (XML parse, max 24) "
                             "and gnomAD (DuckDB) steps (default: CPU cores - 1)")
    parser.add_argument("--skip-gnomad", action="store_true",
                        help="Skip gnomAD download (~300 GB)")
    parser.add_argument("--skip-genome", action="store_true",
                        help="Skip genome FASTA download (~880 MB)")
    parser.add_argument("--skip-vep-cache", action="store_true",
                        help="Skip VEP cache download (~14 GB)")
    parser.add_argument("--with-phylop", action="store_true",
                        help="Download the phyloP100way bigWig (~9.2 GB) for the BP7 "
                             "conservation gate. Off by default; BP7 uses splice-only "
                             "logic without it. The pyBigWig reader is installed by "
                             "default, so the gate activates once this file is present.")
    parser.add_argument("--with-revel", action="store_true",
                        help="Download REVEL (~600 MB zip) and build the per-assembly "
                             "TSV for --insilico-tool revel. Off by default (ESM1b is "
                             "the default in-silico tool).")
    parser.add_argument("--skip-esm1b", action="store_true",
                        help="Skip ESM1b download/build (~1.34 GB)")
    parser.add_argument("--skip-openspliceai", action="store_true",
                        help="Skip the OSAI_MANE model download (all 4 flanking sizes). "
                             "The default --splice-tool is openspliceai, so skip only if "
                             "supplying models another way (e.g. --openspliceai-model-dir). "
                             "The openspliceai CLI is a package dependency, installed with "
                             "the tool regardless of this flag.")
    # MMSplice GTF DISABLED (MMSplice integration is off). Re-enable with:
    # parser.add_argument("--skip-mmsplice-gtf", action="store_true",
    #                     help="Skip MMSplice GTF download/filter (~50 MB)")
    parser.add_argument("--only", nargs="+", default=None, metavar="STEP",
                        help="Run only the named step(s) and skip the rest. Step "
                             "keys: genome, vep, clinvar-vcf, clinvar-sqlite, "
                             "alphamissense, esm1b, revel, openspliceai, "
                             "gnomad-constraint, gnomad, repeatmasker, phylop. "
                             "Each step is idempotent: it checks for existing "
                             "files and downloads only what is missing. Example: "
                             "--only openspliceai")
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

    steps: list[tuple[str, str, object]] = [
        ("genome",            "Reference genome",  lambda: step_genome(asm_dir, assembly, urls, args.genome_fasta, args.skip_genome)),
        ("vep",               "VEP cache",         lambda: step_vep_cache(data_dir, assembly, urls, args.skip_vep_cache)),
        ("clinvar-vcf",       "ClinVar VCF",       lambda: step_clinvar_vcf(asm_dir, assembly, urls)),
        ("clinvar-sqlite",    "ClinVar SQLite",    lambda: step_clinvar_sqlite(asm_dir, assembly, urls, args.workers)),
        ("alphamissense",     "AlphaMissense",     lambda: step_alphamissense(asm_dir, assembly, urls)),
        ("esm1b",             "ESM1b",             lambda: step_esm1b(data_dir, urls, args.skip_esm1b)),
        ("revel",             "REVEL",             lambda: step_revel(asm_dir, assembly, urls, args.with_revel)),
        ("openspliceai",      "OpenSpliceAI",      lambda: step_openspliceai(asm_dir, assembly, args.skip_openspliceai)),
        # MMSplice GTF DISABLED (MMSplice integration is off). Re-enable with:
        # ("mmsplice-gtf",      "MMSplice GTF",      lambda: step_mmsplice_gtf(asm_dir, assembly, urls, args.skip_mmsplice_gtf)),
        ("gnomad-constraint", "gnomAD constraint", lambda: step_gnomad_constraint(asm_dir, assembly, urls)),
        ("gnomad",            "gnomAD DuckDB",     lambda: step_gnomad_duckdb(asm_dir, assembly, urls, args.gnomad_vcf_dir, chroms, args.skip_gnomad, args.workers)),
        ("repeatmasker",      "RepeatMasker",      lambda: step_repeatmasker(asm_dir, assembly, urls)),
        ("phylop",            "phyloP (BP7)",      lambda: step_phylop(asm_dir, assembly, urls, args.with_phylop)),
    ]

    if args.only:
        valid = {key for key, _, _ in steps}
        requested = set(args.only)
        unknown = requested - valid
        if unknown:
            sys.exit(f"[ERROR] Unknown --only step(s): {', '.join(sorted(unknown))}. "
                     f"Valid keys: {', '.join(k for k, _, _ in steps)}")
        steps = [s for s in steps if s[0] in requested]

    ok_steps, failed_steps = [], []
    for _key, name, fn in steps:
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
