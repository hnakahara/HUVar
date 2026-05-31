#!/usr/bin/env python3
"""
scripts/setup_data.py

acmg-classifier 用データファイルの自動セットアップ。
既に存在するファイルはスキップし、不足しているものだけダウンロード・ビルドする。

使い方:
  python scripts/setup_data.py --data-dir ./data

  # 既存ファイルを指定してダウンロードをスキップ
  python scripts/setup_data.py --data-dir ./data \\
      --genome-fasta /db/reference/GRCh38/hg38.fa \\
      --gnomad-vcf-dir /db/gnomad/v4.1/exomes/vcf

  # gnomAD (~300 GB) はスキップして他だけセットアップ
  python scripts/setup_data.py --data-dir ./data --skip-gnomad

  # 特定染色体のみ gnomAD をダウンロード
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
# 定数
# ---------------------------------------------------------------------------

ENSEMBL_RELEASE = 111

SQUIRLS_VERSION = "2309"

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
        "clinvar_vcf_tbi": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi",
        "clinvar_xml": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_xml_old_format/ClinVarFullRelease_00-latest.xml.gz",
        "alphamissense": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz",
        "esm1b_zip": (
            "https://huggingface.co/spaces/ntranoslab/esm_variants/resolve/main/"
            "ALL_hum_isoforms_ESM1b_LLR.zip"
        ),
        "squirls_zip": (
            f"https://squirls.s3.amazonaws.com/squirls-{SQUIRLS_VERSION}-hg38.zip"
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
        "clinvar_vcf_tbi": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi",
        "clinvar_xml": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/RCV_xml_old_format/ClinVarFullRelease_00-latest.xml.gz",
        "alphamissense": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz",
        "esm1b_zip": (
            "https://huggingface.co/spaces/ntranoslab/esm_variants/resolve/main/"
            "ALL_hum_isoforms_ESM1b_LLR.zip"
        ),
        "squirls_zip": (
            f"https://squirls.s3.amazonaws.com/squirls-{SQUIRLS_VERSION}-hg19.zip"
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

# GRCh38(v4.1): chr prefix あり / GRCh37(v2.1.1): なし
GNOMAD_CHROMS: dict[str, list[str]] = {
    "GRCh38": [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]],
    "GRCh37": [str(c) for c in list(range(1, 23)) + ["X", "Y"]],
}

# サーバー上の既存 gnomAD / genome の候補ディレクトリ
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
# ユーティリティ
# ---------------------------------------------------------------------------


def _has(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _require(*cmds: str) -> None:
    missing = [c for c in cmds if not _has(c)]
    if missing:
        sys.exit(
            f"[ERROR] 必要なツールが見つかりません: {', '.join(missing)}\n"
            f"  インストール例: conda install -c bioconda {' '.join(missing)}"
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
        sys.exit("[ERROR] wget または curl が必要です")


def _verify_size(path: Path, min_bytes: int, label: str = "") -> bool:
    """ファイルサイズが期待値より極端に小さければ破損と判定して削除する。"""
    if not path.exists():
        return False
    sz = path.stat().st_size
    if sz < min_bytes:
        print(f"  [WARN] {label or path.name} がサイズ不足 ({sz:,} bytes < {min_bytes:,})")
        print(f"         破損または不完全なため削除して再取得します")
        path.unlink(missing_ok=True)
        return False
    return True


def _verify_sqlite_has_rows(path: Path, table: str = "variants") -> bool:
    """SQLite ファイルにテーブルが存在しレコードがあれば True。空なら削除して False。"""
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
        print(f"  [WARN] {path.name} は空 (0 rows in {table}) — 削除して再作成します")
    except Exception as e:
        print(f"  [WARN] {path.name} は壊れています ({e}) — 削除して再作成します")
    path.unlink(missing_ok=True)
    return False


def _verify_duckdb_has_rows(path: Path, table: str = "variants") -> bool:
    """DuckDB ファイルにテーブルが存在しレコードがあれば True。空なら削除して False。"""
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
        print(f"  [WARN] {path.name} は空 (0 rows in {table}) — 削除して再作成します")
    except Exception as e:
        print(f"  [WARN] {path.name} は壊れています ({e}) — 削除して再作成します")
    path.unlink(missing_ok=True)
    return False


def _find_existing(search_dirs: list[str], patterns: list[str]) -> Path | None:
    for d in search_dirs:
        p = Path(d)
        if not p.is_dir():
            continue
        for pat in patterns:
            for f in p.rglob(pat):
                print(f"  既存ファイルを発見: {f}")
                return f
    return None


def _add_src_to_path() -> None:
    src = Path(__file__).parent.parent / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


# ---------------------------------------------------------------------------
# セットアップ ステップ
# ---------------------------------------------------------------------------


def step_genome(asm_dir: Path, assembly: str, urls: dict, genome_fasta: Path | None, skip: bool) -> bool:
    fa_name = _GENOME_NAMES[assembly]
    dest = asm_dir / "genome" / fa_name
    if dest.exists():
        print(f"  [SKIP] {dest.name}")
        return True
    if skip:
        print("  [SKIP] --skip-genome 指定")
        return True

    _require("samtools")
    dest.parent.mkdir(parents=True, exist_ok=True)

    existing = genome_fasta
    if existing is None:
        patterns = [f"*{assembly}*primary*assembly*.fa", f"*{assembly}*primary*assembly*.fasta",
                    "hg38.fa" if assembly == "GRCh38" else "hg19.fa"]
        existing = _find_existing(_GENOME_SEARCH, patterns)

    if existing and existing.exists():
        if str(existing).endswith(".gz"):
            print(f"  展開中: {existing.name} → {fa_name} (~3 GB)...")
            _run_shell(f"zcat {shlex.quote(str(existing))} > {shlex.quote(str(dest))}")
        else:
            print(f"  シンボリックリンク作成: {existing} → {dest}")
            dest.symlink_to(existing.resolve())
    else:
        gz = dest.with_suffix(".fa.gz")
        _download(urls["genome"], gz, f"Ensembl {assembly} primary assembly FASTA (~880 MB)")
        print("  展開中 (~3 GB)...")
        _run_shell(f"zcat {shlex.quote(str(gz))} > {shlex.quote(str(dest))}")
        gz.unlink(missing_ok=True)

    print("  samtools faidx でインデックス作成中...")
    _run(["samtools", "faidx", str(dest)])
    return True


def step_vep_cache(data_dir: Path, assembly: str, urls: dict, skip: bool) -> bool:
    cache_dir = data_dir / "vep_cache"
    asm_cache = cache_dir / "homo_sapiens_merged" / f"{ENSEMBL_RELEASE}_{assembly}"
    if asm_cache.exists() and any(asm_cache.iterdir()):
        print(f"  [SKIP] VEP キャッシュ既存 ({asm_cache})")
        return True
    if skip:
        print("  [SKIP] --skip-vep-cache 指定")
        return True

    tar = cache_dir / f"homo_sapiens_merged_vep_{ENSEMBL_RELEASE}_{assembly}.tar.gz"
    _download(urls["vep_cache"], tar, f"Ensembl VEP merged cache {assembly} (~14 GB)")
    print("  VEP キャッシュ展開中...")
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
    if not dest.exists():
        _download(urls["clinvar_vcf"], dest, "ClinVar VCF (~120 MB)")
    if not tbi.exists():
        try:
            _download(urls["clinvar_vcf_tbi"], tbi, "ClinVar VCF .tbi")
        except subprocess.CalledProcessError:
            print("  tabix でインデックス再作成中...")
            _run(["tabix", "-p", "vcf", str(dest)])
    return True


def step_clinvar_sqlite(asm_dir: Path, assembly: str, urls: dict) -> bool:
    dest = asm_dir / "clinvar" / f"clinvar_ps1_pm5_{assembly}.sqlite"
    if _verify_sqlite_has_rows(dest):
        print(f"  [SKIP] {dest.name}")
        return True

    xml_gz = asm_dir / "clinvar" / "ClinVarFullRelease.xml.gz"
    # ClinVarFullRelease は ~5 GB。1 GB 未満なら破損とみなす
    _verify_size(xml_gz, min_bytes=1_000_000_000, label="ClinVar XML")
    if not xml_gz.exists():
        _download(urls["clinvar_xml"], xml_gz, "ClinVar XML (~5 GB)")
        if not _verify_size(xml_gz, min_bytes=1_000_000_000, label="ClinVar XML"):
            raise RuntimeError("ClinVar XML のダウンロードに失敗しました")

    _add_src_to_path()
    from acmg_classifier.setup.clinvar_builder import build_clinvar_sqlite  # type: ignore
    build_clinvar_sqlite(xml_gz, dest, assembly=assembly)
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
        # AlphaMissense: col1=CHROM, col2=POS, ヘッダーは # 始まり
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
        print("  [SKIP] --skip-esm1b 指定")
        return True

    out_dir.mkdir(parents=True, exist_ok=True)
    zip_path = out_dir / "ALL_hum_isoforms_ESM1b_LLR.zip"
    _verify_size(zip_path, min_bytes=500_000_000, label="ESM1b zip")
    if not zip_path.exists():
        _download(urls["esm1b_zip"], zip_path, "ESM1b LLR archive (~1.34 GB)")

    print("  ESM1b SQLite 構築中...")
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
        # v2.1.1 は bgzip 圧縮
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

    # 既存 VCF ディレクトリを探す
    if vcf_dir is None:
        found = _find_existing(_GNOMAD_SEARCH, ["gnomad.exomes*.vcf.bgz"])
        if found:
            vcf_dir = found.parent

    staging = asm_dir / "gnomad" / "vcf"

    if vcf_dir is None or not list(vcf_dir.glob("*.vcf.bgz")):
        if skip:
            print(f"  [SKIP] gnomAD VCF ダウンロードをスキップ (--skip-gnomad)")
            print(f"  VCF が準備できたら --gnomad-vcf-dir <path> で再実行してください")
            return False
        total = len(chromosomes)
        print(f"  gnomAD v{ver} exomes ダウンロード ({total} 染色体, 合計 ~300 GB)")
        print("  時間がかかります。Ctrl+C で中断、再実行で再開できます。")
        staging.mkdir(parents=True, exist_ok=True)
        for chrom in chromosomes:
            vcf_url = urls["gnomad_vcf"].format(chrom=chrom)
            tbi_url = urls["gnomad_vcf_tbi"].format(chrom=chrom)
            vcf_f = staging / f"gnomad.exomes.v{ver}.sites.{chrom}.vcf.bgz"
            tbi_f = Path(str(vcf_f) + ".tbi")
            # 既存ファイルが極端に小さければ破損とみなして再取得
            _verify_size(vcf_f, min_bytes=10_000_000, label=f"gnomAD {chrom} VCF")
            _verify_size(tbi_f, min_bytes=1_000, label=f"gnomAD {chrom} .tbi")
            if not vcf_f.exists():
                _download(vcf_url, vcf_f, f"gnomAD {chrom}")
            if not tbi_f.exists():
                _download(tbi_url, tbi_f, f"gnomAD {chrom} .tbi")
        vcf_dir = staging

    n_workers = workers if workers is not None else max(1, (os.cpu_count() or 2) - 1)
    print(f"  gnomAD DuckDB 構築中 (~25 GB, {n_workers} 並列)...")
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
        print("  sorted BED に変換 + bgzip...")
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


def step_squirls(asm_dir: Path, assembly: str, urls: dict, squirls_db: Path | None, skip: bool) -> bool:
    """Download and place the SQUIRLS precomputed splice-score database.

    SQUIRLS distributes a precomputed SQLite DB (~4 GB) as a zip archive.
    The zip contains a single .db file which is extracted to:
      <asm_dir>/squirls/squirls-<version>-<suffix>/

    If --squirls-db points to an existing .db file, it is used directly
    (no download). Use --skip-squirls when using SpliceAI instead.
    """
    suffix = "hg38" if assembly == "GRCh38" else "hg19"
    dest_dir = asm_dir / "squirls" / f"squirls-{SQUIRLS_VERSION}-{suffix}"

    # 既存 .db を確認
    existing_db = next(dest_dir.glob("*.db"), None) if dest_dir.exists() else None
    if existing_db:
        print(f"  [SKIP] {existing_db.name}")
        return True
    if skip:
        print("  [SKIP] --skip-squirls 指定")
        return True

    dest_dir.mkdir(parents=True, exist_ok=True)

    # --squirls-db で既存ファイルを指定された場合はコピー/シンボリックリンク
    if squirls_db and squirls_db.exists():
        dest = dest_dir / squirls_db.name
        print(f"  シンボリックリンク作成: {squirls_db} → {dest}")
        dest.symlink_to(squirls_db.resolve())
        return True

    # zip をダウンロードして展開
    import zipfile
    zip_path = dest_dir / f"squirls-{SQUIRLS_VERSION}-{suffix}.zip"
    _verify_size(zip_path, min_bytes=1_000_000_000, label=f"SQUIRLS {suffix} zip")
    if not zip_path.exists():
        _download(urls["squirls_zip"], zip_path, f"SQUIRLS {suffix} (~4 GB)")

    print("  SQUIRLS zip 展開中...")
    with zipfile.ZipFile(zip_path) as zf:
        db_names = [n for n in zf.namelist() if n.endswith(".db") or n.endswith(".sqlite")]
        if not db_names:
            raise RuntimeError(f"SQUIRLS zip に .db ファイルが見つかりません: {zip_path}")
        for name in db_names:
            zf.extract(name, dest_dir)
            # zipfile が サブディレクトリ付きで展開した場合はフラットに移動
            extracted = dest_dir / name
            if extracted.parent != dest_dir:
                extracted.rename(dest_dir / extracted.name)

    zip_path.unlink(missing_ok=True)
    return True


# ---------------------------------------------------------------------------
# メイン
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="acmg-classifier 用データファイルの自動セットアップ",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--data-dir", type=Path, default=Path("./data"), metavar="PATH",
                        help="データ格納ルートディレクトリ (デフォルト: ./data)")
    parser.add_argument("--assembly", choices=["GRCh38", "GRCh37"], default="GRCh38")
    parser.add_argument("--genome-fasta", type=Path, default=None, metavar="PATH",
                        help="既存のゲノム FASTA パス (ダウンロードをスキップ)")
    parser.add_argument("--gnomad-vcf-dir", type=Path, default=None, metavar="PATH",
                        help="既存の gnomAD *.vcf.bgz ディレクトリ")
    parser.add_argument("--gnomad-chromosomes", nargs="+", default=None, metavar="CHR",
                        help="ダウンロードする染色体 (デフォルト: 全24本)")
    parser.add_argument("--gnomad-workers", type=int, default=None, metavar="N",
                        help="DuckDB 構築の並列ワーカー数 (デフォルト: CPU コア数-1)")
    parser.add_argument("--skip-gnomad", action="store_true",
                        help="gnomAD ダウンロードをスキップ (~300 GB)")
    parser.add_argument("--skip-genome", action="store_true",
                        help="ゲノム FASTA ダウンロードをスキップ (~880 MB)")
    parser.add_argument("--skip-vep-cache", action="store_true",
                        help="VEP キャッシュダウンロードをスキップ (~14 GB)")
    parser.add_argument("--skip-esm1b", action="store_true",
                        help="ESM1b ダウンロード・構築をスキップ (~1.34 GB)")
    parser.add_argument("--skip-squirls", action="store_true",
                        help="SQUIRLS DBダウンロードをスキップ (~4 GB、SpliceAI使用時など)")
    parser.add_argument("--squirls-db", type=Path, default=None, metavar="PATH",
                        help="既存の SQUIRLS *.db ファイルパス (ダウンロードをスキップ)")
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    assembly = args.assembly
    asm_dir = data_dir / assembly
    urls = URLS[assembly]
    chroms = args.gnomad_chromosomes or GNOMAD_CHROMS[assembly]

    sep = "=" * 62
    print(f"\n{sep}")
    print(f"  ACMG Classifier — データセットアップ")
    print(f"  データディレクトリ : {data_dir}")
    print(f"  アセンブリ         : {assembly}")
    print(f"{sep}\n")

    steps: list[tuple[str, object]] = [
        ("Reference genome",  lambda: step_genome(asm_dir, assembly, urls, args.genome_fasta, args.skip_genome)),
        ("VEP cache",         lambda: step_vep_cache(data_dir, assembly, urls, args.skip_vep_cache)),
        ("ClinVar VCF",       lambda: step_clinvar_vcf(asm_dir, assembly, urls)),
        ("ClinVar SQLite",    lambda: step_clinvar_sqlite(asm_dir, assembly, urls)),
        ("AlphaMissense",     lambda: step_alphamissense(asm_dir, assembly, urls)),
        ("ESM1b",             lambda: step_esm1b(data_dir, urls, args.skip_esm1b)),
        ("SQUIRLS",           lambda: step_squirls(asm_dir, assembly, urls, args.squirls_db, args.skip_squirls)),
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
            print("\n[中断しました]")
            sys.exit(130)
        except Exception as e:
            print(f"  [ERROR] {e}")
            failed_steps.append(name)
        print()

    print(sep)
    print(f"  完了: {len(ok_steps)} / {len(steps)} ステップ")
    if failed_steps:
        print(f"  失敗: {', '.join(failed_steps)}")
        sys.exit(1)
    else:
        print("  全ファイル準備完了")
    print(sep + "\n")


if __name__ == "__main__":
    main()
