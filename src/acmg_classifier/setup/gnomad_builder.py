"""Build gnomAD DuckDB from per-chromosome VCF files using cyvcf2.

VCF のパース (cyvcf2) は GIL に縛られるため、染色体ごとにプロセス並列で
Parquet へ書き出し、最後に DuckDB で一括ロードする。
"""
from __future__ import annotations

import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import structlog

from acmg_classifier.utils.progress import progress_bar

log = structlog.get_logger()

_CREATE_TABLE = """
CREATE TABLE IF NOT EXISTS variants (
    chrom TEXT NOT NULL,
    pos INTEGER NOT NULL,
    ref TEXT NOT NULL,
    alt TEXT NOT NULL,
    af DOUBLE,
    an INTEGER,
    ac INTEGER,
    nhomalt INTEGER,
    nhemi INTEGER,
    popmax_af DOUBLE,
    popmax_pop TEXT,
    faf95_popmax DOUBLE,
    filters TEXT
);
"""

_CREATE_DB_INFO = """
CREATE TABLE IF NOT EXISTS db_info (
    gnomad_version TEXT,
    assembly TEXT,
    build_date TEXT
);
"""

_CREATE_INDEX = """
CREATE INDEX IF NOT EXISTS idx_variants ON variants (chrom, pos, ref, alt);
"""

_INSERT_SQL = "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"

_BATCH_SIZE = 50_000


def _to_scalar(val: Any) -> Any:
    """cyvcf2 は multi-allelic で tuple を返すので先頭要素を取り出す。"""
    if isinstance(val, (tuple, list)):
        return val[0] if val else None
    return val


def _to_int(val: Any) -> int | None:
    val = _to_scalar(val)
    if val is None:
        return None
    try:
        return int(val)
    except (ValueError, TypeError):
        return None


def _to_float(val: Any) -> float | None:
    val = _to_scalar(val)
    if val is None:
        return None
    try:
        f = float(val)
    except (ValueError, TypeError):
        return None
    # NaN は NULL 扱い
    if f != f:
        return None
    return f


def _to_str(val: Any) -> str | None:
    val = _to_scalar(val)
    if val is None:
        return None
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="replace")
    return str(val)


def _info(variant, *field_names: str) -> Any:
    """先に見つかった非 None の INFO フィールド値を返す (grpmax/popmax 両対応)。"""
    for f in field_names:
        try:
            val = variant.INFO.get(f)
        except (KeyError, TypeError):
            continue
        if val is not None:
            return val
    return None


def _detect_total_ram_gb() -> float | None:
    """物理メモリ総量 (GB) を best-effort で取得。失敗時は None。"""
    try:
        if sys.platform == "win32":
            import ctypes

            class _MEMORYSTATUSEX(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong),
                    ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]

            stat = _MEMORYSTATUSEX()
            stat.dwLength = ctypes.sizeof(_MEMORYSTATUSEX)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat))
            return stat.ullTotalPhys / (1024 ** 3)
        return (os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")) / (1024 ** 3)
    except Exception:
        return None


def _per_worker_mem_limit(max_workers: int) -> str:
    """各 worker の DuckDB memory_limit を決める。合計で実メモリの ~70% に収める。"""
    total = _detect_total_ram_gb()
    if total is None:
        return "2GB"
    per = max(1.0, total * 0.70 / max_workers)
    return f"{per:.1f}GB"


def build_gnomad_duckdb(
    vcf_dir: Path,
    output_db: Path,
    assembly: str,
    gnomad_version: str,
    max_workers: int | None = None,
) -> None:
    """
    gnomAD VCF (per-chromosome *.vcf.bgz) を DuckDB に取り込む。

    DuckDB コアには read_vcf が無いため cyvcf2 で 1 変異ずつ読み出す。
    染色体ごとにプロセス並列で Parquet へ書き出し、最後に read_parquet で
    一括ロードする。gnomAD v4.x の grpmax と v2.1.1 の popmax を自動で吸収する。
    """
    import duckdb
    from datetime import date

    output_db.parent.mkdir(parents=True, exist_ok=True)
    log.info("building_gnomad_db", output=str(output_db), assembly=assembly,
             version=gnomad_version)

    vcf_files = sorted(vcf_dir.glob("*.vcf.bgz"))
    if not vcf_files:
        raise FileNotFoundError(f"No *.vcf.bgz files found in {vcf_dir}")

    if max_workers is None:
        max_workers = max(1, (os.cpu_count() or 2) - 1)
    max_workers = max(1, min(max_workers, len(vcf_files)))

    mem_limit = _per_worker_mem_limit(max_workers)
    log.info("vcf_files_found", count=len(vcf_files), workers=max_workers,
             mem_per_worker=mem_limit)

    tmp_dir = output_db.parent / "_gnomad_parquet_tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir, ignore_errors=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    parquet_paths: list[Path] = []
    total_rows = 0
    try:
        # ---- 並列フェーズ: VCF → Parquet ----
        jobs = [
            (vcf_path, tmp_dir / f"{vcf_path.stem}.parquet")
            for vcf_path in vcf_files
        ]

        # Per-chromosome progress: each VCF file is one tick. The bar
        # makes it visible whether the bottleneck is one slow chromosome
        # (chr1/chr2 dominate runtime in gnomAD).
        with progress_bar("Loading gnomAD VCFs", total=len(jobs)) as advance:
            if max_workers == 1:
                for vcf_path, parquet_path in jobs:
                    name, rows = _vcf_to_parquet(
                        str(vcf_path), str(parquet_path), str(tmp_dir), mem_limit
                    )
                    parquet_paths.append(parquet_path)
                    total_rows += rows
                    log.info("vcf_loaded", file=name, rows=rows, cumulative=total_rows)
                    advance()
            else:
                with ProcessPoolExecutor(max_workers=max_workers) as pool:
                    future_to_job = {
                        pool.submit(
                            _vcf_to_parquet,
                            str(vcf_path), str(parquet_path), str(tmp_dir), mem_limit,
                        ): parquet_path
                        for vcf_path, parquet_path in jobs
                    }
                    for future in as_completed(future_to_job):
                        parquet_path = future_to_job[future]
                        name, rows = future.result()
                        parquet_paths.append(parquet_path)
                        total_rows += rows
                        log.info("vcf_loaded", file=name, rows=rows, cumulative=total_rows)
                        advance()

        # ---- マージフェーズ: Parquet → DuckDB ----
        log.info("loading_parquet_into_duckdb", files=len(parquet_paths))
        con = duckdb.connect(str(output_db))
        try:
            con.execute(_CREATE_TABLE)
            con.execute(_CREATE_DB_INFO)
            con.execute(
                "INSERT INTO db_info VALUES (?, ?, ?)",
                [gnomad_version, assembly, str(date.today())],
            )
            file_list = [str(p) for p in parquet_paths]
            con.execute(
                "INSERT INTO variants SELECT * FROM read_parquet($files)",
                {"files": file_list},
            )
            log.info("creating_index")
            con.execute(_CREATE_INDEX)
        finally:
            con.close()
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    log.info("gnomad_db_built", path=str(output_db), total_rows=total_rows)


def _vcf_to_parquet(
    vcf_path: str,
    parquet_path: str,
    temp_dir: str,
    mem_limit: str,
) -> tuple[str, int]:
    """1 つの VCF をパースし Parquet へ書き出す (worker プロセスで実行)。

    返り値は (ファイル名, 書き出した行数)。各 worker は専用の一時 DuckDB を
    使い、メモリ超過分は temp_dir へスピルする。
    """
    import duckdb
    from cyvcf2 import VCF

    name = Path(vcf_path).name
    build_db = parquet_path + ".build.duckdb"
    # 前回の残骸を除去
    for stale in (build_db, build_db + ".wal"):
        try:
            os.remove(stale)
        except FileNotFoundError:
            pass

    con = duckdb.connect(build_db)
    n_inserted = 0
    try:
        con.execute("SET threads TO 1")
        con.execute(f"SET memory_limit='{mem_limit}'")
        con.execute("SET preserve_insertion_order=false")
        con.execute(f"SET temp_directory='{temp_dir}'")
        con.execute(_CREATE_TABLE)

        vcf = VCF(vcf_path)
        batch: list[tuple] = []
        try:
            for v in vcf:
                # cyvcf2 returns None for PASS, otherwise the filter string ("AC0", ...)
                filter_val = v.FILTER if v.FILTER is not None else "PASS"

                chrom_raw = v.CHROM
                chrom = chrom_raw[3:] if chrom_raw.startswith("chr") else chrom_raw

                ref = v.REF
                for alt in v.ALT:
                    row = (
                        chrom,
                        int(v.POS),
                        ref,
                        alt,
                        _to_float(_info(v, "AF")),
                        _to_int(_info(v, "AN")),
                        _to_int(_info(v, "AC")),
                        _to_int(_info(v, "nhomalt")),
                        _to_int(_info(v, "nhemi")),
                        _to_float(_info(v, "AF_grpmax", "AF_popmax")),
                        _to_str(_info(v, "grpmax", "popmax")),
                        _to_float(_info(v, "faf95_grpmax", "faf95_popmax")),
                        filter_val,
                    )
                    batch.append(row)

                if len(batch) >= _BATCH_SIZE:
                    con.executemany(_INSERT_SQL, batch)
                    n_inserted += len(batch)
                    batch.clear()

            if batch:
                con.executemany(_INSERT_SQL, batch)
                n_inserted += len(batch)
        finally:
            vcf.close()

        con.execute(
            f"COPY variants TO '{parquet_path}' (FORMAT PARQUET, COMPRESSION ZSTD)"
        )
    finally:
        con.close()
        for stale in (build_db, build_db + ".wal"):
            try:
                os.remove(stale)
            except FileNotFoundError:
                pass

    return name, n_inserted
