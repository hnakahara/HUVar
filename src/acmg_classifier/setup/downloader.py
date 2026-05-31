"""Data download, verification, and status reporting."""
from __future__ import annotations
from pathlib import Path

import structlog

from acmg_classifier.config import Config
from acmg_classifier.models.enums import Assembly, InSilicoTool, SpliceTool

log = structlog.get_logger()

# Per-assembly files that are always required regardless of which in-silico
# tool is configured.
_BASE_FILES_38 = [
    "genome/GRCh38.p14.fa",
    "gnomad/gnomad_v4.1_exomes.duckdb",
    "clinvar/clinvar_GRCh38.vcf.gz",
    "clinvar/clinvar_ps1_pm5_GRCh38.sqlite",
]
_BASE_FILES_37 = [
    "genome/GRCh37.p13.fa",
    "gnomad/gnomad_v2.1.1_exomes.duckdb",
    "clinvar/clinvar_GRCh37.vcf.gz",
    "clinvar/clinvar_ps1_pm5_GRCh37.sqlite",
]
_BASE_FILES = {
    Assembly.GRCH38: _BASE_FILES_38,
    Assembly.GRCH37: _BASE_FILES_37,
}
_AM_FILE = {
    Assembly.GRCH38: "alphamissense/AlphaMissense_hg38.tsv.gz",
    Assembly.GRCH37: "alphamissense/AlphaMissense_hg19.tsv.gz",
}
# ESM1b SQLite is protein-coordinate, so it lives under data_dir/ (not
# assembly_dir/) and is shared across GRCh37/GRCh38.
_ESM1B_REL = "esm1b/esm1b_llr.sqlite"

# SQUIRLS DB directories (assembly-specific, user must place *.db file inside)
_SQUIRLS_DIR = {
    Assembly.GRCH38: "squirls/squirls-2309-hg38",
    Assembly.GRCH37: "squirls/squirls-2309-hg19",
}


def _squirls_db_exists(assembly_dir: Path, assembly: Assembly) -> bool:
    """Return True if at least one SQUIRLS *.db / *.sqlite file is present."""
    d = assembly_dir / _SQUIRLS_DIR[assembly]
    if not d.exists():
        return False
    return any(d.glob("*.db")) or any(d.glob("*.sqlite"))


def validate_data_dir(cfg: Config) -> bool:
    """Check every file the configured pipeline will need is present."""
    ok = True
    for rel in _BASE_FILES.get(cfg.assembly, []):
        p = cfg.assembly_dir / rel
        if not p.exists():
            log.warning("missing_data_file", path=str(p))
            ok = False
        else:
            log.info("data_file_ok", path=str(p))

    if cfg.insilico_tool == InSilicoTool.ESM1B:
        p = cfg.data_dir / _ESM1B_REL
        if not p.exists():
            log.warning("missing_data_file", path=str(p))
            ok = False
        else:
            log.info("data_file_ok", path=str(p))
    else:
        rel = _AM_FILE[cfg.assembly]
        p = cfg.assembly_dir / rel
        if not p.exists():
            log.warning("missing_data_file", path=str(p))
            ok = False
        else:
            log.info("data_file_ok", path=str(p))

    if cfg.splice_tool == SpliceTool.SQUIRLS:
        if not _squirls_db_exists(cfg.assembly_dir, cfg.assembly):
            log.warning(
                "missing_squirls_db",
                path=str(cfg.assembly_dir / _SQUIRLS_DIR[cfg.assembly]),
                hint="Place the SQUIRLS *.db file inside this directory. "
                     "Download from https://squirls.readthedocs.io/",
            )
            ok = False
        else:
            log.info("squirls_db_ok", path=str(cfg.assembly_dir / _SQUIRLS_DIR[cfg.assembly]))

    return ok


def print_status(data_dir: Path) -> None:
    from rich.console import Console
    from rich.table import Table
    console = Console()
    table = Table(title="Local Data Status")
    table.add_column("Assembly")
    table.add_column("File / Directory")
    table.add_column("Status")
    for asm in Assembly:
        asm_dir = data_dir / asm.value
        for rel in _BASE_FILES.get(asm, []) + [_AM_FILE[asm]]:
            p = asm_dir / rel
            status = "[green]OK[/green]" if p.exists() else "[red]MISSING[/red]"
            table.add_row(asm.value, rel, status)
        # SQUIRLS (optional splice tool)
        squirls_ok = _squirls_db_exists(asm_dir, asm)
        status = "[green]OK[/green]" if squirls_ok else "[yellow]OPTIONAL/MISSING[/yellow]"
        table.add_row(asm.value, _SQUIRLS_DIR[asm] + "/*.db", status)
    # ESM1b is assembly-independent; show it once.
    p = data_dir / _ESM1B_REL
    status = "[green]OK[/green]" if p.exists() else "[yellow]OPTIONAL/MISSING[/yellow]"
    table.add_row("(shared)", _ESM1B_REL, status)
    console.print(table)


def run_setup(cfg: Config) -> None:
    """Guide user through data setup steps (downloads are large; user must run manually)."""
    from rich.console import Console
    console = Console()
    console.print("[bold]ACMG Classifier Data Setup[/bold]")
    console.print("Assembly: " + cfg.assembly.value)
    console.print("\nThe following data files are required (~60-65 GB per assembly).")
    console.print("Download and place them in: " + str(cfg.assembly_dir))
    if cfg.splice_tool == SpliceTool.SQUIRLS:
        squirls_dir = cfg.assembly_dir / _SQUIRLS_DIR[cfg.assembly]
        console.print(
            "\n[yellow]SQUIRLS DB:[/yellow] Download the precomputed SQUIRLS database "
            f"(squirls-2309-hg38 or hg19) and place the *.db file in:\n  {squirls_dir}"
        )
        console.print("  Download: https://squirls.readthedocs.io/en/master/setup.html")
    console.print("\nSee the project documentation for download URLs and conversion scripts.")
    validate_data_dir(cfg)
