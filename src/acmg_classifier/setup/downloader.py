"""Data download, verification, and status reporting."""
from __future__ import annotations
import hashlib
import urllib.request
from pathlib import Path

import structlog

from acmg_classifier.config import Config
from acmg_classifier.models.enums import Assembly

log = structlog.get_logger()

_REQUIRED_FILES_38 = [
    "genome/GRCh38.p14.fa",
    "gnomad/gnomad_v4.1_exomes.duckdb",
    "clinvar/clinvar_GRCh38.vcf.gz",
    "clinvar/clinvar_ps1_pm5_GRCh38.sqlite",
    "alphamissense/AlphaMissense_hg38.tsv.gz",
]
_REQUIRED_FILES_37 = [
    "genome/GRCh37.p13.fa",
    "gnomad/gnomad_v2.1.1_exomes.duckdb",
    "clinvar/clinvar_GRCh37.vcf.gz",
    "clinvar/clinvar_ps1_pm5_GRCh37.sqlite",
    "alphamissense/AlphaMissense_hg19.tsv.gz",
]
_REQUIRED_FILES = {
    Assembly.GRCH38: _REQUIRED_FILES_38,
    Assembly.GRCH37: _REQUIRED_FILES_37,
}


def validate_data_dir(cfg: Config) -> bool:
    ok = True
    for rel in _REQUIRED_FILES.get(cfg.assembly, []):
        p = cfg.assembly_dir / rel
        if not p.exists():
            log.warning("missing_data_file", path=str(p))
            ok = False
        else:
            log.info("data_file_ok", path=str(p))
    return ok


def print_status(data_dir: Path) -> None:
    from rich.console import Console
    from rich.table import Table
    console = Console()
    table = Table(title="Local Data Status")
    table.add_column("Assembly")
    table.add_column("File")
    table.add_column("Status")
    for asm in Assembly:
        for rel in _REQUIRED_FILES.get(asm, []):
            p = data_dir / asm.value / rel
            status = "[green]OK[/green]" if p.exists() else "[red]MISSING[/red]"
            table.add_row(asm.value, rel, status)
    console.print(table)


def run_setup(cfg: Config) -> None:
    """Guide user through data setup steps (downloads are large; user must run manually)."""
    from rich.console import Console
    console = Console()
    console.print("[bold]ACMG Classifier Data Setup[/bold]")
    console.print("Assembly: " + cfg.assembly.value)
    console.print("\nThe following data files are required (~60-65 GB per assembly).")
    console.print("Download and place them in: " + str(cfg.assembly_dir))
    console.print("\nSee the project documentation for download URLs and conversion scripts.")
    validate_data_dir(cfg)
