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

# MMSplice gene-annotation GTF (assembly-specific). Must match config.mmsplice_gtf.
# Fetched + protein-coding-filtered by scripts/setup_data.py.
_MMSPLICE_GTF = {
    Assembly.GRCH38: "mmsplice/Homo_sapiens.GRCh38.111.protein_coding.gtf",
    Assembly.GRCH37: "mmsplice/Homo_sapiens.GRCh37.87.protein_coding.gtf",
}


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

    # Splice tool data checks are intentionally absent: SQUIRLS (default) is on
    # hold (no downloadable DB → predictor unavailable → splice scoring skipped),
    # SpliceAI scores are licence-gated and placed manually, and MMSplice is
    # disabled. MMSplice GTF check retained, commented out, for later:
    # if cfg.splice_tool == SpliceTool.MMSPLICE:
    #     p = cfg.mmsplice_gtf
    #     if not p.exists():
    #         log.warning("missing_mmsplice_gtf", path=str(p),
    #                     hint="Run setup to fetch the GTF, and: pip install -e .[mmsplice]")
    #         ok = False
    #     else:
    #         log.info("mmsplice_gtf_ok", path=str(p))

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
        # MMSplice GTF row DISABLED (MMSplice integration is off). Retained:
        # gtf_p = asm_dir / _MMSPLICE_GTF[asm]
        # status = "[green]OK[/green]" if gtf_p.exists() else "[yellow]OPTIONAL/MISSING[/yellow]"
        # table.add_row(asm.value, _MMSPLICE_GTF[asm], status)
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
    # MMSplice setup guidance DISABLED (integration is off). Retained:
    # if cfg.splice_tool == SpliceTool.MMSPLICE:
    #     console.print(
    #         "\n[yellow]MMSplice:[/yellow] install the optional dependency "
    #         "(pip install -e .[mmsplice]) for runtime splice scoring. The "
    #         f"protein-coding GTF is fetched automatically to:\n  {cfg.mmsplice_gtf}"
    #     )
    console.print("\nSee the project documentation for download URLs and conversion scripts.")
    validate_data_dir(cfg)
