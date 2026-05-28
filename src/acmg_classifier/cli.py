from __future__ import annotations
import sys
from pathlib import Path
from typing import Optional

import click
import structlog

from acmg_classifier.models.enums import Assembly, InSilicoTool, SpliceTool

log = structlog.get_logger()


def _make_config(ctx: click.Context, **kwargs):
    """Build Config from CLI options, injecting into Click context."""
    from acmg_classifier.config import Config
    cfg = Config(**{k: v for k, v in kwargs.items() if v is not None})
    ctx.ensure_object(dict)
    ctx.obj["config"] = cfg
    return cfg


@click.group()
@click.version_option()
def cli():
    """ACMG 2015 + ClinGen SVI variant pathogenicity classifier (fully local)."""
    structlog.configure(
        processors=[
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.dev.ConsoleRenderer(),
        ]
    )


# ---------------------------------------------------------------------------
# classify
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("vcf", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path), default=None,
              help="Output TSV path (default: stdout)")
@click.option("--data-dir", type=click.Path(path_type=Path), default=Path("./data"),
              show_default=True)
@click.option("--assembly", type=click.Choice([a.value for a in Assembly]),
              default=None, help="Override assembly detection from VCF header")
@click.option("--insilico-tool",
              type=click.Choice([t.value for t in InSilicoTool]),
              default=InSilicoTool.ALPHAMISSENSE.value, show_default=True)
@click.option("--splice-tool",
              type=click.Choice([t.value for t in SpliceTool]),
              default=SpliceTool.SQUIRLS.value, show_default=True)
@click.option("--spliceai-dir", type=click.Path(path_type=Path), default=None,
              help="Directory containing SpliceAI VCF files (snv + indel). Overrides default data-dir lookup.")
@click.option("--supplement", type=click.Path(exists=True, path_type=Path),
              default=None, help="Manual evidence TSV")
@click.option("--inheritance",
              type=click.Choice(["AD", "AR", "XL", "Unknown"]),
              default="Unknown", show_default=True)
@click.option("--workers", type=int, default=4, show_default=True)
@click.pass_context
def classify(
    ctx: click.Context,
    vcf: Path,
    output: Optional[Path],
    data_dir: Path,
    assembly: Optional[str],
    insilico_tool: str,
    splice_tool: str,
    spliceai_dir: Optional[Path],
    supplement: Optional[Path],
    inheritance: str,
    workers: int,
) -> None:
    """Classify all variants in VCF_FILE and write results to TSV."""
    from acmg_classifier.config import Config
    from acmg_classifier.pipeline.pipeline import run_pipeline

    cfg = Config(
        data_dir=data_dir,
        assembly=Assembly(assembly) if assembly else Assembly.GRCH38,
        insilico_tool=InSilicoTool(insilico_tool),
        splice_tool=SpliceTool(splice_tool),
        spliceai_dir=spliceai_dir,
        workers=workers,
    )
    run_pipeline(vcf, cfg, output_path=output, supplement_path=supplement)


# ---------------------------------------------------------------------------
# explain
# ---------------------------------------------------------------------------

@cli.command()
@click.argument("chrom")
@click.argument("pos", type=int)
@click.argument("ref")
@click.argument("alt")
@click.option("--data-dir", type=click.Path(path_type=Path), default=Path("./data"))
@click.option("--assembly", type=click.Choice([a.value for a in Assembly]),
              default=Assembly.GRCH38.value, show_default=True)
@click.pass_context
def explain(
    ctx: click.Context,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    data_dir: Path,
    assembly: str,
) -> None:
    """Show detailed classification for a single variant (CHROM POS REF ALT)."""
    from acmg_classifier.config import Config
    from acmg_classifier.pipeline.pipeline import run_single

    cfg = Config(data_dir=data_dir, assembly=Assembly(assembly))
    run_single(chrom, pos, ref, alt, cfg)


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

@cli.command()
@click.option("--data-dir", type=click.Path(path_type=Path), default=Path("./data"))
@click.option("--assembly",
              type=click.Choice([a.value for a in Assembly] + ["both"]),
              default=Assembly.GRCH38.value, show_default=True)
def validate(data_dir: Path, assembly: str) -> None:
    """Check that all required data files exist and are accessible."""
    from acmg_classifier.setup.downloader import validate_data_dir

    assemblies = list(Assembly) if assembly == "both" else [Assembly(assembly)]
    ok = True
    for asm in assemblies:
        from acmg_classifier.config import Config
        cfg = Config(data_dir=data_dir, assembly=asm)
        ok &= validate_data_dir(cfg)
    sys.exit(0 if ok else 1)


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------

@cli.command()
@click.option("--data-dir", type=click.Path(path_type=Path), default=Path("./data"))
def status(data_dir: Path) -> None:
    """Show status of local databases (gnomAD version, ClinVar date, etc.)."""
    from acmg_classifier.setup.downloader import print_status
    print_status(data_dir)


# ---------------------------------------------------------------------------
# setup
# ---------------------------------------------------------------------------

@cli.command()
@click.option("--data-dir", type=click.Path(path_type=Path), default=Path("./data"))
@click.option("--assembly",
              type=click.Choice([a.value for a in Assembly] + ["both"]),
              default=Assembly.GRCH38.value, show_default=True)
def setup(data_dir: Path, assembly: str) -> None:
    """Download and build all required local data files."""
    from acmg_classifier.setup.downloader import run_setup

    assemblies = list(Assembly) if assembly == "both" else [Assembly(assembly)]
    for asm in assemblies:
        from acmg_classifier.config import Config
        cfg = Config(data_dir=data_dir, assembly=asm)
        run_setup(cfg)
