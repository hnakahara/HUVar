from __future__ import annotations
import logging
import sys
from pathlib import Path
from typing import Optional

import click
import structlog

from acmg_classifier.models.enums import Assembly, InSilicoTool, SpliceTool, SupplementMode

log = structlog.get_logger()


def _make_config(ctx: click.Context, **kwargs):
    """Build Config from CLI options, injecting into Click context.

    NOTE: currently unused — each subcommand constructs Config inline.
    See docs/cleanup-candidates.md."""
    from acmg_classifier.config import Config
    cfg = Config(**{k: v for k, v in kwargs.items() if v is not None})
    ctx.ensure_object(dict)
    ctx.obj["config"] = cfg
    return cfg


@click.group()
@click.version_option()
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    default="INFO", show_default=True,
    help="Console log verbosity. DEBUG shows per-batch VEP commands and other "
         "low-level detail; the default INFO suppresses them.",
)
def cli(log_level: str):
    """ACMG 2015 + ClinGen SVI variant pathogenicity classifier (fully local).

    Pass --log-level before the subcommand, e.g.
    `acmg-classify --log-level DEBUG classify ...`.
    """
    # Structured console logging is configured here (not at import time) so
    # library users importing acmg_classifier modules don't get our logging
    # config forced onto their own application.
    #
    # make_filtering_bound_logger drops records below the chosen level at the
    # bound-logger layer (cheaply, before processors run). Without it structlog
    # emits every level — which is why DEBUG-level `vep_exec` lines appeared.
    structlog.configure(
        processors=[
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.dev.ConsoleRenderer(),
        ],
        wrapper_class=structlog.make_filtering_bound_logger(
            getattr(logging, log_level.upper())
        ),
    )


# ---------------------------------------------------------------------------
# classify
# ---------------------------------------------------------------------------

def _with_cspec(cfg, cspec):
    """Apply --cspec GENE=CSPEC_ID selections (cspec_overlay) or exit with the choices."""
    if not cspec:
        return cfg
    from acmg_classifier.criteria.cspec_overlay import apply_cspec_selection, parse_cspec_options
    try:
        sel = parse_cspec_options(cspec, cfg.disease_prevalence_multispec_tsv)
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="--cspec")
    return apply_cspec_selection(cfg, sel)


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
              default=InSilicoTool.ESM1B.value, show_default=True)
@click.option("--with-bayesdel", is_flag=True, default=False,
              help="Consult the opt-in BayesDel scores (must be staged via "
                   "setup_data.py --with-bayesdel). Licence-gated: applied only when "
                   "--insilico-tool is revel or alphamissense; ignored under esm1b.")
@click.option("--with-cadd", is_flag=True, default=False,
              help="Consult the opt-in CADD scores (must be staged via "
                   "setup_data.py --with-cadd). Licence-gated: applied only when "
                   "--insilico-tool is revel or alphamissense; ignored under esm1b.")
@click.option("--splice-tool",
              type=click.Choice([SpliceTool.OPENSPLICEAI.value, SpliceTool.SPLICEAI.value]),
              default=SpliceTool.OPENSPLICEAI.value, show_default=True,
              help="Splice predictor. 'openspliceai' (default) runs OSAI_MANE at inference time "
                   "(pip install openspliceai); 'spliceai' uses precomputed VCFs (Illumina licence required).")
@click.option("--openspliceai-model-dir", type=click.Path(path_type=Path), default=None,
              help="Directory containing OSAI_MANE model files. "
                   "Overrides default data/<assembly>/openspliceai/<flanking-size>nt/ lookup.")
@click.option("--openspliceai-flanking-size", type=int, default=80, show_default=True,
              help="Model context length (nt). Default 80 for fast CPU inference; raise to "
                   "400/2000/10000 for higher splice sensitivity. Must match a downloaded "
                   "OSAI_MANE model variant. Also selects the default model subdirectory.")
@click.option("--spliceai-dir", type=click.Path(path_type=Path), default=None,
              help="Directory containing SpliceAI VCF files (snv + indel). Overrides default data-dir lookup.")
@click.option("--supplement", type=click.Path(exists=True, path_type=Path),
              default=None, help="Manual evidence TSV")
@click.option("--supplement-mode",
              type=click.Choice([m.value for m in SupplementMode]),
              default=SupplementMode.MERGE.value, show_default=True,
              help="How --supplement combines with the tool's calls. 'merge': curator "
                   "entries override the strength of any criterion they name and add "
                   "criteria the tool left not-met. 'manual-only': for variants listed in "
                   "the supplement, classify purely from it; variants not listed fall back "
                   "to the tool's automated calls.")
@click.option("--with-tp53-functional", is_flag=True, default=False,
              help="TP53 PS3/BS3 from the VCEP functional flowchart on the NCI TP53 "
                   "Database assay classes (bundled tp53_functional.tsv). Off by default: "
                   "the database terms prohibit use that results in monetization. "
                   "Env: ACMG_USE_TP53_FUNCTIONAL=true.")
@click.option("--cspec", "cspec", multiple=True, metavar="GENE=CSPEC_ID",
              help="Evaluate GENE under one of its disease-specific ClinGen CSpecs "
                   "(genes with several released specifications: RYR1, ACTA1, VWF), e.g. "
                   "--cspec RYR1=congenital_myopathy_ar (repeatable). Without it, the "
                   "conservative merged specification is used. Choices: "
                   "malignant_hyperthermia, congenital_myopathy_ad, congenital_myopathy_ar "
                   "(RYR1); congenital_myopathy_ad, congenital_myopathy_ar (ACTA1); "
                   "vwd_ad, vwd_ar (VWF).")
@click.option("--exclude-self-expert-panel", is_flag=True, default=False,
              help="Benchmark mode: do not import PS3/PS4/PP1/BS2 from the ClinVar "
                   "expert-panel (>=3-star) record of the variant being classified, so the "
                   "evaluation does not read the reference curation itself. "
                   "Env: ACMG_EXCLUDE_SELF_EXPERT_PANEL=true.")
@click.option("--workers", type=int, default=4, show_default=True)
@click.option("--no-progress", is_flag=True, default=False,
              help="Disable progress bars (auto-disabled when stderr is not a TTY).")
@click.option("--limit", type=int, default=None,
              help="Process only the first N classifiable variants. Use with --profile for fast profiling.")
@click.option("--profile", "profile_path", type=click.Path(path_type=Path), default=None,
              help="Write cProfile stats to PATH and print top 30 functions to stderr.")
@click.pass_context
def classify(
    ctx: click.Context,
    vcf: Path,
    output: Optional[Path],
    data_dir: Path,
    assembly: Optional[str],
    insilico_tool: str,
    with_bayesdel: bool,
    with_cadd: bool,
    splice_tool: str,
    openspliceai_model_dir: Optional[Path],
    openspliceai_flanking_size: int,
    spliceai_dir: Optional[Path],
    supplement: Optional[Path],
    supplement_mode: str,
    with_tp53_functional: bool,
    cspec: tuple[str, ...],
    exclude_self_expert_panel: bool,
    workers: int,
    no_progress: bool,
    limit: Optional[int],
    profile_path: Optional[Path],
) -> None:
    """Classify all variants in VCF_FILE and write results to TSV."""
    from acmg_classifier.config import Config
    from acmg_classifier.pipeline.pipeline import run_pipeline
    from acmg_classifier.utils import progress

    if no_progress:
        progress.set_enabled(False)

    cfg = Config(
        data_dir=data_dir,
        assembly=Assembly(assembly) if assembly else Assembly.GRCH38,
        insilico_tool=InSilicoTool(insilico_tool),
        use_bayesdel=with_bayesdel,
        use_cadd=with_cadd,
        splice_tool=SpliceTool(splice_tool),
        openspliceai_model_dir=openspliceai_model_dir,
        openspliceai_flanking_size=openspliceai_flanking_size,
        spliceai_dir=spliceai_dir,
        supplement_mode=SupplementMode(supplement_mode),
        workers=workers,
        # Only override when the flag is given, so the env variable still works.
        **({"exclude_self_expert_panel": True} if exclude_self_expert_panel else {}),
        **({"use_tp53_functional": True} if with_tp53_functional else {}),
    )
    cfg = _with_cspec(cfg, cspec)
    run_pipeline(
        vcf, cfg,
        output_path=output,
        supplement_path=supplement,
        limit=limit,
        profile_path=profile_path,
    )


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
@click.option("--insilico-tool",
              type=click.Choice([t.value for t in InSilicoTool]),
              default=InSilicoTool.ESM1B.value, show_default=True)
@click.option("--with-bayesdel", is_flag=True, default=False,
              help="Consult the opt-in BayesDel scores (staged via setup_data.py "
                   "--with-bayesdel). Applied only under --insilico-tool revel/alphamissense; "
                   "ignored under esm1b.")
@click.option("--with-cadd", is_flag=True, default=False,
              help="Consult the opt-in CADD scores (staged via setup_data.py --with-cadd). "
                   "Applied only under --insilico-tool revel/alphamissense; ignored under esm1b.")
@click.option("--splice-tool",
              type=click.Choice([SpliceTool.OPENSPLICEAI.value, SpliceTool.SPLICEAI.value]),
              default=SpliceTool.OPENSPLICEAI.value, show_default=True,
              help="Splice predictor. 'openspliceai' (default) runs OSAI_MANE at inference time "
                   "(pip install openspliceai); 'spliceai' uses precomputed VCFs (Illumina licence required).")
@click.option("--openspliceai-model-dir", type=click.Path(path_type=Path), default=None,
              help="Directory containing OSAI_MANE model files.")
@click.option("--openspliceai-flanking-size", type=int, default=80, show_default=True,
              help="Model context length (nt). Default 80 for fast CPU inference; raise to "
                   "400/2000/10000 for higher splice sensitivity. Must match a downloaded "
                   "OSAI_MANE model variant.")
@click.option("--spliceai-dir", type=click.Path(path_type=Path), default=None,
              help="Directory containing SpliceAI VCF files (snv + indel).")
@click.option("--evidence", "evidence", multiple=True, metavar="CRITERION:STRENGTH[:NOTE]",
              help="Add manual evidence to this variant (repeatable), e.g. "
                   "--evidence PS3:strong --evidence PM1:moderate:hotspot. "
                   "Criterion (PVS1/PS3/PM1/BA1/...) and strength "
                   "(very_strong/strong/moderate/supporting, or stand_alone for BA1) "
                   "are case-insensitive; an optional 3rd field is a free-text note.")
@click.option("--supplement", type=click.Path(exists=True, path_type=Path),
              default=None, help="Manual evidence TSV (rows matching this variant are applied).")
@click.option("--supplement-mode",
              type=click.Choice([m.value for m in SupplementMode]),
              default=SupplementMode.MERGE.value, show_default=True,
              help="How manual evidence combines with automated calls. 'merge': curator "
                   "entries override the strength of any criterion they name and add "
                   "criteria left not-met. 'manual-only': classify purely from the "
                   "supplied evidence.")
@click.option("--with-tp53-functional", is_flag=True, default=False,
              help="TP53 PS3/BS3 from the VCEP functional flowchart on the NCI TP53 "
                   "Database assay classes (bundled tp53_functional.tsv). Off by default: "
                   "the database terms prohibit use that results in monetization. "
                   "Env: ACMG_USE_TP53_FUNCTIONAL=true.")
@click.option("--cspec", "cspec", multiple=True, metavar="GENE=CSPEC_ID",
              help="Evaluate GENE under one of its disease-specific ClinGen CSpecs "
                   "(genes with several released specifications: RYR1, ACTA1, VWF), e.g. "
                   "--cspec RYR1=congenital_myopathy_ar (repeatable). Without it, the "
                   "conservative merged specification is used. Choices: "
                   "malignant_hyperthermia, congenital_myopathy_ad, congenital_myopathy_ar "
                   "(RYR1); congenital_myopathy_ad, congenital_myopathy_ar (ACTA1); "
                   "vwd_ad, vwd_ar (VWF).")
@click.option("--exclude-self-expert-panel", is_flag=True, default=False,
              help="Benchmark mode: do not import PS3/PS4/PP1/BS2 from the ClinVar "
                   "expert-panel (>=3-star) record of the variant being classified, so the "
                   "evaluation does not read the reference curation itself. "
                   "Env: ACMG_EXCLUDE_SELF_EXPERT_PANEL=true.")
@click.option("--transcript", default=None, metavar="ACCESSION",
              help="Evaluate the criteria on this transcript instead of MANE Select, "
                   "e.g. the gene's MANE Plus Clinical (--transcript NM_033360.4). VCEP "
                   "position rules and ClinVar codon comparisons still use MANE Select "
                   "numbering (reported in the warnings).")
@click.pass_context
def explain(
    ctx: click.Context,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    data_dir: Path,
    assembly: str,
    insilico_tool: str,
    with_bayesdel: bool,
    with_cadd: bool,
    splice_tool: str,
    openspliceai_model_dir: Optional[Path],
    openspliceai_flanking_size: int,
    spliceai_dir: Optional[Path],
    evidence: tuple[str, ...],
    supplement: Optional[Path],
    supplement_mode: str,
    with_tp53_functional: bool,
    cspec: tuple[str, ...],
    exclude_self_expert_panel: bool,
    transcript: Optional[str],
) -> None:
    """Show detailed classification for a single variant (CHROM POS REF ALT).

    Manual evidence can be added inline with --evidence (repeatable) and/or from
    a --supplement TSV; both are combined with the automated calls per
    --supplement-mode.
    """
    from acmg_classifier.config import Config
    from acmg_classifier.io.supplement_reader import parse_inline_evidence, read_supplement
    from acmg_classifier.pipeline.pipeline import run_single
    from acmg_classifier.models.variant import VariantRecord

    cfg = Config(
        data_dir=data_dir,
        assembly=Assembly(assembly),
        insilico_tool=InSilicoTool(insilico_tool),
        use_bayesdel=with_bayesdel,
        use_cadd=with_cadd,
        splice_tool=SpliceTool(splice_tool),
        openspliceai_model_dir=openspliceai_model_dir,
        openspliceai_flanking_size=openspliceai_flanking_size,
        spliceai_dir=spliceai_dir,
        supplement_mode=SupplementMode(supplement_mode),
        **({"exclude_self_expert_panel": True} if exclude_self_expert_panel else {}),
        **({"use_tp53_functional": True} if with_tp53_functional else {}),
    )
    cfg = _with_cspec(cfg, cspec)

    # Assemble manual evidence for this variant: TSV rows keyed to it + inline.
    variant_key = VariantRecord(
        chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=cfg.assembly
    ).key
    entries = []
    if supplement:
        entries.extend(read_supplement(supplement).get(variant_key, []))
    if evidence:
        entries.extend(parse_inline_evidence(list(evidence), variant_key))

    run_single(chrom, pos, ref, alt, cfg, supplement=entries or None, transcript=transcript)


# ---------------------------------------------------------------------------
# somatic (oncogenicity, ClinGen/CGC/VICC SOP)
# ---------------------------------------------------------------------------

def _somatic_options(f):
    """Annotation options shared by `somatic` and `somatic-explain`."""
    opts = [
        click.option("--data-dir", type=click.Path(path_type=Path), default=Path("./data"),
                     show_default=True),
        click.option("--insilico-tool", type=click.Choice([t.value for t in InSilicoTool]),
                     default=InSilicoTool.ESM1B.value, show_default=True,
                     help="Predictor behind OP1 / SBP1 (HUVar PP3 / BP4)."),
        click.option("--splice-tool",
                     type=click.Choice([SpliceTool.OPENSPLICEAI.value, SpliceTool.SPLICEAI.value]),
                     default=SpliceTool.OPENSPLICEAI.value, show_default=True),
        click.option("--openspliceai-model-dir", type=click.Path(path_type=Path), default=None),
        click.option("--openspliceai-flanking-size", type=int, default=80, show_default=True),
        click.option("--spliceai-dir", type=click.Path(path_type=Path), default=None),
        click.option("--with-tp53-functional", is_flag=True, default=False,
                     help="Use the NCI TP53 Database assay classes for OS2/SBS2 in TP53 "
                          "(non-commercial use only; default: calibrated MaveDB data)."),
    ]
    for o in reversed(opts):
        f = o(f)
    return f


def _somatic_config(data_dir, assembly, insilico_tool, splice_tool, openspliceai_model_dir,
                    openspliceai_flanking_size, spliceai_dir, with_tp53_functional, **extra):
    from acmg_classifier.config import Config
    return Config(
        data_dir=data_dir,
        assembly=Assembly(assembly) if assembly else Assembly.GRCH38,
        insilico_tool=InSilicoTool(insilico_tool),
        splice_tool=SpliceTool(splice_tool),
        openspliceai_model_dir=openspliceai_model_dir,
        openspliceai_flanking_size=openspliceai_flanking_size,
        spliceai_dir=spliceai_dir,
        **({"use_tp53_functional": True} if with_tp53_functional else {}),
        **extra,
    )


@cli.command()
@click.argument("vcf", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path), default=None,
              help="Output TSV path (default: stdout)")
@click.option("--assembly", type=click.Choice([a.value for a in Assembly]), default=None,
              help="Override assembly detection from VCF header")
@_somatic_options
@click.option("--supplement", type=click.Path(exists=True, path_type=Path), default=None,
              help="Manual somatic evidence TSV: variant_id, criterion (OVS1..SBP2), "
                   "evidence[, met=yes/no]. Sets or overrides a criterion (e.g. OP2, OS2).")
@click.option("--tumor-type", default=None, metavar="ONCOTREE_CODE",
              help="Tumour type (OncoTree code, e.g. SKCM, LUAD, COADREAD) for the AMP/ASCO/CAP "
                   "2017 tiers. Without it no tier is assigned and evidence in other tumour "
                   "types is listed for reference.")
@click.option("--workers", type=int, default=4, show_default=True)
@click.option("--no-progress", is_flag=True, default=False)
def somatic(vcf, output, assembly, supplement, tumor_type, workers, no_progress, **opts) -> None:
    """Classify the oncogenicity of somatic variants in VCF (Horak et al. 2022).

    Points-based ClinGen/CGC/VICC SOP: Oncogenic / Likely oncogenic / VUS /
    Likely benign / Benign. Clinical-significance tiers (AMP/ASCO/CAP 2017) are
    not assigned. Uses only licence-free data (NCG gene roles, cancerhotspots.org,
    ClinVar, gnomAD, in-silico predictors, MaveDB)."""
    from acmg_classifier.somatic.pipeline import run_somatic_pipeline
    from acmg_classifier.utils import progress

    if no_progress:
        progress.set_enabled(False)
    cfg = _somatic_config(assembly=assembly, workers=workers, **opts)
    _check_tumor_type(cfg, tumor_type)
    run_somatic_pipeline(vcf, cfg, output_path=output, supplement_path=supplement,
                         tumor_type=tumor_type)


def _check_tumor_type(cfg, code):
    if not code:
        return
    from acmg_classifier.somatic.tumor import tumor_type
    try:
        tumor_type(code, cfg)
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="--tumor-type")


@cli.command("somatic-explain")
@click.argument("chrom")
@click.argument("pos", type=int)
@click.argument("ref")
@click.argument("alt")
@click.option("--assembly", type=click.Choice([a.value for a in Assembly]),
              default=Assembly.GRCH38.value, show_default=True)
@_somatic_options
@click.option("--evidence", "evidence", multiple=True, metavar="CRITERION[:not_met][:NOTE]",
              help="Manual somatic evidence (repeatable), e.g. --evidence OP2:retinoblastoma "
                   "--evidence OS2:'PMID 12345' --evidence OM1:not_met.")
@click.option("--supplement", type=click.Path(exists=True, path_type=Path), default=None,
              help="Manual somatic evidence TSV (rows for this variant are applied).")
@click.option("--transcript", default=None, metavar="ACCESSION",
              help="Evaluate on this transcript (e.g. MANE Plus Clinical) instead of MANE Select.")
@click.option("--tumor-type", default=None, metavar="ONCOTREE_CODE",
              help="Tumour type (OncoTree code, e.g. SKCM) for the AMP/ASCO/CAP 2017 tier.")
def somatic_explain(chrom, pos, ref, alt, assembly, evidence, supplement, transcript,
                    tumor_type, **opts) -> None:
    """Show the somatic oncogenicity evaluation of one variant (CHROM POS REF ALT)."""
    from acmg_classifier.models.variant import VariantRecord
    from acmg_classifier.somatic.pipeline import (
        parse_inline_evidence, read_somatic_supplement, run_somatic_single,
    )

    cfg = _somatic_config(assembly=assembly, **opts)
    key = VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=cfg.assembly).key
    manual = []
    if supplement:
        manual.extend(read_somatic_supplement(supplement).get(key, []))
    if evidence:
        try:
            manual.extend(parse_inline_evidence(list(evidence), key))
        except ValueError as exc:
            raise click.BadParameter(str(exc), param_hint="--evidence")
    _check_tumor_type(cfg, tumor_type)
    run_somatic_single(chrom, pos, ref, alt, cfg, manual=manual or None, transcript=transcript,
                       tumor_type=tumor_type)


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
@click.option("--no-progress", is_flag=True, default=False,
              help="Disable progress bars (auto-disabled when stderr is not a TTY).")
def setup(data_dir: Path, assembly: str, no_progress: bool) -> None:
    """Download and build all required local data files."""
    from acmg_classifier.setup.downloader import run_setup
    from acmg_classifier.utils import progress

    if no_progress:
        progress.set_enabled(False)

    assemblies = list(Assembly) if assembly == "both" else [Assembly(assembly)]
    for asm in assemblies:
        from acmg_classifier.config import Config
        cfg = Config(data_dir=data_dir, assembly=asm)
        run_setup(cfg)
