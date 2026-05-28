"""Rich terminal report for the explain command."""
from __future__ import annotations
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich import box

from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.classification import ClassificationResult
from acmg_classifier.models.enums import Pathogenicity
from acmg_classifier.models.variant import VariantRecord

_PATH_COLORS = {
    Pathogenicity.PATHOGENIC: "bold red",
    Pathogenicity.LIKELY_PATHOGENIC: "red",
    Pathogenicity.VUS: "yellow",
    Pathogenicity.LIKELY_BENIGN: "blue",
    Pathogenicity.BENIGN: "bold blue",
}


def print_report(
    result: ClassificationResult,
    variant: VariantRecord,
    annotation: AnnotationData,
) -> None:
    console = Console()
    console.print()
    console.print(Panel(
        "[bold]" + result.variant_id + "[/bold]",
        title="ACMG Variant Classification Report",
        expand=False,
    ))

    cls2015_color = _PATH_COLORS.get(result.classification_2015, "white")
    clsbay_color = _PATH_COLORS.get(result.classification_bayesian, "white")

    summary = Table(box=box.SIMPLE, show_header=False)
    summary.add_column("Field", style="bold")
    summary.add_column("Value")
    summary.add_row("Assembly", variant.assembly.value)
    summary.add_row(
        "ACMG 2015",
        "[" + cls2015_color + "]" + result.classification_2015.value + "[/" + cls2015_color + "]"
        + "  (" + result.classification_2015_rules + ")",
    )
    summary.add_row(
        "Bayesian",
        "[" + clsbay_color + "]" + result.classification_bayesian.value + "[/" + clsbay_color + "]"
        + "  (score=" + str(result.bayesian_score) + ")",
    )
    console.print(summary)

    pc = annotation.primary_consequence
    if pc:
        console.print("[bold]Primary Consequence[/bold]")
        cseq_table = Table(box=box.SIMPLE, show_header=False)
        cseq_table.add_column("Field", style="dim")
        cseq_table.add_column("Value")
        cseq_table.add_row("Gene", pc.gene_symbol)
        cseq_table.add_row("Transcript", pc.transcript_id + (" [MANE Select]" if pc.is_mane_select else ""))
        cseq_table.add_row("Consequence", pc.consequence.value)
        if pc.hgvs_c:
            cseq_table.add_row("HGVS c.", pc.hgvs_c)
        if pc.hgvs_p:
            cseq_table.add_row("HGVS p.", pc.hgvs_p)
        console.print(cseq_table)

    console.print("[bold]Criteria[/bold]")
    crit_table = Table(box=box.SIMPLE)
    crit_table.add_column("Criterion")
    crit_table.add_column("Triggered")
    crit_table.add_column("Strength")
    crit_table.add_column("Evidence")

    for r in result.criteria_results:
        triggered_str = "[green]YES[/green]" if (r.triggered and not r.suppressed) else (
            "[yellow]SUPPRESSED[/yellow]" if r.suppressed else "[dim]no[/dim]"
        )
        crit_table.add_row(
            r.criterion.value,
            triggered_str,
            r.strength.value if r.triggered else "",
            r.evidence[:80],
        )
    console.print(crit_table)

    if result.warnings:
        console.print("[bold yellow]Warnings:[/bold yellow]")
        for w in result.warnings:
            console.print("  [yellow]" + w + "[/yellow]")
    console.print()
