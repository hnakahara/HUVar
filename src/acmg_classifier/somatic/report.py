"""TSV output and terminal report for somatic oncogenicity results."""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Optional

from acmg_classifier.somatic.models import OncCriterion, Oncogenicity, OncogenicityResult

_CRITERIA = list(OncCriterion)

HEADER = (
    ["variant_id", "chrom", "pos", "ref", "alt", "filter", "transcript", "gene", "hgvs_c",
     "hgvs_p", "consequence", "gene_role", "gene_role_source", "gene_role_note",
     "oncogenicity", "oncogenicity_score", "criteria_met", "hotspot",
     "tumor_type", "amp_tier", "amp_level", "amp_evidence", "amp_note", "amp_candidates",
     "gnomad_af", "gnomad_popmax_af", "gnomad_popmax_pop", "splice_tool", "splice_score"]
    + [c.value for c in _CRITERIA]
    + [c.value + "_evidence" for c in _CRITERIA]
    + ["warnings"]
)


def _clean(v) -> str:
    if v is None:
        return ""
    return str(v).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def build_row(r: OncogenicityResult) -> list:
    ann = r.annotation
    gd = getattr(ann, "gnomad", None)
    sp = getattr(ann, "splice", None)
    by = {c.criterion: c for c in r.criteria}

    def flag(c):
        x = by[c]
        if not x.met:
            return ""
        return f"{x.points:+d}" if not x.suppressed else "suppressed"

    return ([r.variant_id, r.chrom, r.pos, r.ref, r.alt, r.filter, r.transcript_id,
             r.gene_symbol, r.hgvs_c, r.hgvs_p, r.consequence, r.gene_role,
             r.gene_role_source, r.gene_role_note, r.classification.value, r.score,
             r.summary(), r.hotspot,
             r.tumor_type, r.amp_tier, r.amp_level, r.amp_evidence, r.amp_note,
             " | ".join(r.amp_candidates),
             getattr(gd, "af", None), getattr(gd, "popmax_af", None),
             getattr(gd, "popmax_pop", None),
             getattr(sp, "tool", None) if sp is not None and sp.is_available else None,
             getattr(sp, "max_delta", None) if sp is not None and sp.is_available else None]
            + [flag(c) for c in _CRITERIA]
            + [by[c].evidence for c in _CRITERIA]
            + ["; ".join(r.warnings)])


def write_somatic_tsv(results: list[OncogenicityResult], output_path: Optional[Path]) -> None:
    fh = open(output_path, "w", newline="", encoding="utf-8") if output_path else sys.stdout
    try:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        for r in results:
            w.writerow([_clean(x) for x in build_row(r)])
    finally:
        if output_path:
            fh.close()


_COLORS = {
    Oncogenicity.ONCOGENIC: "bold red",
    Oncogenicity.LIKELY_ONCOGENIC: "red",
    Oncogenicity.VUS: "yellow",
    Oncogenicity.LIKELY_BENIGN: "blue",
    Oncogenicity.BENIGN: "bold blue",
}


def print_somatic_report(r: OncogenicityResult) -> None:
    from rich import box
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table

    con = Console()
    con.print()
    con.print(Panel(f"[bold]{r.variant_id}[/bold]",
                    title="Somatic Oncogenicity Report (ClinGen/CGC/VICC SOP)", expand=False))
    t = Table(box=box.SIMPLE, show_header=False)
    t.add_column("Field", style="bold")
    t.add_column("Value")
    t.add_row("Gene", f"{r.gene_symbol or '-'} ({r.gene_role}"
                      f"{', ' + r.gene_role_source if r.gene_role_source else ''})")
    t.add_row("Transcript", r.transcript_id or "-")
    t.add_row("HGVS", f"{r.hgvs_c or '-'}  {(r.hgvs_p or '').split(':')[-1]}")
    t.add_row("Consequence", r.consequence or "-")
    col = _COLORS.get(r.classification, "white")
    t.add_row("Oncogenicity", f"[{col}]{r.classification.value}[/{col}]  (score {r.score:+d})")
    t.add_row("Criteria", r.summary() or "none")
    if r.hotspot:
        t.add_row("Hotspot", r.hotspot)
    if r.amp_tier:
        t.add_row("Tumour type", r.tumor_type)
        lvl = f", level {r.amp_level}" if r.amp_level else ""
        t.add_row("AMP/ASCO/CAP", f"[bold]{r.amp_tier}{lvl}[/bold]"
                  + (f"  ({r.amp_note})" if r.amp_note else ""))
    else:
        t.add_row("AMP/ASCO/CAP", r.amp_note or "not assigned")
    con.print(t)
    if r.amp_evidence:
        et = Table(title="AMP/ASCO/CAP evidence (strongest first)", box=box.SIMPLE_HEAD)
        et.add_column("Evidence")
        for e in r.amp_evidence.split(" | "):
            et.add_row(e)
        con.print(et)
    if r.amp_candidates:
        ct0 = Table(title="Evidence in other tumour types (reference; give --tumor-type)",
                    box=box.SIMPLE_HEAD)
        ct0.add_column("Evidence")
        for e in r.amp_candidates:
            ct0.add_row(e)
        con.print(ct0)

    ct = Table(title="Criteria", box=box.SIMPLE_HEAD)
    ct.add_column("Code", style="bold")
    ct.add_column("Points", justify="right")
    ct.add_column("Evidence")
    for c in r.criteria:
        if c.met and not c.suppressed:
            pts = f"[green]{c.points:+d}[/green]"
        elif c.suppressed:
            pts = "[dim]suppr.[/dim]"
        else:
            pts = "[dim]-[/dim]"
        ct.add_row(c.criterion.value, pts, c.evidence or "")
    con.print(ct)
    for w in r.warnings:
        con.print(f"[yellow]! {w}[/yellow]")
    con.print("[dim]Points: Oncogenic >= 10, Likely oncogenic 6-9, VUS 0-5, Likely benign "
              "-1 to -6, Benign <= -7 (Horak et al. 2022). Tiers: AMP/ASCO/CAP 2017 (Li et al.) "
              "from CIViC and ClinVar somatic clinical impact.[/dim]")
