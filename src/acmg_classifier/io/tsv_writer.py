"""Write classification results to TSV."""
from __future__ import annotations
import csv
import re
import sys
from pathlib import Path
from typing import Optional

from acmg_classifier.models.classification import ClassificationResult
from acmg_classifier.models.enums import ACMGCriterion
from acmg_classifier.models.variant import VariantRecord

_ALL_CRITERIA = list(ACMGCriterion)

_ANNOTATION_HEADER = [
    "gnomad_af", "gnomad_ac", "gnomad_an",
    "gnomad_faf95_popmax", "gnomad_popmax_af", "gnomad_popmax_pop",
    "gnomad_pli", "gnomad_loeuf",
    "clinvar_variation_id", "clinvar_significance", "clinvar_stars",
    "alphamissense_score", "alphamissense_classification",
    "splice_tool", "splice_score",
    "in_repeat", "repeat_class",
]

_HEADER = (
    [
        "variant_id", "chrom", "pos", "ref", "alt", "filter",
        "transcript", "gene", "hgvs_c", "hgvs_p",
        "classification_2015", "rules_2015", "bayesian_score", "classification_bayesian",
    ]
    + _ANNOTATION_HEADER
    + [c.value for c in _ALL_CRITERIA]
    + [c.value + "_strength" for c in _ALL_CRITERIA]
    + [c.value + "_evidence" for c in _ALL_CRITERIA]
    + ["warnings"]
)


def write_tsv(results: list[ClassificationResult], output_path: Optional[Path]) -> None:
    fh = open(output_path, "w", newline="", encoding="utf-8") if output_path else sys.stdout
    try:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(_HEADER)
        for result in results:
            writer.writerow([_sanitize_cell(c) for c in _build_row(result)])
    finally:
        if output_path:
            fh.close()


_SKIPPED_HEADER = ["chrom", "pos", "ref", "alt", "vcf_id", "filter", "reason"]
_SKIPPED_REASON = "no ALT allele (ALT='.'): not a variant, excluded from annotation/classification"


def write_skipped(records: list[VariantRecord], output_path: Path) -> None:
    """Write VCF records excluded from classification (e.g. ALT='.' no-variant sites).

    These sites carry no alternate allele, so they cannot be ACMG-classified. This
    file documents which input rows were skipped so the user can reconcile counts.
    """
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(_SKIPPED_HEADER)
        for r in records:
            writer.writerow([
                r.chrom,
                str(r.pos),
                r.ref,
                r.alt or ".",
                r.vcf_id or "",
                r.filter or "",
                _SKIPPED_REASON,
            ])


def _sanitize_cell(value: str) -> str:
    """Collapse newlines/tabs/CR to spaces so each variant stays on one TSV row.

    ClinVar-mined evidence (PS3/PP1) can contain embedded newlines.
    """
    if value is None:
        return ""
    return re.sub(r"\s*[\r\n\t]+\s*", " ", str(value)).strip()


def _strip_transcript_prefix(hgvs: Optional[str]) -> str:
    """NM_015658.4:c.898A>G -> c.898A>G"""
    if not hgvs:
        return ""
    return hgvs.split(":")[-1]


def _fmt(value, fmt=None) -> str:
    if value is None:
        return ""
    if fmt:
        return format(value, fmt)
    return str(value)


def _build_annotation_row(result: ClassificationResult) -> list[str]:
    ann = result.annotation

    if ann is None:
        return [""] * len(_ANNOTATION_HEADER)

    g = ann.gnomad
    gnomad_cols = [
        _fmt(g.af, ".6g") if g else "",
        _fmt(g.ac) if g else "",
        _fmt(g.an) if g else "",
        _fmt(g.faf95_popmax, ".6g") if g else "",
        _fmt(g.popmax_af, ".6g") if g else "",
        g.popmax_pop or "" if g else "",
        _fmt(g.pli, ".4f") if g and g.pli is not None else "",
        _fmt(g.loeuf, ".4f") if g and g.loeuf is not None else "",
    ]

    best_cv = max(ann.clinvar_vcf, key=lambda r: r.star_rating, default=None)
    clinvar_cols = [
        best_cv.variation_id or "" if best_cv else "",
        best_cv.clinical_significance if best_cv else "",
        str(best_cv.star_rating) if best_cv else "",
    ]

    am = ann.alphamissense
    am_cols = [
        _fmt(am.score, ".4f") if am and am.score is not None else "",
        am.classification or "" if am else "",
    ]

    sp = ann.splice
    if sp and sp.is_available:
        score = sp.raw_score if sp.raw_score is not None else sp.max_delta
        splice_cols = [sp.tool, _fmt(score, ".4f") if score is not None else ""]
    else:
        splice_cols = ["", ""]

    rep = ann.repeat
    repeat_cols = [
        "1" if rep and rep.in_repeat else "0",
        rep.repeat_class or "" if rep else "",
    ]

    return gnomad_cols + clinvar_cols + am_cols + splice_cols + repeat_cols


def _build_row(result: ClassificationResult) -> list[str]:
    by_criterion = {r.criterion: r for r in result.criteria_results}

    triggered = [
        "1" if (c in by_criterion and by_criterion[c].triggered and not by_criterion[c].suppressed) else "0"
        for c in _ALL_CRITERIA
    ]
    strengths = [
        by_criterion[c].strength.value if c in by_criterion else ""
        for c in _ALL_CRITERIA
    ]
    evidences = [
        by_criterion[c].evidence if c in by_criterion else ""
        for c in _ALL_CRITERIA
    ]

    return (
        [
            result.variant_id,
            result.chrom,
            str(result.pos),
            result.ref,
            result.alt,
            result.filter or "",
            result.transcript_id or "",
            result.gene_symbol or "",
            _strip_transcript_prefix(result.hgvs_c),
            _strip_transcript_prefix(result.hgvs_p),
            result.classification_2015.value,
            result.classification_2015_rules,
            str(result.bayesian_score),
            result.classification_bayesian.value,
        ]
        + _build_annotation_row(result)
        + triggered
        + strengths
        + evidences
        + ["; ".join(result.warnings)]
    )
