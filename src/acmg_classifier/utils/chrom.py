"""Chromosome-name helpers for matching regardless of 'chr' prefix.

Input variants are normalised to a 'chr'-prefixed form by VariantRecord, but the
underlying databases may be built with either convention ('chr1' or '1'). These
helpers let position-based DB queries match both forms so lookups never fail
silently on a naming mismatch.
"""
from __future__ import annotations


def strip_chr(chrom: str) -> str:
    """'chr1' -> '1'; '1' -> '1'."""
    return chrom[3:] if chrom.startswith("chr") else chrom


def chrom_candidates(chrom: str) -> list[str]:
    """Return both the bare and 'chr'-prefixed forms, e.g. ['1', 'chr1'].

    Use in a SQL ``WHERE chrom IN (?, ?)`` clause so the query matches whichever
    convention the database was built with.
    """
    bare = strip_chr(chrom)
    return [bare, f"chr{bare}"]
