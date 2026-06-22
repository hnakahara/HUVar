"""Disease-specific allele-frequency thresholds (Whiffin/Ware 2017).

BA1/BS1 cutoffs are most defensible when derived from the maximum credible
population allele frequency for the specific disorder, rather than flat
catch-alls. We compute the "maximum credible AF" (a.k.a. disease allele
frequency threshold, DAFT) from gene/disease parameters and map it to:

    BS1 = max(maxAF, 0.0005)        # 0.05% floor guards ultra-rare disorders
    BA1 = min(0.05, 10 x maxAF)     # 5% absolute ceiling for stand-alone benign

Parameters come from ``disease_prevalence.tsv``. Direct ``bs1_threshold`` /
``ba1_threshold`` columns, when present, override the computed values (lets a
VCEP-published cutoff be used verbatim). Genes absent from the table — or with
incomplete parameters — fall back to the historical flat defaults.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from acmg_classifier.models.enums import CriterionStrength

# Flat fallbacks (pre-Whiffin behaviour) when no per-gene data is available.
_DEFAULT_BS1 = 0.005
_DEFAULT_BA1 = 0.05

# Mapping of the computed maximum credible AF to BA1/BS1.
_BS1_FLOOR = 0.0005   # 0.05% — below this BS1 must not fire on ultra-rare genes
_BA1_CAP = 0.05       # 5% — a variant this common is stand-alone benign anyway
_BA1_FACTOR = 10      # BA1 sits an order of magnitude above BS1 (DAFT)


def compute_max_credible_af(
    prevalence: Optional[float],
    allelic_het: Optional[float],
    genetic_het: Optional[float],
    penetrance: Optional[float],
    recessive: bool,
) -> Optional[float]:
    """Maximum credible population allele frequency (Whiffin/Ware 2017).

    ``G`` is the causal-genotype frequency implied by the disorder:
        G = prevalence x genetic_het x allelic_het / penetrance
    For a dominant disorder the causal genotype is heterozygous, so under
    Hardy-Weinberg for a rare allele genotype ~= 2*AF → maxAF = G / 2. For a
    recessive disorder the causal genotype is homozygous (genotype ~= AF^2),
    so maxAF = sqrt(G).

    Returns None when any input is missing or out of range, so callers can
    fall back to a flat default rather than trusting a degenerate value.
    """
    if prevalence is None or penetrance is None:
        return None
    if allelic_het is None or genetic_het is None:
        return None
    # Penetrance must be a positive proportion; the heterogeneity factors and
    # prevalence must be positive proportions too. Reject anything outside
    # (0, 1] (prevalence (0,1)) so a typo can't produce a nonsense threshold.
    if not (0.0 < penetrance <= 1.0):
        return None
    if not (0.0 < allelic_het <= 1.0) or not (0.0 < genetic_het <= 1.0):
        return None
    if not (0.0 < prevalence < 1.0):
        return None

    g = prevalence * genetic_het * allelic_het / penetrance
    if recessive:
        return math.sqrt(g)
    return g / 2.0


@dataclass(frozen=True)
class GeneThresholds:
    """Resolved BA1/BS1 allele-frequency cutoffs for one gene.

    ``af_basis`` selects which gnomAD frequency the BA1/BS1 evaluators compare
    against: "" (default) → overall population FAF95; "males" → male (XY)
    allele frequency, for X-linked genes whose VCEP states the cutoff "in
    males" (RPGR, RS1, ABCD1, SLC6A8, OTC); "popmax" → the grpmax/popmax POINT
    estimate, for VCEPs that define BA1/BS1 on the point allele frequency rather
    than the 95%-CI filtering AF (RUNX1, GAA, MYOC, …) — gated by
    ``Config.popmax_af_basis``.

    ``ba1_hom_count`` is the homozygote+hemizygote count at/above which BA1 fires
    independently of frequency, for VCEPs with that OR-clause (SLC6A8, OTC: ≥10).
    """

    ba1: float
    bs1: float
    af_basis: str = ""
    ba1_hom_count: Optional[int] = None
    # Strength the BS1 cutoff fires at — the spec's BS1 tier (e.g. MYO15A/OTOF
    # fire VeryStrong at >=0.3%); defaults to the ACMG BS1 Strong level.
    bs1_strength: CriterionStrength = CriterionStrength.STRONG
    # Comma-separated bare protein changes the VCEP bars from BS1 regardless of
    # frequency (a recurrent disease allele, e.g. MYOC p.Gln368Ter). The BS1
    # evaluator withholds BS1 for a variant matching one of these.
    bs1_exclude: str = ""


_DEFAULT_THRESHOLDS = GeneThresholds(ba1=_DEFAULT_BA1, bs1=_DEFAULT_BS1)

# TSV bs1_strength value → CriterionStrength (blank → the Strong default).
_BS1_STRENGTH = {
    "VeryStrong": CriterionStrength.VERY_STRONG,
    "Strong": CriterionStrength.STRONG,
    "Moderate": CriterionStrength.MODERATE,
    "Supporting": CriterionStrength.SUPPORTING,
}


def _to_float(row: dict[str, str], key: str) -> Optional[float]:
    raw = (row.get(key) or "").strip()
    if not raw:
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _resolve_row(row: dict[str, str]) -> GeneThresholds:
    """Turn one TSV row into BA1/BS1 cutoffs (override → compute → default)."""
    bs1_override = _to_float(row, "bs1_threshold")
    ba1_override = _to_float(row, "ba1_threshold")

    prevalence = _to_float(row, "prevalence")
    penetrance = _to_float(row, "penetrance")
    # Heterogeneity factors default to 1.0 (a single allele/gene explains all
    # cases) — the most permissive assumption, which yields the *highest*
    # maxAF and therefore the most conservative (hardest-to-fire) benign
    # thresholds when those columns are not supplied.
    allelic_het = _to_float(row, "allelic_het")
    genetic_het = _to_float(row, "genetic_het")
    if prevalence is not None and penetrance is not None:
        if allelic_het is None:
            allelic_het = 1.0
        if genetic_het is None:
            genetic_het = 1.0

    inh = (row.get("inheritance") or "").strip().upper()
    recessive = "AR" in inh or "RECESSIVE" in inh

    max_af = compute_max_credible_af(
        prevalence, allelic_het, genetic_het, penetrance, recessive
    )

    if bs1_override is not None:
        bs1 = bs1_override
    elif max_af is not None:
        bs1 = max(max_af, _BS1_FLOOR)
    else:
        bs1 = _DEFAULT_BS1

    if ba1_override is not None:
        ba1 = ba1_override
    elif max_af is not None:
        ba1 = min(_BA1_CAP, _BA1_FACTOR * max_af)
    else:
        ba1 = _DEFAULT_BA1

    af_basis = (row.get("af_basis") or "").strip().lower()
    bs1_strength = _BS1_STRENGTH.get(
        (row.get("bs1_strength") or "").strip(), CriterionStrength.STRONG
    )
    bs1_exclude = (row.get("bs1_exclude") or "").strip()
    raw_hom = (row.get("ba1_hom_count") or "").strip()
    try:
        ba1_hom_count = int(raw_hom) if raw_hom else None
    except ValueError:
        ba1_hom_count = None
    return GeneThresholds(
        ba1=ba1, bs1=bs1, af_basis=af_basis, bs1_strength=bs1_strength,
        bs1_exclude=bs1_exclude, ba1_hom_count=ba1_hom_count,
    )


class DiseaseThresholds:
    """Per-gene BA1/BS1 cutoffs loaded once from ``disease_prevalence.tsv``.

    Shared by the BA1 and BS1 evaluators so both read a single source. A
    missing file is treated as "no per-gene data" (every gene gets the flat
    defaults) rather than a hard error, so minimal setups still work.
    """

    def __init__(self, tsv_path: Path) -> None:
        self._by_gene: dict[str, GeneThresholds] = self._load(tsv_path)

    @staticmethod
    def _load(tsv_path: Path) -> dict[str, GeneThresholds]:
        if not tsv_path.exists():
            return {}
        out: dict[str, GeneThresholds] = {}
        # UTF-8 explicitly: the curated table carries non-ASCII (VCEP names,
        # ≥/µ in notes). Relying on the platform default (cp932 on Windows)
        # would raise UnicodeDecodeError at load time.
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                out[gene] = _resolve_row(row)
        return out

    def get(self, gene: Optional[str]) -> GeneThresholds:
        """BA1/BS1 cutoffs for *gene*, or the flat defaults if unknown."""
        if not gene:
            return _DEFAULT_THRESHOLDS
        return self._by_gene.get(gene, _DEFAULT_THRESHOLDS)
