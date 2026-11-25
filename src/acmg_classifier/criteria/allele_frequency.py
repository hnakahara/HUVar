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
import re

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
    # gnomAD sample subset the BA1/BS1 cutoff is judged on: "" (default) → the
    # overall/full release; "non_cancer" → the non-cancer subset's popmax FAF95
    # (ENIGMA BRCA1/2). Sourced from an explicit ``af_subset`` column, falling
    # back to the row's ``pm2_subset`` so genes already flagged non-cancer for
    # PM2 (BRCA1/2) get the same basis for BA1/BS1 without a TSV edit. The BA1/BS1
    # evaluators degrade to the overall FAF95 when the non-cancer value is absent.
    af_subset: str = ""
    ba1_hom_count: Optional[int] = None
    # Strength the BS1 cutoff fires at — the spec's BS1 tier (e.g. MYO15A/OTOF
    # fire VeryStrong at >=0.3%); defaults to the ACMG BS1 Strong level.
    bs1_strength: CriterionStrength = CriterionStrength.STRONG
    # Lower cutoff of the VCEP's BS1_Supporting tier (``bs1_supporting_threshold``;
    # e.g. PTEN 0.0000043, BRCA1/2 0.00002, hearing-loss AR genes 0.0007): a
    # frequency at or above it but below ``bs1`` gives BS1 at Supporting.
    bs1_supporting: Optional[float] = None
    # Comma-separated bare protein changes the VCEP bars from BS1 regardless of
    # frequency (a recurrent disease allele, e.g. MYOC p.Gln368Ter). The BS1
    # evaluator withholds BS1 for a variant matching one of these.
    bs1_exclude: str = ""
    # Same format for BA1 (LGMD "BS1 exception" variants are barred from both).
    ba1_exclude: str = ""


# Point-estimate BA1 / BS1 rules that also need enough observations in the
# population the frequency comes from (general "VCEP Specifications" text):
#   CDH1 "subpopulation must >= 2,000 alleles and have a minimum of five variant
#   alleles present"; MYOC / CYP1B1 "present in >= 5 alleles in any ... population
#   dataset of at least 2,000 observed alleles"; RUNX1 ">= 2,000 alleles tested and
#   variant present in >= 5 alleles"; F8 / F9 "subpopulation must have a minimum of
#   five variant alleles present".   gene -> (min allele count, min allele number)
MIN_GRPMAX_ALLELES: dict[str, tuple[int, int]] = {
    "CDH1": (5, 2000), "MYOC": (5, 2000), "CYP1B1": (5, 2000), "RUNX1": (5, 2000),
    "F8": (5, 0), "F9": (5, 0),
}


def min_alleles_block(gene: str | None, gd) -> str:
    """Reason BA1 / BS1 must be withheld because the GrpMax population carries too
    few variant alleles / too few alleles for the gene's VCEP, or "". Unknown
    counts (a gnomAD DB without the grpmax AC / AN columns) do not block."""
    rule = MIN_GRPMAX_ALLELES.get(gene or "")
    if rule is None or gd is None:
        return ""
    min_ac, min_an = rule
    ac = gd.ac_grpmax if gd.ac_grpmax is not None else gd.ac
    if ac is not None and ac < min_ac:
        return f"{gene} VCEP: only {ac} variant allele(s) (>= {min_ac} required)"
    if min_an and gd.an_grpmax is not None and gd.an_grpmax < min_an:
        return (f"{gene} VCEP: only {gd.an_grpmax} alleles in the GrpMax population "
                f"(>= {min_an} required)")
    return ""


# ENIGMA BRCA1 / BRCA2 BA1 and BS1 (general text): "Do not apply if read depth
# <20". gene -> minimum gnomAD mean read depth around the variant.
AF_MIN_DEPTH: dict[str, float] = {"BRCA1": 20.0, "BRCA2": 20.0}


def depth_block(coverage, gene: str | None, variant) -> str:
    """Reason BA1 / BS1 must be withheld because gnomAD's mean read depth at the
    variant is below the gene's VCEP minimum, or "" (also when the coverage DB is
    unavailable or has no value there: the gate is skipped, as for PM2)."""
    need = AF_MIN_DEPTH.get(gene or "")
    if need is None or coverage is None or not getattr(coverage, "available", False):
        return ""
    end = variant.pos + max(0, len(variant.ref) - 1)
    dp = coverage.mean_depth(variant.chrom, variant.pos, end)
    if dp is None or dp >= need:
        return ""
    return f"{gene} VCEP: gnomAD mean read depth {dp:.0f} < {need:g}"


def filtering_af(ac: int | None, an: int | None, confidence: float = 0.9999) -> float | None:
    """Filtering allele frequency at *confidence* (Whiffin et al. 2017): the
    highest true allele frequency at which observing ``ac`` or more alleles out of
    ``an`` still has probability >= 1 - *confidence* (Poisson). ``None`` when the
    counts are unknown; 0.0 for ``ac`` = 0. gnomAD's FAF95 is the same quantity at
    95%; the CDH1 and F8 / F9 VCEPs state BA1 / BS1 at "99.99% CI"."""
    import math
    if ac is None or an is None or an <= 0 or ac <= 0:
        return None if (ac is None or an is None or an <= 0) else 0.0
    alpha = 1.0 - confidence

    def p_ge(lam: float) -> float:          # P(X >= ac | lambda)
        term, cdf = math.exp(-lam), 0.0
        for k in range(ac):
            if k:
                term *= lam / k
            cdf += term
        return max(0.0, 1.0 - cdf)

    lo, hi = 0.0, float(ac)
    while p_ge(hi) < alpha:                  # p_ge grows with lambda
        lo, hi = hi, hi * 2 + 1
    # largest lambda with P(X >= ac) <= alpha
    for _ in range(80):
        mid = (lo + hi) / 2
        if p_ge(mid) <= alpha:
            lo = mid
        else:
            hi = mid
    return lo / an


def variant_listed(listing: str, hgvs_c: str | None, hgvs_p: str | None) -> str:
    """Return the entry of a comma-separated VCEP variant list (bs1_exclude /
    ba1_exclude) that matches this variant, or "". Entries are protein changes
    ("p.Gln368Ter") or MANE transcript changes ("c.2697+1G>A", "c.191dup");
    a transcript prefix and the inserted/deleted bases of a dup/del are ignored,
    and "*" equals "Ter"."""
    if not listing:
        return ""
    c = _norm_c(hgvs_c)
    p = _norm_p(hgvs_p)
    for tok in listing.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok.startswith("c.") and c and _norm_c(tok) == c:
            return tok
        if tok.startswith("p.") and p and _norm_p(tok) == p:
            return tok
    return ""


_P_RE = re.compile(r"p\.\(?([A-Za-z]{3}\d+(?:[A-Za-z]{3}|Ter|\*|=)(?:fs(?:Ter|\*)?\d*)?)\)?")


def _norm_p(h: str | None) -> str:
    if not h:
        return ""
    m = _P_RE.search(h)
    return m.group(1).replace("*", "Ter").lower() if m else ""


def _norm_c(h: str | None) -> str:
    if not h:
        return ""
    h = h.split(":", 1)[1] if ":" in h else h
    h = h.strip()
    # c.191dupA -> c.191dup ; c.100_102delAAA -> c.100_102del
    h = re.sub(r"(dup|del)[ACGTN]+$", r"\1", h)
    return h.lower()


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
    # Explicit af_subset wins; otherwise inherit the PM2 subset so BRCA1/2 (which
    # carry pm2_subset=non_cancer) get non-cancer BA1/BS1 without a TSV change.
    af_subset = (row.get("af_subset") or "").strip().lower()
    if not af_subset:
        af_subset = (row.get("pm2_subset") or "").strip().lower()
    if af_subset != "non_cancer":
        af_subset = ""
    bs1_strength = _BS1_STRENGTH.get(
        (row.get("bs1_strength") or "").strip(), CriterionStrength.STRONG
    )
    bs1_exclude = (row.get("bs1_exclude") or "").strip()
    ba1_exclude = (row.get("ba1_exclude") or "").strip()
    raw_hom = (row.get("ba1_hom_count") or "").strip()
    try:
        ba1_hom_count = int(raw_hom) if raw_hom else None
    except ValueError:
        ba1_hom_count = None
    try:
        bs1_supporting = float((row.get("bs1_supporting_threshold") or "").strip())
    except ValueError:
        bs1_supporting = None
    return GeneThresholds(
        ba1=ba1, bs1=bs1, af_basis=af_basis, af_subset=af_subset,
        bs1_supporting=bs1_supporting,
        bs1_strength=bs1_strength, bs1_exclude=bs1_exclude, ba1_exclude=ba1_exclude,
        ba1_hom_count=ba1_hom_count,
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
