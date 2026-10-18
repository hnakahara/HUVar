"""Per-gene BP1 / BP3 applicability from ClinGen VCEP specs (``disease_prevalence.tsv``).

* ``bp1`` / ``bp1_target`` — BP1 ("variant type unlikely to be pathogenic in this
  gene") applicability and the TARGET consequence. Most VCEPs decline BP1; those
  that apply it target ``missense`` (PALB2, APC, BRCA1/2) or — for gain-of-function
  RASopathy genes where loss-of-function is benign — ``truncating``.
* ``bp3`` — BP3 (in-frame indel in a repetitive region) applicability. A VCEP that
  declined BP3 resolves to ``not_applicable`` and suppresses the heuristic.
* ``bp7_phylop`` — the per-gene phyloP policy for BP7. A numeric cutoff defines
  "not highly conserved" (most VCEPs use the global default phyloP100way 2.0;
  some tighten it to 0.1 / 0.2 / 0 — the neurodevelopmental and coagulation
  panels, VHL, RPGR — or 1.5 for the platelet GP genes). The sentinel ``na``
  means the VCEP declared conservation NON-informative (TP53, GALT, BRCA1/2, the
  SCID T-/B-cell genes), so BP7 skips the conservation gate entirely. Blank when
  the gene has no VCEP policy, so the evaluator keeps its global default.
* ``bp7_intronic`` — the BP7 intronic applicability range. ``noncanonical`` (the
  RASopathy / PIK3 panels) admits any intronic position except the canonical
  +/-1,2 sites; blank keeps the Walker deep-intronic (+7/-21) default.
* ``pp3_scope`` — the VCEP PP3 scope read by ``PP3Evaluator``: ``not_applicable``
  (PIK3/mTOR GOF panel), ``no_missense`` (SCID genes: synonymous/intronic splice
  impact only), ``splice_only`` (APC, PALB2, CDH1: no protein-level predictor).
* ``bp4_splice_cutoff`` / ``bp7_splice_cutoff`` — the per-gene SpliceAI
  "no-impact" cutoff for the BP4 splice branch and the BP7 splice gate. A
  variant is splice-benign when ``SpliceAI <= cutoff``. Most VCEPs use the Walker
  default 0.10; some tighten it (the LGMD panels 0.05) or loosen it (RUNX1 0.20,
  GAA/TP53 0.2, HBA2/HBB 0.3, RPGR-BP7 0.2). Blank → the global 0.10 default.

The BP counterpart to :class:`~acmg_classifier.criteria.pp2_genes.PP2Applicability`.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

APPLICABLE = "applicable"
BP4_SCOPES = ("not_applicable", "splice_only", "missense_only", "no_synonymous")
# PP3 scope (pp3_scope): see PP3Evaluator.
PP3_SCOPES = ("not_applicable", "no_missense", "splice_only")
NOT_APPLICABLE = "not_applicable"


class BPApplicability:
    """VCEP BP1/BP3 applicability per gene, loaded once from the TSV. A missing
    file or column degrades to "no VCEP data" (every gene resolves to "")."""

    def __init__(self, tsv_path: Path) -> None:
        self._bp1: dict[str, str] = {}
        self._bp1_target: dict[str, str] = {}
        self._bp1_exclude: dict[str, list[tuple[int, int]]] = {}
        self._bp1_strong: set[str] = set()
        self._bp1_no_splice: set[str] = set()
        self._bp3: dict[str, str] = {}
        self._bp3_regions: dict[str, list[tuple[int, int]]] = {}
        self._bp7_phylop: dict[str, float] = {}
        self._bp7_no_conservation: set[str] = set()
        self._bp7_intronic: dict[str, str] = {}
        self._bp7_intronic_cutoffs: dict[str, tuple[int, int]] = {}
        self._bp4_splice_cutoff: dict[str, float] = {}
        self._bp4_scope: dict[str, str] = {}
        self._pp3_scope: dict[str, str] = {}
        self._bp4_missense_splice_max: dict[str, float] = {}
        self._bp4_exclude_c: dict[str, list[tuple[str, int, int]]] = {}
        self._bp7_splice_cutoff: dict[str, float] = {}
        self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        if not tsv_path.exists():
            return
        with tsv_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gene = (row.get("gene_symbol") or "").strip()
                if not gene:
                    continue
                bp1 = (row.get("bp1") or "").strip().lower()
                if bp1 in (APPLICABLE, NOT_APPLICABLE):
                    self._bp1[gene] = bp1
                target = (row.get("bp1_target") or "").strip().lower()
                if target in ("missense", "truncating", "broad"):
                    self._bp1_target[gene] = target
                excl = _parse_ranges(row.get("bp1_exclude") or "")
                if excl:
                    self._bp1_exclude[gene] = excl
                if (row.get("bp1_strength") or "").strip().lower() == "strong":
                    self._bp1_strong.add(gene)
                if (row.get("bp1_no_splice") or "").strip().lower() == "yes":
                    self._bp1_no_splice.add(gene)
                bp3 = (row.get("bp3") or "").strip().lower()
                if bp3 in (APPLICABLE, NOT_APPLICABLE):
                    self._bp3[gene] = bp3
                regions = _parse_ranges(row.get("bp3_regions") or "")
                if regions:
                    self._bp3_regions[gene] = regions
                raw_phylop = (row.get("bp7_phylop") or "").strip()
                if raw_phylop.lower() == "na":
                    # Conservation declared non-informative → BP7 skips the gate.
                    self._bp7_no_conservation.add(gene)
                elif raw_phylop:
                    try:
                        # "0" is a valid cutoff (RPGR: only accelerated positions
                        # are BP7-eligible), distinct from "" (no VCEP cutoff).
                        self._bp7_phylop[gene] = float(raw_phylop)
                    except ValueError:
                        pass
                intronic = (row.get("bp7_intronic") or "").strip().lower()
                if intronic == "noncanonical":
                    self._bp7_intronic[gene] = intronic
                elif "donor:" in intronic or "acceptor:" in intronic:
                    cutoffs = _parse_intronic_cutoffs(intronic)
                    if cutoffs is not None:
                        self._bp7_intronic[gene] = "parametric"
                        self._bp7_intronic_cutoffs[gene] = cutoffs
                scope = (row.get("bp4_scope") or "").strip().lower()
                if scope in BP4_SCOPES:
                    self._bp4_scope[gene] = scope
                p3 = (row.get("pp3_scope") or "").strip().lower()
                if p3 in PP3_SCOPES:
                    self._pp3_scope[gene] = p3
                excl_c = _parse_c_positions(row.get("bp4_exclude_c") or "")
                if excl_c:
                    self._bp4_exclude_c[gene] = excl_c
                for col, store in (
                    ("bp4_missense_splice_max", self._bp4_missense_splice_max),
                    ("bp4_splice_cutoff", self._bp4_splice_cutoff),
                    ("bp7_splice_cutoff", self._bp7_splice_cutoff),
                ):
                    raw = (row.get(col) or "").strip()
                    if raw:
                        try:
                            store[gene] = float(raw)
                        except ValueError:
                            pass

    def bp4_splice_cutoff(self, gene: str | None) -> float | None:
        """Per-gene SpliceAI no-impact cutoff for the BP4 splice branch, or
        ``None`` (the evaluator then uses the global 0.10 default)."""
        if not gene:
            return None
        return self._bp4_splice_cutoff.get(gene)

    def bp4_scope(self, gene: str | None) -> str:
        """VCEP BP4 scope: ``not_applicable`` / ``splice_only`` (no protein-based
        predictors; BP4 not for missense) / ``missense_only`` (no splice-based BP4)
        / ``no_synonymous`` / "" (no restriction)."""
        if not gene:
            return ""
        return self._bp4_scope.get(gene, "")

    def pp3_scope(self, gene: str | None) -> str:
        """VCEP PP3 scope: ``not_applicable`` (PP3 never applied) /
        ``no_missense`` (PP3 only for synonymous/intronic splice impact; never for
        missense) / ``splice_only`` (no protein-level predictor; a missense meets
        PP3 only through predicted splice impact) / "" (no restriction)."""
        if not gene:
            return ""
        return self._pp3_scope.get(gene, "")

    def bp4_missense_splice_max(self, gene: str | None) -> float | None:
        """Maximum splice-predictor score allowed for a missense BP4 (GP9: 0;
        RUNX1: 0.20), or ``None`` (the default splice-impact gate applies)."""
        if not gene:
            return None
        return self._bp4_missense_splice_max.get(gene)

    def bp4_excluded_position(self, gene: str | None, hgvs_c: str | None) -> bool:
        """True if the variant's MANE c. position is one where the VCEP bars BP4
        (PTEN: variants that may impact the intron 1 donor/acceptor)."""
        rules = self._bp4_exclude_c.get(gene or "")
        if not rules or not hgvs_c:
            return False
        pos = _c_position(hgvs_c)
        if pos is None:
            return False
        base, offset = pos
        for kind, a, b in rules:
            if kind == "+" and offset > 0 and base == a:
                return True
            if kind == "-" and offset < 0 and base == a:
                return True
            if kind == "range" and offset == 0 and a <= base <= b:
                return True
        return False

    def bp7_splice_cutoff(self, gene: str | None) -> float | None:
        """Per-gene SpliceAI no-impact cutoff for the BP7 splice gate, or
        ``None`` (the evaluator then uses the global 0.10 default)."""
        if not gene:
            return None
        return self._bp7_splice_cutoff.get(gene)

    def bp1(self, gene: str | None) -> str:
        """VCEP BP1 status: ``applicable`` / ``not_applicable`` / ""."""
        if not gene:
            return ""
        return self._bp1.get(gene, "")

    def bp1_target(self, gene: str | None) -> str:
        """BP1 target consequence for *gene*: ``missense`` / ``truncating`` / ""."""
        if not gene:
            return ""
        return self._bp1_target.get(gene, "")

    def bp1_excluded(self, gene: str | None, position: int | None) -> bool:
        """True if *position* falls in a region where BP1 must NOT apply (APC's
        β-catenin repeat; the BRCA1/2 clinically-important functional domains)."""
        if not gene or position is None:
            return False
        return any(a <= position <= b for a, b in self._bp1_exclude.get(gene, ()))

    def bp1_is_strong(self, gene: str | None) -> bool:
        """True if BP1 is applied at Strong for *gene* (BRCA1/2)."""
        return bool(gene) and gene in self._bp1_strong

    def bp1_requires_no_splice(self, gene: str | None) -> bool:
        """True if BP1 requires no predicted splice impact (BRCA1/2: SpliceAI<=0.1)."""
        return bool(gene) and gene in self._bp1_no_splice

    def bp3(self, gene: str | None) -> str:
        """VCEP BP3 status: ``applicable`` / ``not_applicable`` / ""."""
        if not gene:
            return ""
        return self._bp3.get(gene, "")

    def bp3_in_region(self, gene: str | None, position: int | None) -> bool | None:
        """For a gene whose BP3 is region-restricted (RPGR ORF15, FOXG1 poly-AA),
        whether *position* is inside an allowed region. ``None`` when the gene has
        no region restriction (the generic repeat heuristic applies)."""
        regions = self._bp3_regions.get(gene or "")
        if not regions:
            return None
        return position is not None and any(a <= position <= b for a, b in regions)

    def bp7_phylop(self, gene: str | None) -> float | None:
        """The VCEP's per-gene phyloP "highly conserved" cutoff for BP7, or
        ``None`` (use the global default). A position with phyloP >= this value
        is highly conserved and blocks BP7."""
        if not gene:
            return None
        return self._bp7_phylop.get(gene)

    def bp7_conservation_na(self, gene: str | None) -> bool:
        """True if the VCEP declared conservation NON-informative for BP7 (TP53,
        GALT, BRCA1/2, the SCID T-/B-cell genes). The evaluator then skips the
        phyloP conservation gate entirely — BP7 rests on splice + distance."""
        return bool(gene) and gene in self._bp7_no_conservation

    def bp7_intronic_mode(self, gene: str | None) -> str:
        """BP7 intronic range mode for *gene*: ``"noncanonical"`` (any intronic
        position except the canonical +/-1,2 sites — RASopathy / PIK3 panels),
        ``"parametric"`` (explicit per-gene donor/acceptor cutoffs — see
        :meth:`bp7_intronic_cutoffs`, e.g. the Cardiomyopathy panel's -4/+7), or
        ``""`` (the Walker deep-intronic +7/-21 default)."""
        if not gene:
            return ""
        return self._bp7_intronic.get(gene, "")

    def bp7_intronic_cutoffs(self, gene: str | None) -> tuple[int, int] | None:
        """The ``(donor_min, acceptor_max)`` BP7 cutoffs for a ``parametric``
        gene — an intronic variant is eligible when ``dist >= donor_min`` or
        ``dist <= acceptor_max`` (e.g. Cardiomyopathy ``(7, -4)``). ``None`` when
        the gene uses a named mode or the default."""
        if not gene:
            return None
        return self._bp7_intronic_cutoffs.get(gene)


def _parse_intronic_cutoffs(raw: str) -> tuple[int, int] | None:
    """Parse a parametric ``bp7_intronic`` cell ("donor:7,acceptor:-4") into
    ``(donor_min, acceptor_max)``. Order-independent; returns ``None`` when both
    cutoffs cannot be parsed. ``donor_min`` defaults to the Walker +7 and
    ``acceptor_max`` to -21 when only one side is given."""
    donor_min: int | None = None
    acceptor_max: int | None = None
    for part in raw.split(","):
        key, _, val = part.strip().partition(":")
        key = key.strip()
        try:
            n = int(val.strip())
        except ValueError:
            continue
        if key == "donor":
            donor_min = n
        elif key == "acceptor":
            acceptor_max = n
    if donor_min is None and acceptor_max is None:
        return None
    return (donor_min if donor_min is not None else 7,
            acceptor_max if acceptor_max is not None else -21)


def _parse_ranges(raw: str) -> list[tuple[int, int]]:
    """Parse ";"-joined "a-b" residue ranges (e.g. "1021-1035;2-101")."""
    out: list[tuple[int, int]] = []
    for part in raw.split(";"):
        part = part.strip()
        if not part or "-" not in part:
            continue
        a, _, b = part.partition("-")
        try:
            out.append((int(a), int(b)))
        except ValueError:
            continue
    return out


_C_POS_RE = re.compile(r"c\.(\d+)([+-]\d+)?")


def _c_position(hgvs_c: str) -> tuple[int, int] | None:
    """(base, intronic offset) of the first position of a MANE c. change:
    'c.79+7A>G' -> (79, 7); 'c.80-3C>T' -> (80, -3); 'c.78C>T' -> (78, 0).
    None for UTR ('c.-12', 'c.*5') or unparsable changes."""
    h = hgvs_c.split(":", 1)[1] if ":" in hgvs_c else hgvs_c
    m = _C_POS_RE.match(h.strip())
    if not m:
        return None
    return int(m.group(1)), int(m.group(2) or 0)


def _parse_c_positions(raw: str) -> list[tuple[str, int, int]]:
    """Parse bp4_exclude_c tokens: 'c.N+' (intronic, donor side of base N),
    'c.N-' (intronic, acceptor side of base N), 'c.A_B' (exonic bases A..B)."""
    out: list[tuple[str, int, int]] = []
    for tok in raw.split(","):
        tok = tok.strip()
        m = re.fullmatch(r"c\.(\d+)([+-])", tok)
        if m:
            out.append((m.group(2), int(m.group(1)), int(m.group(1))))
            continue
        m = re.fullmatch(r"c\.(\d+)(?:_(\d+))?", tok)
        if m:
            a = int(m.group(1))
            out.append(("range", a, int(m.group(2) or a)))
    return out
