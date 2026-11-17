"""Optional exon-aware refinement of the VCEP PVS1 canonical-splice strength.

The flat per-gene splice strength in :mod:`acmg_classifier.pvs1.vcep_pvs1`
(``VERY_STRONG`` for most genes) over-calls canonical ±1,2 splice variants whose
skipped exon is in-frame and/or non-critical — several VCEP decision trees score
those at Strong or Moderate (e.g. DICER1 exon 10 → Strong; exons 5/15/18/22 →
Moderate; CDKL5 exon 17 → Moderate / exon 18 → Strong).

This module lets a reviewer encode those per-(gene, skipped-exon) strengths in a
small TSV that *overrides* the flat default. It is deliberately OPT-IN and
data-driven: with no TSV (or no matching row) the caller keeps the flat default,
so behaviour is unchanged until a reviewer supplies a verified table.

The exon a canonical splice variant skips is derived from VEP's ``intron = n/N``
field per the VCEP convention (also stated verbatim in the ACADVL/GAMT trees):

    * a DONOR (+1/+2) variant in intron *n* skips the exon 5' of it  → exon *n*
    * an ACCEPTOR (-1/-2) variant in intron *n* skips the exon 3' of it → exon *n+1*

Reviewers should populate the override TSV using the coordinate table emitted by
``scripts/build_vcep_pvs1_exons.py`` (which lists, per coding exon, whether
skipping is in-frame and what % of the protein it removes) so the exon numbering
matches the RefSeq transcript VEP annotates against — avoiding the off-by-one
trap from genes with a non-coding exon 1 (HNF4A, CDKL5, MECP2, DICER1, …).

Override TSV columns: ``gene``, ``skipped_exon``, ``strength`` (one of
very_strong / strong / moderate / supporting / na), optional ``note``.
"""
from __future__ import annotations

import csv
from pathlib import Path

from acmg_classifier.models.enums import ConsequenceType, CriterionStrength

_STR_BY_NAME = {
    "very_strong": CriterionStrength.VERY_STRONG,
    "pvs1": CriterionStrength.VERY_STRONG,
    "strong": CriterionStrength.STRONG,
    "pvs1_strong": CriterionStrength.STRONG,
    "moderate": CriterionStrength.MODERATE,
    "pvs1_moderate": CriterionStrength.MODERATE,
    "supporting": CriterionStrength.SUPPORTING,
    "pvs1_supporting": CriterionStrength.SUPPORTING,
    "na": CriterionStrength.NOT_MET,
    "n/a": CriterionStrength.NOT_MET,
    "not_met": CriterionStrength.NOT_MET,
}


def skipped_exon(consequence: ConsequenceType, intron: str | None) -> int | None:
    """Exon predicted to be skipped by a canonical ±1,2 splice variant, from the
    VEP ``intron = n/N`` field. Donor → exon n; acceptor → exon n+1. Returns
    None when the intron index cannot be parsed."""
    if not intron:
        return None
    try:
        n = int(intron.split("/")[0])
    except (ValueError, IndexError):
        return None
    if consequence == ConsequenceType.SPLICE_DONOR:
        return n
    if consequence == ConsequenceType.SPLICE_ACCEPTOR:
        return n + 1
    return None


class SpliceExonOverrides:
    """Per-(gene, skipped-exon) PVS1 splice-strength overrides loaded from a TSV.
    A missing file leaves the table empty, so every lookup misses and the caller
    falls back to the flat per-gene default (behaviour unchanged)."""

    def __init__(self, tsv_path: Path | None) -> None:
        self._table: dict[tuple[str, int], CriterionStrength] = {}
        if tsv_path is not None and tsv_path.exists():
            self._load(tsv_path)

    def _load(self, tsv_path: Path) -> None:
        with tsv_path.open(encoding="utf-8") as fh:
            # Drop leading '#' comment lines so the real column header (not a
            # comment) is what csv.DictReader uses for field names.
            lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
            for row in csv.DictReader(lines, delimiter="\t"):
                gene = (row.get("gene") or "").strip()
                exon_raw = (row.get("skipped_exon") or "").strip()
                name = (row.get("strength") or "").strip().lower()
                if not gene or gene.startswith("#") or not exon_raw:
                    continue
                strength = _STR_BY_NAME.get(name)
                if strength is None:
                    continue
                try:
                    exon = int(exon_raw)
                except ValueError:
                    continue
                self._table[(gene, exon)] = strength

    def __bool__(self) -> bool:
        return bool(self._table)

    def lookup(self, gene: str, exon: int | None) -> CriterionStrength | None:
        if exon is None:
            return None
        return self._table.get((gene, exon))


# ---------------------------------------------------------------------------
# Reading-frame effect of skipping the exon (ClinGen SVI PVS1 tree, splice arm)
# ---------------------------------------------------------------------------

from typing import NamedTuple, Optional


class ExonSkip(NamedTuple):
    """Coding consequence of skipping one exon of a MANE transcript."""
    exon: int                 # skipped exon (VEP numbering, all exons counted)
    n_exons: int
    coding_nt: int            # coding bases in the skipped exon
    in_frame: bool            # coding_nt divisible by 3
    pct_protein: float        # coding_nt / CDS length (0-1)
    has_start: bool           # exon holds the initiation codon
    has_stop: bool            # exon holds the stop codon
    penultimate: bool         # exon is the last-but-one exon of the transcript
    pct_from_exon: float      # CDS fraction from this exon's first coding base to the stop


def exon_skip_effect(consequence: ConsequenceType, intron: str | None,
                     transcript_id: str | None) -> Optional[ExonSkip]:
    """Frame / size of the exon a canonical +/-1,2 variant is predicted to skip,
    computed from the MANE exon structure (``resources/shared/mane_exons.tsv``,
    same numbering as VEP). None when the transcript structure or the intron
    number is unavailable."""
    exon = skipped_exon(consequence, intron)
    if exon is None:
        return None
    from acmg_classifier.pvs1.nmd_predictor import transcript_structure
    st = transcript_structure(transcript_id)
    if st is None or not (1 <= exon <= st.n_exons):
        return None
    start = sum(st.exon_lengths[:exon - 1]) + 1
    end = start + st.exon_lengths[exon - 1] - 1
    lo, hi = max(start, st.cds_start), min(end, st.cds_end)
    coding = hi - lo + 1 if hi >= lo else 0
    cds_len = st.cds_end - st.cds_start + 1
    return ExonSkip(
        exon=exon, n_exons=st.n_exons, coding_nt=coding,
        in_frame=coding > 0 and coding % 3 == 0,
        pct_protein=coding / cds_len if cds_len else 0.0,
        has_start=start <= st.cds_start <= end,
        has_stop=start <= st.cds_end <= end,
        penultimate=exon == st.n_exons - 1,
        pct_from_exon=((st.cds_end - lo + 1) / cds_len) if coding and cds_len else 0.0,
    )


def _ten_pct(frac: float) -> CriterionStrength:
    return CriterionStrength.STRONG if frac > 0.10 else CriterionStrength.MODERATE


def svi_exon_skip_strength(eff: Optional[ExonSkip],
                           inframe_strength: CriterionStrength | None = None,
                           ) -> Optional[tuple[CriterionStrength, str]]:
    """PVS1 strength for a canonical +/-1,2 splice variant by the SVI tree's
    exon-skipping branch (Abou Tayoun et al. 2018, "role of region unknown"):

    * in-frame skip -> PVS1_Strong if the exon encodes >10% of the protein, else
      PVS1_Moderate (``inframe_strength`` replaces this when a VCEP fixes it);
    * skip of the exon holding the stop codon, or out-of-frame skip of the
      penultimate exon (frameshift that reaches the last exon, so NMD is not
      expected) -> 10% rule on the C-terminal part from the skipped exon.

    Returns None — keep the caller's default (frameshift predicted to undergo
    NMD -> PVS1) — for out-of-frame skips of other exons, for exons holding the
    initiation codon (start-loss, a different branch) and when the exon
    structure is unknown."""
    if eff is None or eff.coding_nt == 0 or eff.has_start:
        return None
    pct = f"{eff.pct_protein * 100:.1f}%"
    if eff.has_stop:
        frac = eff.pct_from_exon
        s = _ten_pct(frac)
        return s, (f"skips exon {eff.exon}/{eff.n_exons} holding the stop codon "
                   f"({frac * 100:.1f}% of the protein altered, NMD not expected; "
                   "10% rule)")
    if eff.in_frame:
        if inframe_strength is not None:
            return inframe_strength, (f"in-frame skip of exon {eff.exon}/{eff.n_exons} "
                                      f"({eff.coding_nt} nt, {pct} of the protein)")
        s = _ten_pct(eff.pct_protein)
        return s, (f"in-frame skip of exon {eff.exon}/{eff.n_exons} ({eff.coding_nt} nt, "
                   f"{pct} of the protein; role of region unknown, 10% rule)")
    if eff.penultimate:
        frac = eff.pct_from_exon
        s = _ten_pct(frac)
        return s, (f"out-of-frame skip of penultimate exon {eff.exon}/{eff.n_exons}: "
                   f"frameshift into the last exon, NMD not expected "
                   f"({frac * 100:.1f}% of the protein altered; 10% rule)")
    return None

