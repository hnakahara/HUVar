"""SQLite query layer for ClinVar PS1/PM5/PS4 amino-acid-level lookups."""
from __future__ import annotations
import re
import sqlite3
from pathlib import Path
from typing import Optional

import structlog

from acmg_classifier.models.annotation import ClinVarRecord
from acmg_classifier.utils.chrom import chrom_candidates, strip_chr

log = structlog.get_logger()

_P_LP = frozenset({
    "Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic",
})


def _protein_change_only(hgvs_p: Optional[str]) -> Optional[str]:
    """Strip transcript prefix: 'NP_000298.6:p.Gly341Ala' -> 'p.Gly341Ala'."""
    if not hgvs_p:
        return None
    return hgvs_p.split(":", 1)[-1]


def query_same_aa_change(
    db_path: Path,
    gene_symbol: str,
    hgvs_p: Optional[str],
    exclude_chrom: Optional[str] = None,
    exclude_pos: Optional[int] = None,
    exclude_ref: Optional[str] = None,
    exclude_alt: Optional[str] = None,
    min_stars: int = 1,
) -> list[ClinVarRecord]:
    """
    Return ClinVar records with same AA change (PS1).

    Excludes rows matching exclude_* (PS1 requires a *different* nucleotide change
    producing the same AA change; otherwise the hit is the variant itself).
    """
    p_change = _protein_change_only(hgvs_p)
    if not db_path.exists() or not p_change:
        return []

    excl_chrom = strip_chr(exclude_chrom) if exclude_chrom else exclude_chrom

    try:
        con = sqlite3.connect(str(db_path))
        rows = con.execute(
            """
            SELECT variation_id, clinical_significance, review_status, star_rating,
                   gene_symbol, hgvs_c, hgvs_p, amino_acid_change, chrom, pos, ref, alt
            FROM variants
            WHERE gene_symbol = ?
              AND hgvs_p LIKE ?
              AND star_rating >= ?
              AND clinical_significance IN (
                  'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'
              )
            """,
            (gene_symbol, "%:" + p_change, min_stars),
        ).fetchall()
        con.close()
    except Exception as exc:
        log.error("clinvar_sqlite_error", op="same_aa", error=str(exc))
        return []

    results: list[ClinVarRecord] = []
    # Self-match exclusion: PS1 requires a DIFFERENT nucleotide change
    # producing the same AA. If a ClinVar row matches the input variant by
    # CHROM:POS:REF:ALT it is the variant itself and would let PS1 "vote
    # for itself". We strip the chr prefix on both sides because ClinVar
    # historically stored chromosomes both with and without the prefix.
    for r in rows:
        if (excl_chrom is not None and exclude_pos is not None
                and exclude_ref is not None and exclude_alt is not None
                and strip_chr(str(r[8])) == excl_chrom and r[9] == exclude_pos
                and r[10] == exclude_ref and r[11] == exclude_alt):
            continue
        results.append(ClinVarRecord(
            variation_id=str(r[0]),
            clinical_significance=r[1],
            review_status=r[2],
            star_rating=r[3],
            gene_symbol=r[4],
            hgvs_c=r[5],
            hgvs_p=r[6],
            amino_acid_change=r[7],
        ))
    return results


def query_same_codon_different_aa(
    db_path: Path,
    gene_symbol: str,
    codon_position: Optional[int],
    hgvs_p: Optional[str],
    min_stars: int = 1,
) -> list[ClinVarRecord]:
    """Return ClinVar records at same codon with different AA change (PM5)."""
    p_change = _protein_change_only(hgvs_p)
    if not db_path.exists() or codon_position is None:
        return []
    try:
        con = sqlite3.connect(str(db_path))
        rows = con.execute(
            """
            SELECT variation_id, clinical_significance, review_status, star_rating,
                   gene_symbol, hgvs_c, hgvs_p, amino_acid_change
            FROM variants
            WHERE gene_symbol = ?
              AND codon_position = ?
              AND (hgvs_p IS NULL OR hgvs_p NOT LIKE ?)
              AND star_rating >= ?
              AND clinical_significance IN (
                  'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'
              )
            """,
            (gene_symbol, codon_position, "%:" + (p_change or ""), min_stars),
        ).fetchall()
        con.close()
    except Exception as exc:
        log.error("clinvar_sqlite_error", op="same_codon", error=str(exc))
        return []

    return [
        ClinVarRecord(
            variation_id=str(r[0]),
            clinical_significance=r[1],
            review_status=r[2],
            star_rating=r[3],
            gene_symbol=r[4],
            hgvs_c=r[5],
            hgvs_p=r[6],
            amino_acid_change=r[7],
        )
        for r in rows
    ]


def _sum_column(db_path: Path, column: str, chrom: str, pos: int, ref: str, alt: str) -> int:
    """Return SUM(column) across SCV rows for a variant; 0 if column/DB missing.

    Shared helper for the three text-mined evidence counters (PS3, PS4, PP1).
    Each counter column is populated by the clinvar_builder when scanning
    SCV free text — see local_db/clinvar_builder. Returning 0 on missing
    column (OperationalError) lets us read older databases built before a
    given counter was added, without forcing a full rebuild."""
    if not db_path.exists():
        return 0
    # chrom_candidates returns both "1" and "chr1" forms because the ClinVar
    # build pipeline sometimes loaded one or the other depending on source
    # file. Searching both forms avoids relying on a normalised schema.
    c1, c2 = chrom_candidates(chrom)
    try:
        con = sqlite3.connect(str(db_path))
        try:
            row = con.execute(
                f"""
                SELECT COALESCE(SUM({column}), 0)
                FROM variants
                WHERE chrom IN (?, ?) AND pos = ? AND ref = ? AND alt = ?
                """,
                (c1, c2, pos, ref, alt),
            ).fetchone()
            return int(row[0]) if row and row[0] is not None else 0
        finally:
            con.close()
    except sqlite3.OperationalError:
        # Old ClinVar build without this column — backward-compatible no-op.
        return 0
    except Exception as exc:
        log.error("clinvar_sqlite_error", op=column, error=str(exc))
        return 0


def query_affected_cases(db_path: Path, chrom: str, pos: int, ref: str, alt: str) -> int:
    """Total SCV submissions reporting this variant in affected cases (PS4)."""
    return _sum_column(db_path, "affected_cases", chrom, pos, ref, alt)


def query_functional_evidence(db_path: Path, chrom: str, pos: int, ref: str, alt: str) -> int:
    """Total SCV submissions describing damaging functional studies (PS3)."""
    return _sum_column(db_path, "functional_evidence", chrom, pos, ref, alt)


def query_segregation_evidence(db_path: Path, chrom: str, pos: int, ref: str, alt: str) -> int:
    """Total SCV submissions describing cosegregation with disease (PP1)."""
    return _sum_column(db_path, "segregation_evidence", chrom, pos, ref, alt)


def query_hotspot_cluster(
    db_path: Path,
    gene_symbol: str,
    protein_position: Optional[int],
    window: int = 25,
    min_path_variants: int = 3,
) -> tuple[bool, str]:
    """
    Return (is_hotspot, evidence) by counting P/LP variants in a +/-window aa window (PM1).

    A cluster of >=3 P/LP variants within 25 amino acids is considered a hotspot.
    """
    if not db_path.exists() or protein_position is None:
        return False, "No protein position or DB unavailable"
    try:
        con = sqlite3.connect(str(db_path))
        count = con.execute(
            """
            SELECT COUNT(DISTINCT amino_acid_change)
            FROM variants
            WHERE gene_symbol = ?
              AND codon_position BETWEEN ? AND ?
              AND star_rating >= 1
              AND clinical_significance IN (
                  'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'
              )
            """,
            (gene_symbol, protein_position - window, protein_position + window),
        ).fetchone()[0]
        # PM1 requires the region to be "without benign variation". If any B/LB
        # variant (>=1 star) lies in the same window, PM1 is not applicable.
        benign_count = con.execute(
            """
            SELECT COUNT(DISTINCT amino_acid_change)
            FROM variants
            WHERE gene_symbol = ?
              AND codon_position BETWEEN ? AND ?
              AND star_rating >= 1
              AND clinical_significance IN (
                  'Benign', 'Likely benign', 'Benign/Likely benign'
              )
            """,
            (gene_symbol, protein_position - window, protein_position + window),
        ).fetchone()[0]
        con.close()
    except Exception as exc:
        log.error("clinvar_sqlite_error", op="hotspot", error=str(exc))
        return False, str(exc)

    if benign_count > 0:
        return False, (
            "Benign variation present within " + str(window) + " aa of position "
            + str(protein_position) + " in " + gene_symbol + " ("
            + str(benign_count) + " B/LB variant(s)) — PM1 not applicable"
        )

    if count >= min_path_variants:
        return True, (
            "Hotspot cluster: " + str(count) + " P/LP variants within "
            + str(window) + " aa of position " + str(protein_position)
            + " in " + gene_symbol
        )
    return False, "No hotspot cluster (" + str(count) + " P/LP variants in window)"


# --- PP2: missense is a common disease mechanism with low benign missense rate ---
_PP2_PATH = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")
_PP2_BENIGN = ("Benign", "Likely benign", "Benign/Likely benign")
_PP2_MIN_PATH = 5            # missense must be a recurrent pathogenic mechanism
_PP2_MAX_BENIGN_FRAC = 0.10  # gene must have a low rate of benign missense
_PP2_MIN_MIS_Z = 3.09        # gnomAD missense Z-score qualifying a constrained gene

# Match a missense protein change (p.Val377Ile); the trailing AA must not be
# Ter (stop). We use the 3-letter HGVS form because ClinVar normalises to it;
# the regex deliberately rejects synonymous (=) and frameshift (fs) syntax.
_MISSENSE_RE = re.compile(r"p\.[A-Z][a-z]{2}\d+([A-Z][a-z]{2})")


def _is_missense_p(hgvs_p: Optional[str]) -> bool:
    """Recognise true missense (not stop-gain, frameshift, or synonymous).
    Used to filter ClinVar P/LP rows so PP2 / PVS1-cap counts include only
    missense changes — counting truncating variants here would mix
    mechanisms and break the missense-dominant heuristic."""
    if not hgvs_p:
        return False
    m = _MISSENSE_RE.search(hgvs_p)
    return bool(m) and m.group(1) != "Ter"


def query_pp2_eligible(
    db_path: Path,
    gene_symbol: Optional[str],
    mis_z: float | None = None,
    min_stars: int = 1,
) -> tuple[bool, str]:
    """PP2 eligibility (ClinVar + gnomAD missense constraint).

    The gene qualifies when missense is a recurrent pathogenic mechanism
    (>= _PP2_MIN_PATH P/LP missense, >=1 star) AND at least one of:
      (a) ClinVar benign missense rate <= _PP2_MAX_BENIGN_FRAC, OR
      (b) gnomAD missense Z-score >= _PP2_MIN_MIS_Z (constrained gene).

    The Z-score branch matches Franklin's PP2 logic and rescues genes with a
    slightly elevated benign-missense rate that are nonetheless strongly
    constrained against missense (e.g. clean tumour-suppressor / kinase genes).
    """
    if not db_path.exists() or not gene_symbol:
        return False, "No gene or DB unavailable"
    try:
        con = sqlite3.connect(str(db_path))
        rows = con.execute(
            """
            SELECT hgvs_p, clinical_significance
            FROM variants
            WHERE gene_symbol = ?
              AND star_rating >= ?
              AND clinical_significance IN (?,?,?,?,?,?)
            """,
            (gene_symbol, min_stars, *_PP2_PATH, *_PP2_BENIGN),
        ).fetchall()
        con.close()
    except Exception as exc:
        log.error("clinvar_sqlite_error", op="pp2", error=str(exc))
        return False, str(exc)

    path = sum(1 for hp, sig in rows if sig in _PP2_PATH and _is_missense_p(hp))
    benign = sum(1 for hp, sig in rows if sig in _PP2_BENIGN and _is_missense_p(hp))
    total = path + benign

    if path < _PP2_MIN_PATH:
        return False, (
            f"{gene_symbol}: only {path} P/LP missense (<{_PP2_MIN_PATH}) — PP2 not applicable"
        )

    frac = benign / total if total else 0.0
    benign_ok = frac <= _PP2_MAX_BENIGN_FRAC
    z_ok = mis_z is not None and mis_z >= _PP2_MIN_MIS_Z

    if benign_ok:
        return True, (
            f"{gene_symbol}: {path} P/LP missense, benign missense rate {frac:.0%} "
            f"({benign}/{total}) — missense is a common disease mechanism"
        )
    if z_ok:
        return True, (
            f"{gene_symbol}: {path} P/LP missense, missense Z={mis_z:.2f} "
            f">= {_PP2_MIN_MIS_Z} (constrained against missense; benign rate {frac:.0%} "
            f"{benign}/{total} is permitted via Z-score branch)"
        )
    # Neither branch qualifies — report both negative signals so it is clear why.
    z_note = f", missense Z={mis_z:.2f} < {_PP2_MIN_MIS_Z}" if mis_z is not None else ", missense Z unavailable"
    return False, (
        f"{gene_symbol}: benign missense rate {frac:.0%} > {_PP2_MAX_BENIGN_FRAC:.0%} "
        f"({benign}/{total}){z_note} — PP2 not applicable"
    )


# --- PVS1: is loss-of-function an established disease mechanism for the gene? ---
def query_pathogenic_null_count(
    db_path: Path,
    gene_symbol: Optional[str],
    min_stars: int = 1,
) -> int:
    """Count P/LP null (nonsense/frameshift) variants reported in ClinVar for a gene.

    Used by PVS1's LoF-mechanism gate (cf. Franklin's "pathogenic null variants
    reported"). Counts >=min_stars Pathogenic/Likely-pathogenic records whose
    protein change is a stop-gain or frameshift (HGVS p. contains Ter / fs / *).
    Canonical-splice nulls without a protein HGVS are not counted, but
    nonsense+frameshift counts are sufficient to establish the signal.
    """
    if not db_path.exists() or not gene_symbol:
        return 0
    try:
        con = sqlite3.connect(str(db_path))
        try:
            row = con.execute(
                """
                SELECT COUNT(*)
                FROM variants
                WHERE gene_symbol = ?
                  AND star_rating >= ?
                  AND clinical_significance IN (?,?,?)
                  AND (hgvs_p LIKE '%Ter%' OR hgvs_p LIKE '%fs%' OR hgvs_p LIKE '%*%')
                """,
                (gene_symbol, min_stars, *_PP2_PATH),
            ).fetchone()
            return int(row[0]) if row and row[0] is not None else 0
        finally:
            con.close()
    except Exception as exc:
        log.error("clinvar_sqlite_error", op="pvs1_null", error=str(exc))
        return 0


def query_pathogenic_missense_count(
    db_path: Path,
    gene_symbol: Optional[str],
    min_stars: int = 1,
) -> int:
    """Count P/LP missense variants reported in ClinVar for a gene.

    Used by PVS1's strength cap (ClinGen SVI): when a gene has few P/LP null
    variants but many P/LP missense, the disease mechanism is missense-dominant
    and a new null variant may not act via the same mechanism, so PVS1 strength
    is limited to Moderate.
    """
    if not db_path.exists() or not gene_symbol:
        return 0
    try:
        con = sqlite3.connect(str(db_path))
        try:
            rows = con.execute(
                """
                SELECT hgvs_p FROM variants
                WHERE gene_symbol = ?
                  AND star_rating >= ?
                  AND clinical_significance IN (?,?,?)
                """,
                (gene_symbol, min_stars, *_PP2_PATH),
            ).fetchall()
        finally:
            con.close()
    except Exception as exc:
        log.error("clinvar_sqlite_error", op="pvs1_miss", error=str(exc))
        return 0
    return sum(1 for (hp,) in rows if _is_missense_p(hp))
