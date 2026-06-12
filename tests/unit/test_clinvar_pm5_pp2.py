"""Precision-focused tests for the PM5 / PP2 ClinVar eligibility queries.

PM5: the same-codon comparator must be a pathogenic *missense* (not a
truncating change that merely shares the residue number). PP2: the gene-level
eligibility thresholds were tightened to curb over-assignment, including a
benign-rate ceiling on the gnomAD missense-Z rescue branch.
"""
import sqlite3
from pathlib import Path

from acmg_classifier.local_db.clinvar_sqlite import (
    query_same_codon_different_aa,
    query_pp2_eligible,
)

_SCHEMA = """
CREATE TABLE variants (
    variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER
);
"""


def _db(tmp_path: Path, rows: list[tuple]) -> Path:
    """Create a minimal ClinVar sqlite with the given variant rows."""
    p = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, gene_symbol, hgvs_p, "
        "amino_acid_change, codon_position, clinical_significance, "
        "review_status, star_rating) VALUES (?,?,?,?,?,?,?,?)",
        rows,
    )
    con.commit()
    con.close()  # release the writer before the read-only immutable open
    return p


def _row(vid, gene, hgvs_p, aa, codon, sig, stars):
    return (vid, gene, hgvs_p, aa, codon, sig, "criteria provided", stars)


class TestPM5MissenseComparator:
    def test_missense_comparator_is_a_hit(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "GENE", "NM_1:p.Arg175His", "R175H", 175, "Pathogenic", 1),
        ])
        hits = query_same_codon_different_aa(db, "GENE", 175, "NM_1:p.Arg175Gln")
        assert [h.variation_id for h in hits] == ["1"]

    def test_truncating_comparator_excluded(self, tmp_path):
        # A pathogenic nonsense at the same residue must NOT establish PM5.
        db = _db(tmp_path, [
            _row("2", "GENE", "NM_1:p.Arg175Ter", "R175*", 175, "Pathogenic", 2),
        ])
        hits = query_same_codon_different_aa(db, "GENE", 175, "NM_1:p.Arg175Gln")
        assert hits == []

    def test_mixed_keeps_only_missense(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "GENE", "NM_1:p.Arg175His", "R175H", 175, "Pathogenic", 1),
            _row("2", "GENE", "NM_1:p.Arg175Ter", "R175*", 175, "Pathogenic", 2),
        ])
        hits = query_same_codon_different_aa(db, "GENE", 175, "NM_1:p.Arg175Gln")
        assert [h.variation_id for h in hits] == ["1"]

    def test_same_aa_change_excluded_ps1_not_pm5(self, tmp_path):
        # Identical AA change is PS1's job, not PM5 — must be filtered out.
        db = _db(tmp_path, [
            _row("1", "GENE", "NM_1:p.Arg175Gln", "R175Q", 175, "Pathogenic", 2),
        ])
        hits = query_same_codon_different_aa(db, "GENE", 175, "NM_1:p.Arg175Gln")
        assert hits == []


class TestPM5TranscriptCollision:
    """The same codon-proximity guard PS1 uses: a comparator that merely shares
    codon_position but sits far away is a transcript-numbering collision."""

    _SCHEMA2 = (
        "CREATE TABLE variants (variation_id TEXT, chrom TEXT, pos INTEGER, "
        "ref TEXT, alt TEXT, gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, "
        "amino_acid_change TEXT, codon_position INTEGER, "
        "clinical_significance TEXT, review_status TEXT, star_rating INTEGER)"
    )

    def _db2(self, tmp_path, rows):
        p = tmp_path / "cv2.sqlite"
        con = sqlite3.connect(p)
        con.execute(self._SCHEMA2)
        con.executemany(
            "INSERT INTO variants (variation_id, chrom, pos, ref, alt, "
            "gene_symbol, hgvs_p, amino_acid_change, codon_position, "
            "clinical_significance, review_status, star_rating) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            rows,
        )
        con.commit()
        con.close()
        return p

    def test_far_codon_collision_rejected(self, tmp_path):
        db = self._db2(tmp_path, [
            ("c", "10", 87880000, "C", "A", "PTEN", "NP_x:p.Pro38His", "P38H",
             38, "Pathogenic", "criteria provided", 2),
        ])
        hits = query_same_codon_different_aa(
            db, "PTEN", 38, "NP_y:p.Pro38Leu",
            query_chrom="chr10", query_pos=87894058,
        )
        assert hits == []

    def test_same_codon_within_window_found(self, tmp_path):
        db = self._db2(tmp_path, [
            ("c", "10", 87894059, "C", "G", "PTEN", "NP_x:p.Pro38Arg", "P38R",
             38, "Pathogenic", "criteria provided", 2),
        ])
        hits = query_same_codon_different_aa(
            db, "PTEN", 38, "NP_y:p.Pro38Leu",
            query_chrom="chr10", query_pos=87894058,
        )
        assert [h.variation_id for h in hits] == ["c"]

    def test_null_pos_comparator_kept(self, tmp_path):
        db = self._db2(tmp_path, [
            ("c", None, None, None, None, "PTEN", "NP_x:p.Pro38Arg", "P38R",
             38, "Pathogenic", "criteria provided", 2),
        ])
        hits = query_same_codon_different_aa(
            db, "PTEN", 38, "NP_y:p.Pro38Leu",
            query_chrom="chr10", query_pos=87894058,
        )
        assert [h.variation_id for h in hits] == ["c"]


def _pp2_rows(gene, n_path, n_benign):
    rows = []
    for i in range(n_path):
        rows.append(_row(f"P{i}", gene, f"NM_1:p.Val{100 + i}Ile", f"V{100+i}I",
                         100 + i, "Pathogenic", 1))
    for i in range(n_benign):
        rows.append(_row(f"B{i}", gene, f"NM_1:p.Ala{300 + i}Thr", f"A{300+i}T",
                         300 + i, "Benign", 1))
    return rows


class TestPP2Eligibility:
    def test_enough_path_and_clean_is_eligible(self, tmp_path):
        db = _db(tmp_path, _pp2_rows("CLEAN", 10, 0))
        ok, _ = query_pp2_eligible(db, "CLEAN")
        assert ok

    def test_below_min_path_not_eligible(self, tmp_path):
        # 9 P/LP missense < tightened _PP2_MIN_PATH (10).
        db = _db(tmp_path, _pp2_rows("FEWPATH", 9, 0))
        ok, evidence = query_pp2_eligible(db, "FEWPATH")
        assert not ok
        assert "only 9 P/LP missense" in evidence

    def test_benign_rate_above_5pct_needs_z(self, tmp_path):
        # 10 path + 1 benign → frac ~9% > 5% ; no Z → not eligible.
        db = _db(tmp_path, _pp2_rows("MIDBENIGN", 10, 1))
        ok, _ = query_pp2_eligible(db, "MIDBENIGN", mis_z=None)
        assert not ok

    def test_z_rescue_within_ceiling(self, tmp_path):
        # Same gene, frac ~9% <= 15% Z-ceiling, constrained → eligible via Z.
        db = _db(tmp_path, _pp2_rows("ZRESCUE", 10, 1))
        ok, evidence = query_pp2_eligible(db, "ZRESCUE", mis_z=4.0)
        assert ok
        assert "Z-score branch" in evidence

    def test_z_rescue_blocked_by_benign_ceiling(self, tmp_path):
        # 10 path + 3 benign → frac ~23% > 15% ceiling: Z must NOT rescue.
        db = _db(tmp_path, _pp2_rows("DIRTY", 10, 3))
        ok, evidence = query_pp2_eligible(db, "DIRTY", mis_z=5.0)
        assert not ok
        assert "Z-rescue ceiling" in evidence
