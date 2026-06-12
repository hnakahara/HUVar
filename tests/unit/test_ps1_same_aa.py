"""PS1 same-amino-acid ClinVar lookup (query_same_aa_change).

Regression tests for the false negatives where a same-AA pathogenic sibling
exists in ClinVar but was not found, because the old query matched a
``hgvs_p LIKE '%:p.Gly175Arg'`` substring:
  * ClinVar stores predicted protein changes parenthesised (``:p.(Gly175Arg)``),
    so the un-parenthesised LIKE pattern never matched them; and
  * when the stored hgvs_p carried a non-MANE transcript's numbering the
    residue number disagreed with VEP's MANE annotation.
The query now matches the MANE-anchored ``amino_acid_change`` column ('G175R').
"""
import sqlite3
from pathlib import Path

from acmg_classifier.local_db.clinvar_sqlite import query_same_aa_change

_SCHEMA = """
CREATE TABLE variants (
    variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER
);
"""


def _db(tmp_path: Path, rows: list[tuple]) -> Path:
    p = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, "
        "hgvs_p, amino_acid_change, codon_position, clinical_significance, "
        "review_status, star_rating) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        rows,
    )
    con.commit()
    con.close()
    return p


def _row(vid, chrom, pos, ref, alt, gene, hgvs_p, aa, codon, sig, stars):
    return (vid, chrom, pos, ref, alt, gene, hgvs_p, aa, codon, sig,
            "criteria provided, single submitter", stars)


class TestQuerySameAaChange:
    def test_parenthesised_comparator_is_found(self, tmp_path):
        # The comparator's hgvs_p is in ClinVar's predicted (parenthesised)
        # form, which the old LIKE pattern could never match. Matching on
        # amino_acid_change finds it.
        db = _db(tmp_path, [
            _row("1", "7", 44150026, "C", "G", "GCK",
                 "NP_000153.1:p.(Gly175Arg)", "G175R", 175, "Pathogenic", 1),
        ])
        hits = query_same_aa_change(
            db, "GCK", "NP_000153.1:p.Gly175Arg",
            exclude_chrom="chr7", exclude_pos=44150025,
            exclude_ref="C", exclude_alt="T",
        )
        assert [h.variation_id for h in hits] == ["1"]

    def test_different_nucleotide_sibling_found(self, tmp_path):
        # MECP2 p.Leu136Phe arises from two different alts at the same position;
        # each is the other's PS1 sibling. The C>A sibling supports PS1 for C>G.
        db = _db(tmp_path, [
            _row("s", "X", 154032212, "C", "A", "MECP2",
                 "NP_004983.1:p.Leu136Phe", "L136F", 136, "Pathogenic", 2),
        ])
        hits = query_same_aa_change(
            db, "MECP2", "ENSP0:p.Leu136Phe",
            exclude_chrom="chrX", exclude_pos=154032212,
            exclude_ref="C", exclude_alt="G",
        )
        assert [h.variation_id for h in hits] == ["s"]

    def test_self_excluded(self, tmp_path):
        # Only the variant's own record exists (no different-nucleotide sibling,
        # e.g. RPE65 p.Gly484Val). Self-exclusion by coordinate -> no PS1 hit.
        db = _db(tmp_path, [
            _row("self", "1", 68429927, "C", "A", "RPE65",
                 "NP_000320.1:p.Gly484Val", "G484V", 484, "Pathogenic", 2),
        ])
        hits = query_same_aa_change(
            db, "RPE65", "NP_000320.1:p.Gly484Val",
            exclude_chrom="chr1", exclude_pos=68429927,
            exclude_ref="C", exclude_alt="A",
        )
        assert hits == []

    def test_different_aa_not_matched(self, tmp_path):
        # A pathogenic G175Gln at the same codon must NOT satisfy PS1 for G175Arg
        # (that is PM5's job). This also guards the old code[:1] collapse, which
        # would have made Arg and Gln share a key.
        db = _db(tmp_path, [
            _row("x", "7", 44150025, "C", "A", "GCK",
                 "NP_000153.1:p.Gly175Gln", "G175Q", 175, "Pathogenic", 2),
        ])
        hits = query_same_aa_change(db, "GCK", "NP_000153.1:p.Gly175Arg")
        assert hits == []

    def test_benign_and_low_star_excluded(self, tmp_path):
        db = _db(tmp_path, [
            _row("vus", "7", 44150030, "G", "A", "GCK",
                 "NP_000153.1:p.Gly175Arg", "G175R", 175,
                 "Uncertain significance", 2),
            _row("zero", "7", 44150031, "G", "T", "GCK",
                 "NP_000153.1:p.Gly175Arg", "G175R", 175, "Pathogenic", 0),
        ])
        hits = query_same_aa_change(db, "GCK", "NP_000153.1:p.Gly175Arg")
        assert hits == []

    def test_transcript_numbering_collision_rejected(self, tmp_path):
        # The reported PTEN false positive: a comparator that shares the
        # amino_acid_change STRING 'P38L' but sits far from the candidate is the
        # SAME residue number on a DIFFERENT transcript (MANE NM_000314.8 codon 38
        # vs the long isoform NM_001304718, where MANE codon 38 == isoform codon
        # 211). It is a different residue and must NOT satisfy PS1.
        db = _db(tmp_path, [
            _row("collision", "10", 87880000, "C", "T", "PTEN",
                 "NP_001291647.1:p.Pro38Leu", "P38L", 38, "Pathogenic", 2),
        ])
        hits = query_same_aa_change(
            db, "PTEN", "NP_000305.3:p.Pro38Leu",
            exclude_chrom="chr10", exclude_pos=87894058,
            exclude_ref="C", exclude_alt="T",
        )
        assert hits == []

    def test_same_codon_sibling_within_window_found(self, tmp_path):
        # A genuine different-nucleotide sibling lies in the same codon (<=2 bp),
        # so the proximity guard keeps it.
        db = _db(tmp_path, [
            _row("sib", "10", 87894059, "C", "G", "PTEN",
                 "NP_000305.3:p.Pro38Leu", "P38L", 38, "Pathogenic", 2),
        ])
        hits = query_same_aa_change(
            db, "PTEN", "NP_000305.3:p.Pro38Leu",
            exclude_chrom="chr10", exclude_pos=87894058,
            exclude_ref="C", exclude_alt="T",
        )
        assert [h.variation_id for h in hits] == ["sib"]

    def test_null_position_comparator_kept(self, tmp_path):
        # A comparator with no recorded position cannot be disproven as
        # same-codon → kept (preserves rows the build could not coordinate).
        db = _db(tmp_path, [
            _row("nopos", None, None, "C", "G", "GCK",
                 "NP_000153.1:p.Gly175Arg", "G175R", 175, "Pathogenic", 2),
        ])
        hits = query_same_aa_change(
            db, "GCK", "NP_000153.1:p.Gly175Arg",
            exclude_chrom="chr7", exclude_pos=44150025,
            exclude_ref="C", exclude_alt="T",
        )
        assert [h.variation_id for h in hits] == ["nopos"]
