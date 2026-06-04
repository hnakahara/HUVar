"""PS1 splicing extension: same-splice-site (position-based) matching.

PS1's amino-acid rule cannot fire for intronic/splice variants (no protein
change). The ClinGen SVI splicing extension recognises a DIFFERENT nucleotide
change at the SAME splice-site position as having the same predicted effect.
"""
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator
from acmg_classifier.local_db.clinvar_sqlite import query_same_splice_site
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ACMGCriterion, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_SCHEMA = """
CREATE TABLE variants (
    variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER
);
"""


def _db(tmp_path, rows):
    p = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, "
        "hgvs_c, hgvs_p, clinical_significance, review_status, star_rating) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        rows,
    )
    con.commit()
    con.close()
    return p


def _row(vid, chrom, pos, ref, alt, gene, hgvs_c, sig, stars=1):
    return (vid, chrom, pos, ref, alt, gene, hgvs_c, None, sig, "criteria provided", stars)


class TestQuerySameSpliceSite:
    def test_same_pos_diff_alt_pathogenic_hit(self, tmp_path):
        # HNF1A c.526+1: G>C is P/LP; query for the G>T variant finds it.
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "C", "HNF1A", "c.526+1G>C", "Pathogenic"),
        ])
        hits = query_same_splice_site(db, "HNF1A", "chr12", 120989033, "G", "T")
        assert [h.variation_id for h in hits] == ["1"]

    def test_self_excluded(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "T", "HNF1A", "c.526+1G>T", "Pathogenic"),
        ])
        hits = query_same_splice_site(db, "HNF1A", "chr12", 120989033, "G", "T")
        assert hits == []

    def test_benign_not_counted(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "C", "HNF1A", "c.526+1G>C", "Benign"),
        ])
        assert query_same_splice_site(db, "HNF1A", "chr12", 120989033, "G", "T") == []

    def test_star_threshold(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "C", "HNF1A", "c.526+1G>C", "Pathogenic", stars=0),
        ])
        assert query_same_splice_site(db, "HNF1A", "chr12", 120989033, "G", "T", min_stars=1) == []

    def test_other_gene_excluded(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "C", "OTHER", "c.1G>C", "Pathogenic"),
        ])
        assert query_same_splice_site(db, "HNF1A", "chr12", 120989033, "G", "T") == []


def _cfg(db):
    cfg = MagicMock()
    cfg.clinvar_sqlite = db
    return cfg


def _ann(consequence, gene="HNF1A"):
    return AnnotationData(
        gnomad=GnomADData(),
        consequences=[ConsequenceInfo(
            transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
            consequence=consequence, biotype="protein_coding", is_mane_select=True,
        )],
    )


def _variant(pos=120989033, ref="G", alt="T"):
    return VariantRecord(chrom="chr12", pos=pos, ref=ref, alt=alt, assembly=Assembly.GRCH38)


class TestPs1SpliceEvaluator:
    def test_splice_donor_fires_on_same_site(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "C", "HNF1A", "c.526+1G>C", "Pathogenic"),
        ])
        r = PS1Evaluator(_cfg(db)).evaluate(_variant(), _ann(ConsequenceType.SPLICE_DONOR))
        assert r.triggered
        assert "same splice-site position" in r.evidence

    def test_splice_region_noncanonical_supported(self, tmp_path):
        # c.xxx+5 (splice_region) recovers PS1 when a P/LP variant shares the position.
        db = _db(tmp_path, [
            _row("1", "12", 999, "A", "G", "HNF1A", "c.1+5A>G", "Likely pathogenic"),
        ])
        r = PS1Evaluator(_cfg(db)).evaluate(
            _variant(pos=999, ref="A", alt="T"), _ann(ConsequenceType.SPLICE_REGION)
        )
        assert r.triggered

    def test_no_same_site_hit_not_met(self, tmp_path):
        db = _db(tmp_path, [])
        r = PS1Evaluator(_cfg(db)).evaluate(_variant(), _ann(ConsequenceType.SPLICE_DONOR))
        assert not r.triggered
        assert "same-splice-site" in r.evidence

    def test_deep_intron_uses_position_path(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 500, "A", "C", "HNF1A", "c.1-12A>C", "Pathogenic"),
        ])
        r = PS1Evaluator(_cfg(db)).evaluate(
            _variant(pos=500, ref="A", alt="G"), _ann(ConsequenceType.INTRON)
        )
        assert r.triggered
