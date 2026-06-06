"""PS1 splicing extension: same-splice-site (position-based) matching.

PS1's amino-acid rule cannot fire for intronic/splice variants (no protein
change). The ClinGen SVI splicing extension recognises a DIFFERENT nucleotide
change at the SAME splice-site position as having the same predicted effect.
"""
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator, _ps1_strength
from acmg_classifier.local_db.clinvar_sqlite import query_same_splice_site
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo, ClinVarRecord
from acmg_classifier.models.enums import Assembly, ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord


def _rec(sig: str) -> ClinVarRecord:
    return ClinVarRecord(
        variation_id="1", clinical_significance=sig, review_status="x", star_rating=2,
    )


class TestPS1Strength:
    """ClinGen SVI: Pathogenic comparator -> Strong; Likely-pathogenic only -> Moderate."""

    def test_pathogenic_comparator_strong(self):
        assert _ps1_strength([_rec("Pathogenic")]) == CriterionStrength.STRONG

    def test_lp_only_comparator_moderate(self):
        assert _ps1_strength([_rec("Likely pathogenic")]) == CriterionStrength.MODERATE

    def test_p_lp_aggregate_is_strong(self):
        # 'Pathogenic/Likely pathogenic' contains a P assertion -> Strong.
        assert _ps1_strength([_rec("Pathogenic/Likely pathogenic")]) == CriterionStrength.STRONG

    def test_any_pathogenic_among_lp_wins_strong(self):
        hits = [_rec("Likely pathogenic"), _rec("Pathogenic")]
        assert _ps1_strength(hits) == CriterionStrength.STRONG

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


# Standard PS1 spec fixture: HNF1A extends to canonical splice, MSH2 only to
# non-canonical, GAA is missense-only, CDH1 declines PS1 entirely.
_PS1_TSV = (
    "gene_symbol\tps1\tps1_splice\n"
    "HNF1A\tapplicable\tcanonical\n"
    "MSH2\tapplicable\tnoncanonical\n"
    "GAA\tapplicable\t\n"
    "CDH1\tnot_applicable\t\n"
)


def _cfg(db, tmp_path):
    cfg = MagicMock()
    cfg.clinvar_sqlite = db
    tsv = tmp_path / "disease_prevalence.tsv"
    tsv.write_text(_PS1_TSV, encoding="utf-8")
    cfg.disease_prevalence_tsv = tsv
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
        r = PS1Evaluator(_cfg(db, tmp_path)).evaluate(_variant(), _ann(ConsequenceType.SPLICE_DONOR))
        assert r.triggered
        assert "same splice-site position" in r.evidence

    def test_splice_region_noncanonical_supported(self, tmp_path):
        # c.xxx+5 (splice_region) recovers PS1 when a P/LP variant shares the position.
        db = _db(tmp_path, [
            _row("1", "12", 999, "A", "G", "HNF1A", "c.1+5A>G", "Likely pathogenic"),
        ])
        r = PS1Evaluator(_cfg(db, tmp_path)).evaluate(
            _variant(pos=999, ref="A", alt="T"), _ann(ConsequenceType.SPLICE_REGION)
        )
        assert r.triggered

    def test_no_same_site_hit_not_met(self, tmp_path):
        db = _db(tmp_path, [])
        r = PS1Evaluator(_cfg(db, tmp_path)).evaluate(_variant(), _ann(ConsequenceType.SPLICE_DONOR))
        assert not r.triggered
        assert "same-splice-site" in r.evidence

    def test_deep_intron_uses_position_path(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 500, "A", "C", "HNF1A", "c.1-12A>C", "Pathogenic"),
        ])
        r = PS1Evaluator(_cfg(db, tmp_path)).evaluate(
            _variant(pos=500, ref="A", alt="G"), _ann(ConsequenceType.INTRON)
        )
        assert r.triggered


class TestPs1NotApplicable:
    """A VCEP may decline PS1 entirely for its gene (CDH1)."""

    def test_not_applicable_blocks_missense(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "16", 68800000, "C", "G", "CDH1", "c.1A>G", "Pathogenic"),
        ])
        tsv = tmp_path / "disease_prevalence.tsv"
        tsv.write_text("gene_symbol\tps1\nCDH1\tnot_applicable\n", encoding="utf-8")
        cfg = MagicMock()
        cfg.clinvar_sqlite = db
        cfg.disease_prevalence_tsv = tsv
        r = PS1Evaluator(cfg).evaluate(_variant(), _ann(ConsequenceType.MISSENSE, gene="CDH1"))
        assert not r.triggered
        assert "not applicable" in r.evidence


class TestPs1SpliceNonCanonicalOnly:
    """InSiGHT MMR genes restrict PS1-splice to non-canonical nucleotides; a
    canonical ±1/±2 variant must not receive PS1 (it is PVS1 territory)."""

    def test_canonical_blocked(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "C", "MSH2", "c.1+1G>C", "Pathogenic"),
        ])
        cfg = _cfg(db, tmp_path)
        r = PS1Evaluator(cfg).evaluate(_variant(), _ann(ConsequenceType.SPLICE_DONOR, gene="MSH2"))
        assert not r.triggered
        assert "non-canonical" in r.evidence

    def test_noncanonical_still_fires(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 999, "A", "G", "MSH2", "c.1+5A>G", "Pathogenic"),
        ])
        cfg = _cfg(db, tmp_path)
        r = PS1Evaluator(cfg).evaluate(
            _variant(pos=999, ref="A", alt="T"), _ann(ConsequenceType.SPLICE_REGION, gene="MSH2")
        )
        assert r.triggered

    def test_unrestricted_gene_canonical_fires(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "C", "HNF1A", "c.1+1G>C", "Pathogenic"),
        ])
        cfg = _cfg(db, tmp_path)  # restriction not for HNF1A
        r = PS1Evaluator(cfg).evaluate(_variant(), _ann(ConsequenceType.SPLICE_DONOR, gene="HNF1A"))
        assert r.triggered


class TestPs1SpliceMissenseOnly:
    """Genes whose PS1 is the original missense-only ACMG rule (no splice
    extension — e.g. GAA, the HCM genes) must not receive PS1 for a splice
    variant even when a same-site P/LP comparator exists."""

    def test_missense_only_gene_splice_withheld(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "17", 80110841, "G", "C", "GAA", "c.1551+1G>C", "Pathogenic"),
        ])
        cfg = _cfg(db, tmp_path)  # GAA -> ps1_splice blank (missense-only)
        r = PS1Evaluator(cfg).evaluate(
            VariantRecord(chrom="chr17", pos=80110841, ref="G", alt="T", assembly=Assembly.GRCH38),
            _ann(ConsequenceType.SPLICE_DONOR, gene="GAA"),
        )
        assert not r.triggered
        assert "missense-only" in r.evidence

    def test_gene_absent_from_tsv_is_missense_only(self, tmp_path):
        # A gene with no VCEP PS1 spec defaults to missense-only for splice.
        db = _db(tmp_path, [
            _row("1", "12", 120989033, "G", "C", "NOVCEP", "c.1+1G>C", "Pathogenic"),
        ])
        cfg = _cfg(db, tmp_path)
        r = PS1Evaluator(cfg).evaluate(_variant(), _ann(ConsequenceType.SPLICE_DONOR, gene="NOVCEP"))
        assert not r.triggered and "missense-only" in r.evidence
