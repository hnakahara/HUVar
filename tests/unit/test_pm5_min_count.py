"""PM5 minimum-comparator-count rule (ACVRL1/ENG, HHT VCEP GN135/136).

The HHT VCEP applies PM5_Strong when >=2 DIFFERENT same-codon missense changes
have been determined likely pathogenic or pathogenic, and PM5_Moderate for one.
LP comparators count like P (pm5_lp cleared); there is no Supporting level. The
splicing caveat is a note in the evidence, not an exclusion.
"""
import sqlite3
from unittest.mock import MagicMock

from acmg_classifier.criteria.pm5_genes import PM5Spec
from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_SCHEMA = """
CREATE TABLE variants (
    variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER
);
"""
_AA1_TO_AA3 = {"R": "Arg", "H": "His", "Q": "Gln", "C": "Cys", "K": "Lys"}


def _db(tmp_path, rows):
    p = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, gene_symbol, hgvs_p, amino_acid_change, "
        "codon_position, clinical_significance, review_status, star_rating) "
        "VALUES (?,?,?,?,?,?,?,?)",
        rows,
    )
    con.commit()
    con.close()
    return p


def _row(vid, gene, hgvs_p, aa, sig, codon=175, stars=2):
    return (vid, gene, hgvs_p, aa, codon, sig, "criteria provided", stars)


def _cfg(tmp_path, clinvar):
    tsv = tmp_path / "dp.tsv"
    tsv.write_text(
        "gene_symbol\tpm5_grantham\tpm5_max\tpm5_lp\tpm5_min_count\n"
        "ACVRL1\t\tStrong\t\t2\n",
        encoding="utf-8",
    )
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = clinvar
    cfg.pm5_min_stars = 1
    return cfg


def _ann(gene, amino_acids, codon=175):
    ref, alt = (p.strip() for p in amino_acids.split("/"))
    hgvs_p = f"NM:p.{_AA1_TO_AA3.get(ref, ref)}{codon}{_AA1_TO_AA3.get(alt, alt)}"
    return AnnotationData(gnomad=GnomADData(), consequences=[ConsequenceInfo(
        transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        is_mane_select=True, protein_position=codon, codon_position=codon,
        amino_acids=amino_acids, hgvs_p=hgvs_p,
    )])


def _snv():
    return VariantRecord(chrom="chr1", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38)


class TestLoader:
    def test_min_count_parsed(self, tmp_path):
        s = PM5Spec(_cfg(tmp_path, None).disease_prevalence_tsv)
        assert s.min_count("ACVRL1") == 2
        assert s.min_count("OTHER") == 1     # default
        assert not s.requires_pathogenic("ACVRL1")   # pm5_lp cleared -> LP allowed


class TestEvaluator:
    def test_two_distinct_lp_p_strong(self, tmp_path):
        # Candidate R->C; two distinct same-codon comparators (one P, one LP).
        db = _db(tmp_path, [
            _row("1", "ACVRL1", "NM:p.Arg175His", "R175H", "Pathogenic"),
            _row("2", "ACVRL1", "NM:p.Arg175Gln", "R175Q", "Likely pathogenic"),
        ])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("ACVRL1", "R/C"))
        assert r.triggered and r.strength == CriterionStrength.STRONG

    def test_two_distinct_lp_only_strong(self, tmp_path):
        # Both comparators only Likely pathogenic -> still >=2 distinct -> Strong.
        db = _db(tmp_path, [
            _row("1", "ACVRL1", "NM:p.Arg175His", "R175H", "Likely pathogenic"),
            _row("2", "ACVRL1", "NM:p.Arg175Gln", "R175Q", "Likely pathogenic"),
        ])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("ACVRL1", "R/C"))
        assert r.triggered and r.strength == CriterionStrength.STRONG

    def test_single_comparator_moderate(self, tmp_path):
        # One distinct same-codon comparator -> PM5_Moderate (HHT VCEP).
        db = _db(tmp_path, [
            _row("1", "ACVRL1", "NM:p.Arg175His", "R175H", "Pathogenic"),
        ])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("ACVRL1", "R/C"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE
        assert "1 distinct" in r.evidence

    def test_single_lp_comparator_moderate_not_supporting(self, tmp_path):
        # LP counts like P: one LP comparator -> Moderate, never Supporting.
        db = _db(tmp_path, [
            _row("1", "ACVRL1", "NM:p.Arg175His", "R175H", "Likely pathogenic"),
            # same change submitted twice -> still one distinct change
            _row("3", "ACVRL1", "NM:p.Arg175His", "R175H", "Pathogenic"),
        ])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("ACVRL1", "R/C"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_no_comparator_not_met(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "ACVRL1", "NM:p.Arg175His", "R175H", "Uncertain significance"),
        ])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("ACVRL1", "R/C"))
        assert not r.triggered

    def test_splice_caveat_is_a_note_only(self, tmp_path):
        from acmg_classifier.models.annotation import SpliceScore
        db = _db(tmp_path, [
            _row("1", "ACVRL1", "NM:p.Arg175His", "R175H", "Pathogenic"),
            _row("2", "ACVRL1", "NM:p.Arg175Gln", "R175Q", "Pathogenic"),
        ])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        ann = _ann("ACVRL1", "R/C")
        ann.splice = _splice(SpliceScore, 0.45)
        r = ev.evaluate(_snv(), ann)
        assert r.triggered and r.strength == CriterionStrength.STRONG
        assert "VCEP caveat" in r.evidence and "max_delta=0.45" in r.evidence
        ann.splice = _splice(SpliceScore, 0.05)
        r = ev.evaluate(_snv(), ann)
        assert r.triggered and "caveat" not in r.evidence
        ann.splice = None
        r = ev.evaluate(_snv(), ann)
        assert r.triggered and "splice prediction unavailable" in r.evidence


def _splice(cls, delta):
    return cls(tool="openspliceai", max_delta=delta, is_available=True)
