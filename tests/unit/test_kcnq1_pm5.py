"""KCNQ1 PM5 (LQTS VCEP): residue conserved in KCNQ1-5, P or LP comparator ->
Moderate, both variants SpliceAI < 0.2."""
import sqlite3
from unittest.mock import MagicMock

from acmg_classifier.criteria import pm5_conditions as cond
from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, GnomADData, SpliceScore,
)
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength as S
from acmg_classifier.models.variant import VariantRecord

_SCHEMA = """CREATE TABLE variants (variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT,
    alt TEXT, gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER)"""
_3 = {"R": "Arg", "W": "Trp", "Q": "Gln", "H": "His", "P": "Pro", "L": "Leu", "E": "Glu"}


def _cfg(tmp_path, rows):
    db = tmp_path / "cv.sqlite"
    con = sqlite3.connect(db)
    con.execute(_SCHEMA)
    con.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit()
    con.close()
    tsv = tmp_path / "dp.tsv"
    tsv.write_text("gene_symbol\tpm5_grantham\tpm5_max\tpm5_lp\tpm5_min_count\n"
                   "KCNQ1\t\tModerate\tno\t\n", encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = db
    cfg.pm5_min_stars = 1
    return cfg


def _row(vid, pos, ref, alt, sig):
    return (vid, None, None, None, None, "KCNQ1", None, f"NP_000209.2:p.{_3[ref]}{pos}{_3[alt]}",
            f"{ref}{pos}{alt}", pos, sig, "criteria provided", 2)


def _ann(pos, ref, alt, splice=0.01):
    return AnnotationData(
        gnomad=GnomADData(), splice=SpliceScore(tool="openspliceai", max_delta=splice,
                                                is_available=True),
        consequences=[ConsequenceInfo(
            transcript_id="NM_000218.3", gene_id="ENSG", gene_symbol="KCNQ1",
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
            is_mane_select=True, protein_position=pos, codon_position=pos,
            amino_acids=f"{ref}/{alt}", hgvs_p=f"NP_000209.2:p.{_3[ref]}{pos}{_3[alt]}")])


_V = VariantRecord(chrom="chr11", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38)


def test_table():
    assert cond.paralog_conserved("KCNQ1", 243, "R")[0] is True
    assert cond.paralog_conserved("KCNQ1", 357, "Q")[0] is False
    assert cond.paralog_conserved("KCNQ1", 243, "Q")[0] is None   # residue mismatch
    assert cond.paralog_conserved("BRCA1", 243, "R") == (None, "")


def test_conserved_lp_comparator_moderate(tmp_path):
    cfg = _cfg(tmp_path, [_row("1", 243, "R", "H", "Likely pathogenic")])
    r = PM5Evaluator(cfg).evaluate(_V, _ann(243, "R", "W"))
    assert r.triggered and r.strength == S.MODERATE and "conserved in KCNQ1-5" in r.evidence


def test_not_conserved_not_met(tmp_path):
    cfg = _cfg(tmp_path, [_row("1", 357, "Q", "P", "Pathogenic")])
    r = PM5Evaluator(cfg).evaluate(_V, _ann(357, "Q", "L"))
    assert not r.triggered and "not conserved" in r.evidence


def test_candidate_splice(tmp_path):
    cfg = _cfg(tmp_path, [_row("1", 243, "R", "H", "Pathogenic")])
    r = PM5Evaluator(cfg).evaluate(_V, _ann(243, "R", "W", splice=0.35))
    assert not r.triggered and "splicing" in r.evidence
