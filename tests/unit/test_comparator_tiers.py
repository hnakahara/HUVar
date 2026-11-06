"""PS1 / PM5 strengths that depend on how many P / LP comparators there are
(MYOC, ACADVL, VWF, F8/F9, LGMD genes, RASopathy genes): criteria/comparator_tiers.py."""
import sqlite3
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.comparator_tiers import (
    PM5_RULES, PS1_RULES, count_distinct, tier_strength,
)
from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength as S
from acmg_classifier.models.variant import VariantRecord

N = None


@pytest.mark.parametrize("gene,n_p,n_lp,want", [
    # MYOC: Strong 2P / 1P+2LP / 4LP; Moderate 1P / 2LP; Supporting 1LP
    ("MYOC", 2, 0, S.STRONG), ("MYOC", 1, 2, S.STRONG), ("MYOC", 0, 4, S.STRONG),
    ("MYOC", 1, 0, S.MODERATE), ("MYOC", 0, 2, S.MODERATE), ("MYOC", 0, 3, S.MODERATE),
    ("MYOC", 0, 1, S.SUPPORTING), ("MYOC", 0, 0, N),
    # ACADVL table: LP 0.5, P 1.0; Strong >= 2, Moderate >= 1, Supporting 0.5
    ("ACADVL", 2, 0, S.STRONG), ("ACADVL", 1, 2, S.STRONG), ("ACADVL", 1, 1, S.MODERATE),
    ("ACADVL", 0, 2, S.MODERATE), ("ACADVL", 0, 1, S.SUPPORTING),
    # LGMD: Strong 2P / 3LP; Moderate 1P / 2LP; Supporting 1LP
    ("DYSF", 2, 0, S.STRONG), ("CAPN3", 0, 3, S.STRONG), ("SGCA", 1, 0, S.MODERATE),
    ("ANO5", 0, 2, S.MODERATE), ("SGCB", 0, 1, S.SUPPORTING),
    # VWF / F8 / F9: Moderate 1P / 2LP; Supporting 1LP; no Strong
    ("VWF", 3, 0, S.MODERATE), ("VWF", 0, 2, S.MODERATE), ("VWF", 0, 1, S.SUPPORTING),
    ("F8", 0, 2, S.MODERATE), ("F9", 0, 1, S.SUPPORTING),
    # RASopathy: [L]P alike; Moderate for one; Strong capped at Moderate
    ("KRAS", 0, 1, S.MODERATE), ("PTPN11", 1, 0, S.MODERATE), ("LZTR1", 2, 3, S.MODERATE),
])
def test_pm5_tiers(gene, n_p, n_lp, want):
    assert tier_strength(PM5_RULES[gene], n_p, n_lp)[0] == want


def test_rasopathy_cap_note():
    st, text = tier_strength(PM5_RULES["KRAS"], 1, 1)
    assert st == S.MODERATE and ">=5 probands" in text and "capped at Moderate" in text
    st, text = tier_strength(PM5_RULES["KRAS"], 1, 0)
    assert st == S.MODERATE and "probands" not in text


@pytest.mark.parametrize("gene", ["DYSF", "SGCG", "F8", "F9", "ACADVL"])
def test_ps1_tiers(gene):
    rule = PS1_RULES[gene]
    assert tier_strength(rule, 1, 0)[0] == S.STRONG
    assert tier_strength(rule, 0, 2)[0] == S.STRONG
    assert tier_strength(rule, 0, 1)[0] == S.MODERATE


def test_count_distinct_best_classification():
    hits = [SimpleNamespace(clinical_significance=s, k=k) for s, k in [
        ("Likely pathogenic", "a"), ("Pathogenic", "a"), ("Likely pathogenic", "b"),
        ("Pathogenic/Likely pathogenic", "c"), ("Uncertain significance", "d")]]
    assert count_distinct(hits, lambda h: h.k) == (2, 1)


# --- evaluators -------------------------------------------------------------

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
        "INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, hgvs_p, "
        "amino_acid_change, codon_position, clinical_significance, review_status, "
        "star_rating) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit()
    con.close()
    return p


def _cfg(tmp_path, clinvar, tsv_rows=""):
    tsv = tmp_path / "dp.tsv"
    tsv.write_text("gene_symbol\tpm5_grantham\tpm5_max\tpm5_lp\tpm5_min_count\tps1\n"
                   + tsv_rows, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = clinvar
    cfg.pm5_min_stars = 1
    cfg.ps1_min_stars = 1
    return cfg


def _ann(gene, aa="G/R", hgvs_p="NP_0:p.Gly175Arg"):
    return AnnotationData(gnomad=GnomADData(), consequences=[ConsequenceInfo(
        transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        is_mane_select=True, protein_position=175, codon_position=175,
        amino_acids=aa, hgvs_p=hgvs_p)])


def _snv():
    return VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)


def _pm5_row(vid, gene, change, sig):
    return (vid, None, None, None, None, gene, f"NP_0:p.Gly175{change}", f"G175{change}",
            175, sig, "criteria provided", 2)


def test_pm5_myoc_two_lp_moderate(tmp_path):
    db = _db(tmp_path, [_pm5_row("1", "MYOC", "Glu", "Likely pathogenic"),
                        _pm5_row("2", "MYOC", "Val", "Likely pathogenic")])
    r = PM5Evaluator(_cfg(tmp_path, db, "MYOC\t\tStrong\t\t\t\n")).evaluate(
        _snv(), _ann("MYOC"))
    assert r.triggered and r.strength == S.MODERATE and "0 P + 2 LP" in r.evidence


def test_pm5_rasopathy_lp_counts_and_strong_capped(tmp_path):
    # pm5_lp=no in the TSV no longer drops LP comparators for the RASopathy genes.
    db = _db(tmp_path, [_pm5_row("1", "KRAS", "Glu", "Likely pathogenic")])
    cfg = _cfg(tmp_path, db, "KRAS\t\tStrong\tno\t\t\n")
    r = PM5Evaluator(cfg).evaluate(_snv(), _ann("KRAS"))
    assert r.triggered and r.strength == S.MODERATE
    (tmp_path / "x").mkdir()
    db = _db(tmp_path / "x", [_pm5_row("1", "KRAS", "Glu", "Pathogenic"),
                              _pm5_row("2", "KRAS", "Val", "Pathogenic")])
    cfg = _cfg(tmp_path, db, "KRAS\t\tStrong\tno\t\t\n")
    r = PM5Evaluator(cfg).evaluate(_snv(), _ann("KRAS"))
    assert r.triggered and r.strength == S.MODERATE and "probands" in r.evidence


def test_pm5_vwf_tsv_cap_kept(tmp_path):
    db = _db(tmp_path, [_pm5_row("1", "VWF", "Glu", "Pathogenic"),
                        _pm5_row("2", "VWF", "Val", "Pathogenic")])
    r = PM5Evaluator(_cfg(tmp_path, db, "VWF\t\tModerate\t\t\t\n")).evaluate(
        _snv(), _ann("VWF"))
    assert r.triggered and r.strength == S.MODERATE


def _ps1_row(vid, gene, pos, alt, sig):
    return (vid, "chr1", pos, "C", alt, gene, "NP_0:p.Gly175Arg", "G175R", 175, sig,
            "criteria provided", 2)


@pytest.fixture()
def no_vcep_filter(monkeypatch):
    # The LGMD VCEPs require VCEP-classified comparators (vcep_comparator_rules.tsv);
    # these tests exercise the counting, so every ClinVar comparator is accepted.
    from acmg_classifier.criteria import vcep_comparators
    monkeypatch.setattr(vcep_comparators, "filter_comparators",
                        lambda gene, crit, hits, asm: (hits, ""))


def test_ps1_dysf_two_lp_strong(tmp_path, no_vcep_filter):
    db = _db(tmp_path, [_ps1_row("1", "DYSF", 100, "G", "Likely pathogenic"),
                        _ps1_row("2", "DYSF", 102, "A", "Likely pathogenic")])
    r = PS1Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("DYSF"))
    assert r.triggered and r.strength == S.STRONG and "0 P + 2 LP" in r.evidence


def test_ps1_f8_one_lp_moderate_and_default_gene_unchanged(tmp_path, no_vcep_filter):
    db = _db(tmp_path, [_ps1_row("1", "F8", 100, "G", "Likely pathogenic"),
                        _ps1_row("2", "OTHER", 100, "G", "Likely pathogenic"),
                        _ps1_row("3", "OTHER", 102, "A", "Likely pathogenic")])
    cfg = _cfg(tmp_path, db)
    r = PS1Evaluator(cfg).evaluate(_snv(), _ann("F8"))
    assert r.triggered and r.strength == S.MODERATE
    # Genes without a count rule keep the SVI tiering: LP only -> Moderate.
    r = PS1Evaluator(cfg).evaluate(_snv(), _ann("OTHER"))
    assert r.triggered and r.strength == S.MODERATE
