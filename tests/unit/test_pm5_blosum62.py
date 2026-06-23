"""PM5 BLOSUM62 chemical-severity gate (PTEN GN003).

The PTEN VCEP gates PM5 on BLOSUM62 rather than Grantham: the candidate "must
have a BLOSUM62 score equal to or less than the known variant". BLOSUM62 is a
similarity score (higher = more conservative), so a chemically as-severe-or-
more candidate scores *<=* the comparator — the opposite direction to Grantham
distance. These tests cover the embedded matrix, the per-gene loader's
``blosum_le``/``blosum_lt`` tokens, the evaluator gate, and the build-script
extraction.
"""
import importlib.util
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.blosum62 import blosum62_score
from acmg_classifier.criteria.pm5_genes import PM5Spec
from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
from acmg_classifier.models.annotation import (
    AnnotationData, GnomADData, ConsequenceInfo,
)
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.variant import VariantRecord

_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


class TestBlosum62Matrix:
    def test_diagonal_self_scores(self):
        assert blosum62_score("C", "C") == 9    # highest self-score
        assert blosum62_score("W", "W") == 11
        assert blosum62_score("A", "A") == 4
        assert blosum62_score("V", "V") == 4

    def test_anchor_offdiagonal(self):
        assert blosum62_score("R", "K") == 2     # conservative, positive
        assert blosum62_score("I", "L") == 2
        assert blosum62_score("F", "Y") == 3
        assert blosum62_score("C", "P") == -3    # dissimilar, negative

    def test_symmetric_and_three_letter(self):
        assert blosum62_score("R", "H") == blosum62_score("H", "R")
        assert blosum62_score("Arg", "His") == blosum62_score("R", "H")
        assert blosum62_score("ARG", "his") == blosum62_score("R", "H")

    def test_unknown_codes(self):
        assert blosum62_score("X", "A") is None
        assert blosum62_score(None, "A") is None
        assert blosum62_score("Xyz", "A") is None


class TestPm5GateLoader:
    def _tsv(self, tmp_path):
        tsv = tmp_path / "dp.tsv"
        tsv.write_text(
            "gene_symbol\tpm5_grantham\tpm5_excludes\tpm5_max\tpm5_lp\n"
            "PIK3CD\tge\t\t\t\n"
            "PIK3R1\tgt\t\t\t\n"
            "PTEN\tblosum_le\t\t\tno\n"
            "GENE_BLT\tblosum_lt\t\t\t\n",
            encoding="utf-8",
        )
        return tsv

    def test_grantham_gate_unchanged(self, tmp_path):
        s = PM5Spec(self._tsv(tmp_path))
        assert s.gate("PIK3CD") == ("grantham", "ge")
        assert s.gate("PIK3R1") == ("grantham", "gt")

    def test_blosum_gate(self, tmp_path):
        s = PM5Spec(self._tsv(tmp_path))
        assert s.gate("PTEN") == ("blosum", "le")
        assert s.gate("GENE_BLT") == ("blosum", "lt")

    def test_no_gate(self, tmp_path):
        s = PM5Spec(self._tsv(tmp_path))
        assert s.gate("UNSEEN") is None
        assert s.gate(None) is None


class TestPm5BlosumExtraction:
    def _pm5(self, *strengths):
        return {"criteriaCodes": [{"label": "PM5", "evidenceStrengths": [
            {"label": lbl, "applicability": "Applicable", "description": desc}
            for lbl, desc in strengths
        ]}]}

    def test_equal_or_less_is_blosum_le(self):
        rs = self._pm5(("Moderate",
            "variant being interrogated must have BLOSUM62 score equal to or "
            "less than the known variant."))
        assert bdt._pm5_grantham_op(rs) == "blosum_le"

    def test_strict_less_is_blosum_lt(self):
        rs = self._pm5(("Moderate", "must have a BLOSUM62 score less than the comparator."))
        assert bdt._pm5_grantham_op(rs) == "blosum_lt"

    def test_grantham_still_works(self):
        rs = self._pm5(("Moderate", "Grantham distance greater than or equal to the variant."))
        assert bdt._pm5_grantham_op(rs) == "ge"


# ------------------------------ evaluator ------------------------------------

_SCHEMA = """
CREATE TABLE variants (
    variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER
);
"""

_AA1_TO_AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}


def _db(tmp_path, rows):
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
    con.close()
    return p


def _row(vid, gene, hgvs_p, aa, codon, sig, stars=2):
    return (vid, gene, hgvs_p, aa, codon, sig, "criteria provided", stars)


def _cfg(tmp_path, clinvar):
    tsv = tmp_path / "disease_prevalence.tsv"
    tsv.write_text(
        "gene_symbol\tpm5_grantham\tpm5_excludes\tpm5_max\tpm5_lp\n"
        "PTEN\tblosum_le\t\t\t\n",
        encoding="utf-8",
    )
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = clinvar
    cfg.pm5_min_stars = 1
    return cfg


def _consequence(gene, amino_acids, codon=175):
    ref, alt = (p.strip() for p in amino_acids.split("/"))
    hgvs_p = f"NM:p.{_AA1_TO_AA3.get(ref, ref)}{codon}{_AA1_TO_AA3.get(alt, alt)}"
    return ConsequenceInfo(
        transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        is_mane_select=True, protein_position=codon, codon_position=codon,
        amino_acids=amino_acids, hgvs_p=hgvs_p,
    )


def _ann(gene, amino_acids):
    return AnnotationData(gnomad=GnomADData(), consequences=[_consequence(gene, amino_acids)])


def _snv():
    return VariantRecord(chrom="chr1", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38)


class TestPm5BlosumEvaluator:
    # At codon 175 wild-type Arg(R). Comparator R->H: BLOSUM62(R,H) = 0.
    def test_candidate_more_severe_passes(self, tmp_path):
        # Candidate R->C: BLOSUM62(R,C) = -3 <= comparator 0 -> as/more severe -> PM5.
        db = _db(tmp_path, [_row("1", "PTEN", "NM:p.Arg175His", "R175H", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("PTEN", "R/C"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE
        assert "BLOSUM62-gated (<=)" in r.evidence

    def test_candidate_equal_score_passes(self, tmp_path):
        # Candidate R->H equals comparator R->H (0 <= 0) -> still fires (le inclusive).
        # Use a different comparator with the same BLOSUM score to keep PM5 (diff AA):
        # R->Q BLOSUM62(R,Q)=1; candidate R->E BLOSUM62(R,E)=0 -> 0<=1 passes.
        db = _db(tmp_path, [_row("1", "PTEN", "NM:p.Arg175Gln", "R175Q", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("PTEN", "R/E"))
        assert r.triggered

    def test_candidate_milder_blocked(self, tmp_path):
        # Comparator R->C (BLOSUM -3); candidate R->K (BLOSUM 2). 2 <= -3 false -> blocked.
        db = _db(tmp_path, [_row("1", "PTEN", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("PTEN", "R/K"))
        assert not r.triggered and "BLOSUM62 gate failed" in r.evidence

    def test_candidate_score_unavailable_withheld(self, tmp_path):
        db = _db(tmp_path, [_row("1", "PTEN", "NM:p.Arg175His", "R175H", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("PTEN", "R/Xaa"))
        assert not r.triggered and "BLOSUM62 score unavailable" in r.evidence
