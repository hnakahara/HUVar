"""PM5 Grantham-distance gating: matrix, cspec extraction, loader, evaluator.

A subset of ClinGen VCEPs (PIK3CD, PIK3R1, RYR1, VHL, HNF1A, …) require PM5 to
clear a Grantham-distance test against the same-codon comparator. These tests
cover the embedded Grantham 1974 matrix, the ``pm5_grantham`` cspec extraction,
the per-gene loader, the benign-at-codon caveat, and the evaluator's gate and
strength assignment.
"""
import importlib.util
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.grantham import grantham_distance
from acmg_classifier.criteria.pm5_genes import PM5Grantham
from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
from acmg_classifier.criteria.registry import CriteriaRegistry, _PM5_EXCLUSIONS
from acmg_classifier.local_db.clinvar_sqlite import has_benign_at_codon
from acmg_classifier.models.annotation import (
    AnnotationData, GnomADData, ConsequenceInfo,
)
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    Assembly, ACMGCriterion, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.variant import VariantRecord


_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


class TestGranthamMatrix:
    def test_anchor_values(self):
        # Canonical Grantham 1974 Table 2 anchors.
        assert grantham_distance("I", "L") == 5      # most similar
        assert grantham_distance("C", "W") == 215    # most different
        assert grantham_distance("R", "H") == 29
        assert grantham_distance("S", "R") == 110

    def test_symmetric_and_three_letter(self):
        assert grantham_distance("R", "H") == grantham_distance("H", "R")
        assert grantham_distance("Arg", "His") == 29
        assert grantham_distance("ARG", "his") == 29  # case-insensitive

    def test_identity_and_unknown(self):
        assert grantham_distance("A", "A") == 0
        assert grantham_distance("X", "A") is None
        assert grantham_distance(None, "A") is None
        assert grantham_distance("Xyz", "A") is None


class TestPm5GranthamExtraction:
    def _pm5_code(self, *strengths):
        return {"label": "PM5", "evidenceStrengths": [
            {"label": lbl, "applicability": "Applicable", "description": desc}
            for lbl, desc in strengths
        ]}

    def test_inclusive_ge(self):
        rs = {"criteriaCodes": [self._pm5_code(
            ("Moderate", "Variant must have a Grantham distance greater than or "
                         "equal to the previously classified pathogenic variant."),
        )]}
        assert bdt._pm5_grantham_op(rs) == "ge"

    def test_strict_gt_higher(self):
        rs = {"criteriaCodes": [self._pm5_code(
            ("Moderate", "The variant of interest must have a higher Grantham "
                         "score than the Pathogenic comparison variant."),
        )]}
        assert bdt._pm5_grantham_op(rs) == "gt"

    def test_strict_gt_less_than(self):
        rs = {"criteriaCodes": [self._pm5_code(
            ("Moderate", "Grantham score for alternate pathogenic variant must "
                         "be less than for variant being assessed."),
        )]}
        assert bdt._pm5_grantham_op(rs) == "gt"

    def test_ge_wins_when_mixed(self):
        # A gene whose strengths mix an "equal" clause and a strict clause keeps
        # the inclusive operator (the spec's primary rule).
        rs = {"criteriaCodes": [self._pm5_code(
            ("Moderate", "Grantham distance greater than or equal to the variant."),
            ("Supporting", "pathogenic but has a greater Grantham distance."),
        )]}
        assert bdt._pm5_grantham_op(rs) == "ge"

    def test_no_grantham_is_empty(self):
        rs = {"criteriaCodes": [self._pm5_code(
            ("Moderate", "Novel missense at a residue with a known pathogenic change."),
        )]}
        assert bdt._pm5_grantham_op(rs) == ""

    def test_no_pm5_code_is_empty(self):
        assert bdt._pm5_grantham_op({"criteriaCodes": [{"label": "PM1"}]}) == ""


class TestPm5GranthamLoader:
    def _tsv(self, tmp_path):
        tsv = tmp_path / "dp.tsv"
        tsv.write_text(
            "gene_symbol\tpm5_grantham\n"
            "PIK3CD\tge\n"
            "PIK3R1\tgt\n"
            "BRCA1\t\n"
            "BADVAL\tnonsense\n",
            encoding="utf-8",
        )
        return tsv

    def test_operators_loaded(self, tmp_path):
        g = PM5Grantham(self._tsv(tmp_path))
        assert g.operator("PIK3CD") == "ge"
        assert g.operator("PIK3R1") == "gt"
        assert g.operator("BRCA1") == ""       # blank cell
        assert g.operator("BADVAL") == ""      # invalid value ignored
        assert g.operator("UNSEEN") == ""
        assert g.operator(None) == ""

    def test_missing_file_is_empty(self, tmp_path):
        g = PM5Grantham(tmp_path / "nope.tsv")
        assert g.operator("PIK3CD") == ""


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
        "INSERT INTO variants (variation_id, gene_symbol, hgvs_p, "
        "amino_acid_change, codon_position, clinical_significance, "
        "review_status, star_rating) VALUES (?,?,?,?,?,?,?,?)",
        rows,
    )
    con.commit()
    con.close()
    return p


def _row(vid, gene, hgvs_p, aa, codon, sig, stars=1):
    return (vid, gene, hgvs_p, aa, codon, sig, "criteria provided", stars)


class TestHasBenignAtCodon:
    def test_benign_missense_detected(self, tmp_path):
        db = _db(tmp_path, [_row("1", "G", "NM:p.Arg175Ser", "R175S", 175, "Benign")])
        assert has_benign_at_codon(db, "G", 175) is True

    def test_truncating_benign_ignored(self, tmp_path):
        db = _db(tmp_path, [_row("1", "G", "NM:p.Arg175Ter", "R175*", 175, "Likely benign")])
        assert has_benign_at_codon(db, "G", 175) is False

    def test_no_benign_returns_false(self, tmp_path):
        db = _db(tmp_path, [_row("1", "G", "NM:p.Arg175His", "R175H", 175, "Pathogenic")])
        assert has_benign_at_codon(db, "G", 175) is False


def _cfg(tmp_path, clinvar: Path):
    tsv = tmp_path / "disease_prevalence.tsv"
    tsv.write_text(
        "gene_symbol\tpm5_grantham\n"
        "PIK3CD\tge\n"
        "RYR1\tgt\n"
        "BRCA1\t\n",
        encoding="utf-8",
    )
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = clinvar
    return cfg


_AA1_TO_AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}


def _consequence(gene, amino_acids, codon=175, hgvs_p=None):
    # Derive the candidate's own hgvs_p from its REF/ALT pair so it never
    # collides with the comparator (a same-AA hit would be excluded as PS1).
    if hgvs_p is None and amino_acids and "/" in amino_acids:
        ref, alt = (p.strip() for p in amino_acids.split("/"))
        hgvs_p = f"NM:p.{_AA1_TO_AA3[ref]}{codon}{_AA1_TO_AA3[alt]}"
    return ConsequenceInfo(
        transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        is_mane_select=True, protein_position=codon, codon_position=codon,
        amino_acids=amino_acids, hgvs_p=hgvs_p,
    )


def _ann(gene, amino_acids, **kw):
    return AnnotationData(
        gnomad=GnomADData(), consequences=[_consequence(gene, amino_acids, **kw)]
    )


def _snv():
    return VariantRecord(chrom="chr1", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38)


class TestPm5GranthamEvaluator:
    def test_ge_pass_pathogenic_is_moderate(self, tmp_path):
        # candidate R->C (180) >= comparator R->H (29), comparator Pathogenic.
        db = _db(tmp_path, [_row("1", "PIK3CD", "NM:p.Arg175His", "R175H", 175, "Pathogenic")])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_snv(), _ann("PIK3CD", "R/C"))
        assert r.triggered
        assert r.strength == CriterionStrength.MODERATE
        assert "Grantham-gated" in r.evidence

    def test_ge_fail_when_candidate_milder(self, tmp_path):
        # candidate R->H (29) NOT >= comparator R->C (180).
        db = _db(tmp_path, [_row("1", "PIK3CD", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic")])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_snv(), _ann("PIK3CD", "R/H"))
        assert not r.triggered
        assert "Grantham gate failed" in r.evidence

    def test_likely_pathogenic_comparator_is_supporting(self, tmp_path):
        db = _db(tmp_path, [_row("1", "PIK3CD", "NM:p.Arg175His", "R175H", 175, "Likely pathogenic")])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_snv(), _ann("PIK3CD", "R/C"))
        assert r.triggered
        assert r.strength == CriterionStrength.SUPPORTING

    def test_gt_rejects_equal_distance(self, tmp_path):
        # RYR1 is strict-greater: candidate R->I (97) NOT > comparator R->F (97).
        db = _db(tmp_path, [_row("1", "RYR1", "NM:p.Arg175Phe", "R175F", 175, "Pathogenic")])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_snv(), _ann("RYR1", "R/I"))
        assert not r.triggered

    def test_gt_passes_when_strictly_greater(self, tmp_path):
        # candidate R->C (180) > comparator R->F (97).
        db = _db(tmp_path, [_row("1", "RYR1", "NM:p.Arg175Phe", "R175F", 175, "Pathogenic")])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_snv(), _ann("RYR1", "R/C"))
        assert r.triggered

    def test_benign_at_codon_blocks(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "PIK3CD", "NM:p.Arg175His", "R175H", 175, "Pathogenic"),
            _row("2", "PIK3CD", "NM:p.Arg175Ser", "R175S", 175, "Benign"),
        ])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_snv(), _ann("PIK3CD", "R/C"))
        assert not r.triggered
        assert "benign variant known" in r.evidence

    def test_candidate_distance_unavailable_withheld(self, tmp_path):
        db = _db(tmp_path, [_row("1", "PIK3CD", "NM:p.Arg175His", "R175H", 175, "Pathogenic")])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_snv(), _ann("PIK3CD", None))
        assert not r.triggered
        assert "Grantham distance unavailable" in r.evidence

    def test_non_gated_gene_keeps_plain_pm5(self, tmp_path):
        # BRCA1 has no Grantham gate: a same-codon different-AA hit meets PM5 at
        # the default Moderate strength regardless of Grantham distances.
        db = _db(tmp_path, [_row("1", "BRCA1", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic")])
        ev = PM5Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_snv(), _ann("BRCA1", "R/H"))  # candidate milder, still met
        assert r.triggered
        assert r.strength == CriterionStrength.MODERATE
        assert "same codon" in r.evidence


class TestRegistryPm5Exclusions:
    def _registry(self):
        return object.__new__(CriteriaRegistry)

    def test_table_lists_runx1_and_dicer1(self):
        assert _PM5_EXCLUSIONS["RUNX1"] == (ACMGCriterion.PM1,)
        assert _PM5_EXCLUSIONS["DICER1"] == (ACMGCriterion.PM1, ACMGCriterion.PS1)

    def test_pm5_suppressed_when_pm1_present(self):
        reg = self._registry()
        results = [
            CriteriaResult.met(ACMGCriterion.PM5),
            CriteriaResult.met(ACMGCriterion.PM1),
        ]
        reg._apply_pm5_exclusions(results, _ann("RUNX1", "R/C"))
        pm5 = next(r for r in results if r.criterion == ACMGCriterion.PM5)
        assert pm5.suppressed
        assert "not with PM1" in pm5.evidence

    def test_dicer1_suppressed_by_ps1(self):
        reg = self._registry()
        results = [
            CriteriaResult.met(ACMGCriterion.PM5),
            CriteriaResult.met(ACMGCriterion.PS1),
        ]
        reg._apply_pm5_exclusions(results, _ann("DICER1", "R/C"))
        assert results[0].suppressed

    def test_pm5_kept_without_clash(self):
        reg = self._registry()
        results = [CriteriaResult.met(ACMGCriterion.PM5)]
        reg._apply_pm5_exclusions(results, _ann("RUNX1", "R/C"))
        assert not results[0].suppressed

    def test_non_listed_gene_untouched(self):
        reg = self._registry()
        results = [
            CriteriaResult.met(ACMGCriterion.PM5),
            CriteriaResult.met(ACMGCriterion.PM1),
        ]
        reg._apply_pm5_exclusions(results, _ann("PIK3CD", "R/C"))
        assert not results[0].suppressed
