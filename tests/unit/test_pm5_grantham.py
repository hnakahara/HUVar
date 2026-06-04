"""PM5 precision controls: Grantham gate, comparator quality, strength, exclusions.

A subset of ClinGen VCEPs require PM5 to clear a Grantham-distance test against
the same-codon comparator; all VCEPs anchor PM5 on an *established* comparator
and several forbid combining PM5 with PM1/PS1 or cap it at Supporting. These
tests cover the embedded Grantham 1974 matrix, the ``pm5_*`` cspec extraction,
the per-gene loader, the comparator-quality / benign / strength gates, and the
registry exclusion pass.
"""
import importlib.util
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.grantham import grantham_distance
from acmg_classifier.criteria.pm5_genes import PM5Spec
from acmg_classifier.criteria.pathogenic.pm5 import PM5Evaluator
from acmg_classifier.criteria.registry import CriteriaRegistry
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


# --------------------------- cspec extraction --------------------------------

def _pm5_code(*strengths):
    return {"label": "PM5", "evidenceStrengths": [
        {"label": lbl, "applicability": "Applicable", "description": desc}
        for lbl, desc in strengths
    ]}


class TestPm5GranthamExtraction:
    def test_inclusive_ge(self):
        rs = {"criteriaCodes": [_pm5_code(
            ("Moderate", "Grantham distance greater than or equal to the variant."),
        )]}
        assert bdt._pm5_grantham_op(rs) == "ge"

    def test_strict_gt_higher(self):
        rs = {"criteriaCodes": [_pm5_code(
            ("Moderate", "must have a higher Grantham score than the comparison variant."),
        )]}
        assert bdt._pm5_grantham_op(rs) == "gt"

    def test_strict_gt_less_than(self):
        rs = {"criteriaCodes": [_pm5_code(
            ("Moderate", "Grantham score for alternate variant must be less than for "
                         "variant being assessed."),
        )]}
        assert bdt._pm5_grantham_op(rs) == "gt"

    def test_ge_wins_when_mixed(self):
        rs = {"criteriaCodes": [_pm5_code(
            ("Moderate", "Grantham distance greater than or equal to the variant."),
            ("Supporting", "pathogenic but has a greater Grantham distance."),
        )]}
        assert bdt._pm5_grantham_op(rs) == "ge"

    def test_no_grantham_is_empty(self):
        rs = {"criteriaCodes": [_pm5_code(
            ("Moderate", "Novel missense at a residue with a known pathogenic change."),
        )]}
        assert bdt._pm5_grantham_op(rs) == ""


class TestPm5ExcludesExtraction:
    def test_pm1_exclusion(self):
        rs = {"criteriaCodes": [_pm5_code(("Moderate", "PM5 should not be combined with PM1."))]}
        assert bdt._pm5_excludes(rs) == "PM1"

    def test_pm1_if_applied(self):
        rs = {"criteriaCodes": [_pm5_code(("Moderate", "PM5 cannot be used if PM1 was applied."))]}
        assert bdt._pm5_excludes(rs) == "PM1"

    def test_pm1_and_ps1(self):
        rs = {"criteriaCodes": [_pm5_code(
            ("Moderate", "This rule cannot be applied in combination with PM1 or PS1."),
        )]}
        assert bdt._pm5_excludes(rs) == "PM1,PS1"

    def test_not_mutually_exclusive_is_not_an_exclusion(self):
        rs = {"criteriaCodes": [_pm5_code(("Moderate", "Not mutually exclusive with PM1."))]}
        assert bdt._pm5_excludes(rs) == ""

    def test_no_pm1_mention(self):
        rs = {"criteriaCodes": [_pm5_code(("Moderate", "Different missense at the residue."))]}
        assert bdt._pm5_excludes(rs) == ""


class TestPm5MaxExtraction:
    def test_supporting_only(self):
        rs = {"criteriaCodes": [_pm5_code(("Supporting", "PM5_Supporting for an LP comparator."))]}
        assert bdt._pm5_max(rs) == "Supporting"

    def test_moderate_present(self):
        rs = {"criteriaCodes": [_pm5_code(
            ("Supporting", "..."), ("Moderate", "..."),
        )]}
        assert bdt._pm5_max(rs) == "Moderate"

    def test_no_applicable_pm5(self):
        rs = {"criteriaCodes": [{"label": "PM5", "evidenceStrengths": [
            {"label": "Moderate", "applicability": "Not Applicable for this VCEP"},
        ]}]}
        assert bdt._pm5_max(rs) == ""


# ------------------------------ loader ---------------------------------------

class TestPm5SpecLoader:
    def _tsv(self, tmp_path):
        tsv = tmp_path / "dp.tsv"
        tsv.write_text(
            "gene_symbol\tpm5_grantham\tpm5_excludes\tpm5_max\n"
            "PIK3CD\tge\t\t\n"
            "PIK3R1\tgt\t\t\n"
            "RUNX1\tge\tPM1\t\n"
            "DICER1\tge\tPM1,PS1\t\n"
            "ATM\t\t\tSupporting\n"
            "MYH7\t\tPM1\t\n",
            encoding="utf-8",
        )
        return tsv

    def test_operators(self, tmp_path):
        s = PM5Spec(self._tsv(tmp_path))
        assert s.operator("PIK3CD") == "ge"
        assert s.operator("PIK3R1") == "gt"
        assert s.operator("MYH7") == ""
        assert s.operator("UNSEEN") == ""
        assert s.operator(None) == ""

    def test_excludes(self, tmp_path):
        s = PM5Spec(self._tsv(tmp_path))
        assert s.excludes("RUNX1") == (ACMGCriterion.PM1,)
        assert s.excludes("DICER1") == (ACMGCriterion.PM1, ACMGCriterion.PS1)
        assert s.excludes("MYH7") == (ACMGCriterion.PM1,)
        assert s.excludes("PIK3CD") == ()
        assert s.excludes(None) == ()

    def test_max_strength(self, tmp_path):
        s = PM5Spec(self._tsv(tmp_path))
        assert s.max_strength("ATM") == CriterionStrength.SUPPORTING
        assert s.max_strength("MYH7") is None
        assert s.max_strength("UNSEEN") is None

    def test_missing_file_is_empty(self, tmp_path):
        s = PM5Spec(tmp_path / "nope.tsv")
        assert s.operator("PIK3CD") == ""
        assert s.excludes("RUNX1") == ()
        assert s.max_strength("ATM") is None


# --------------------------- ClinVar fixtures --------------------------------

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


def _row(vid, gene, hgvs_p, aa, codon, sig, stars=2):
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

    def test_star_threshold(self, tmp_path):
        db = _db(tmp_path, [_row("1", "G", "NM:p.Arg175Ser", "R175S", 175, "Benign", stars=1)])
        assert has_benign_at_codon(db, "G", 175, min_stars=2) is False
        assert has_benign_at_codon(db, "G", 175, min_stars=1) is True


# ------------------------------ evaluator ------------------------------------

def _cfg(tmp_path, clinvar: Path, min_stars: int = 1):
    tsv = tmp_path / "disease_prevalence.tsv"
    tsv.write_text(
        "gene_symbol\tpm5_grantham\tpm5_excludes\tpm5_max\n"
        "PIK3CD\tge\t\t\n"
        "RYR1\tgt\t\t\n"
        "ATM\t\t\tSupporting\n"
        "BRCA1\t\t\t\n",
        encoding="utf-8",
    )
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = clinvar
    cfg.pm5_min_stars = min_stars
    return cfg


_AA1_TO_AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}


def _consequence(gene, amino_acids, codon=175, hgvs_p=None):
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
        db = _db(tmp_path, [_row("1", "PIK3CD", "NM:p.Arg175His", "R175H", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("PIK3CD", "R/C"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE
        assert "Grantham-gated" in r.evidence

    def test_ge_fail_when_candidate_milder(self, tmp_path):
        db = _db(tmp_path, [_row("1", "PIK3CD", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("PIK3CD", "R/H"))
        assert not r.triggered and "Grantham gate failed" in r.evidence

    def test_likely_pathogenic_comparator_is_supporting(self, tmp_path):
        db = _db(tmp_path, [_row("1", "PIK3CD", "NM:p.Arg175His", "R175H", 175, "Likely pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("PIK3CD", "R/C"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_gt_rejects_equal_distance(self, tmp_path):
        db = _db(tmp_path, [_row("1", "RYR1", "NM:p.Arg175Phe", "R175F", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("RYR1", "R/I"))
        assert not r.triggered

    def test_gt_passes_when_strictly_greater(self, tmp_path):
        db = _db(tmp_path, [_row("1", "RYR1", "NM:p.Arg175Phe", "R175F", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("RYR1", "R/C"))
        assert r.triggered

    def test_candidate_distance_unavailable_withheld(self, tmp_path):
        db = _db(tmp_path, [_row("1", "PIK3CD", "NM:p.Arg175His", "R175H", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("PIK3CD", None))
        assert not r.triggered and "unavailable" in r.evidence


class TestPm5ComparatorQuality:
    def test_default_min_stars_is_one(self):
        # Default keeps single-submitter (1-star) P/LP comparators — often the
        # only legitimate PM5 anchor. ACMG_PM5_MIN_STARS=2 opts into the stricter
        # expert/multi-submitter-only policy.
        from acmg_classifier.config import Config
        assert Config().pm5_min_stars == 1

    def test_benign_at_codon_blocks_any_gene(self, tmp_path):
        db = _db(tmp_path, [
            _row("1", "BRCA1", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic"),
            _row("2", "BRCA1", "NM:p.Arg175Ser", "R175S", 175, "Benign"),
        ])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("BRCA1", "R/H"))
        assert not r.triggered and "Benign variant known" in r.evidence

    def test_min_stars_excludes_single_submitter(self, tmp_path):
        # A 1-star comparator does not anchor PM5 at the default min_stars=2.
        db = _db(tmp_path, [_row("1", "BRCA1", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic", stars=1)])
        r = PM5Evaluator(_cfg(tmp_path, db, min_stars=2)).evaluate(_snv(), _ann("BRCA1", "R/H"))
        assert not r.triggered and "star" in r.evidence

    def test_min_stars_override_admits_single_submitter(self, tmp_path):
        db = _db(tmp_path, [_row("1", "BRCA1", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic", stars=1)])
        r = PM5Evaluator(_cfg(tmp_path, db, min_stars=1)).evaluate(_snv(), _ann("BRCA1", "R/H"))
        assert r.triggered


class TestPm5StrengthForPlainGenes:
    def test_pathogenic_comparator_moderate(self, tmp_path):
        db = _db(tmp_path, [_row("1", "BRCA1", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("BRCA1", "R/H"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE
        assert "same codon" in r.evidence

    def test_likely_pathogenic_only_drops_to_supporting(self, tmp_path):
        db = _db(tmp_path, [_row("1", "BRCA1", "NM:p.Arg175Cys", "R175C", 175, "Likely pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("BRCA1", "R/H"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_supporting_only_gene_capped(self, tmp_path):
        # ATM's VCEP caps PM5 at Supporting even with a Pathogenic comparator.
        db = _db(tmp_path, [_row("1", "ATM", "NM:p.Arg175Cys", "R175C", 175, "Pathogenic")])
        r = PM5Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("ATM", "R/H"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING


class TestRegistryPm5Exclusions:
    def _registry(self, tmp_path):
        tsv = tmp_path / "dp.tsv"
        tsv.write_text(
            "gene_symbol\tpm5_excludes\n"
            "RUNX1\tPM1\n"
            "DICER1\tPM1,PS1\n"
            "MYH7\tPM1\n"
            "PIK3CD\t\n",
            encoding="utf-8",
        )
        reg = object.__new__(CriteriaRegistry)
        reg._pm5 = PM5Spec(tsv)
        return reg

    def test_pm5_suppressed_when_pm1_present(self, tmp_path):
        reg = self._registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PM5),
            CriteriaResult.met(ACMGCriterion.PM1),
        ]
        reg._apply_pm5_exclusions(results, _ann("RUNX1", "R/C"))
        pm5 = next(r for r in results if r.criterion == ACMGCriterion.PM5)
        assert pm5.suppressed and "not with PM1" in pm5.evidence

    def test_rasopathy_gene_generalised(self, tmp_path):
        # MYH7 (Cardiomyopathy VCEP) also forbids PM5+PM1 — not hardcoded.
        reg = self._registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PM5),
            CriteriaResult.met(ACMGCriterion.PM1),
        ]
        reg._apply_pm5_exclusions(results, _ann("MYH7", "R/C"))
        assert results[0].suppressed

    def test_dicer1_suppressed_by_ps1(self, tmp_path):
        reg = self._registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PM5),
            CriteriaResult.met(ACMGCriterion.PS1),
        ]
        reg._apply_pm5_exclusions(results, _ann("DICER1", "R/C"))
        assert results[0].suppressed

    def test_pm5_kept_without_clash(self, tmp_path):
        reg = self._registry(tmp_path)
        results = [CriteriaResult.met(ACMGCriterion.PM5)]
        reg._apply_pm5_exclusions(results, _ann("RUNX1", "R/C"))
        assert not results[0].suppressed

    def test_non_listed_gene_untouched(self, tmp_path):
        reg = self._registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PM5),
            CriteriaResult.met(ACMGCriterion.PM1),
        ]
        reg._apply_pm5_exclusions(results, _ann("PIK3CD", "R/C"))
        assert not results[0].suppressed
