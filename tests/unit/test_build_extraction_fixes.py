"""Regression tests for cspec-extraction bug fixes in build_disease_thresholds.py.

Covers:
  * #9  _af_basis must flag "males" only for an in-males FREQUENCY rule, not for
        a "hemizygotes" COUNT clause (OTC/SLC6A8) or prevalence note (ABCD1).
  * #10 BP3 applicability must not be negated by a "not applicable" sub-clause
        when a positive "can be applied to <region>" clause is present (VHL),
        and "AA14-AA48" residue ranges must be parsed.
"""
import importlib.util
from pathlib import Path

_BUILD = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BUILD)
b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b)


def _rs(label: str, desc: str) -> dict:
    """A minimal rule set with one applicable criterion code carrying *desc*."""
    return {
        "criteriaCodes": [{
            "label": label,
            "evidenceStrengths": [
                {"label": "Strong", "applicability": "Applicable", "description": desc},
            ],
        }],
    }


class TestAfBasisMales:
    def test_in_males_frequency_is_males(self):
        rs = _rs("BS1", "Allele frequency in males is greater than 0.0001.")
        assert b._af_basis(rs) == "males"

    def test_hemizygote_count_is_not_males(self):
        # OTC-style: "(male) hemizygotes" is a COUNT rule, not a frequency basis.
        rs = _rs("BA1", "AF above 1.0% Grpmax FAF OR >=10 (female) homozygotes "
                        "or (male) hemizygotes in gnomAD.")
        assert b._af_basis(rs) == ""

    def test_prevalence_in_hemizygotes_is_not_males(self):
        # ABCD1-style: "Total Grpmax FAF" with a prevalence note mentioning
        # hemizygotes — overall FAF, not males.
        rs = _rs("BA1", "Use a Total Grpmax FAF cutoff of >=0.00017 "
                        "(prevalence in hemizygotes is 1 in 5000).")
        assert b._af_basis(rs) == ""


class TestBP3VHL:
    def test_positive_clause_keeps_applicable(self):
        rs = _rs("BP3", "BP3 can be applied to the GXEEX repeat (AA14-AA48). "
                        "Otherwise the rest of the gene has no repeats and BP3 "
                        "is not applicable.")
        assert b._bp3_applicability(rs) == "applicable"

    def test_aa_prefixed_range_parsed(self):
        assert b._bp_residue_ranges("repeat motif (AA14-AA48)") == "14-48"

    def test_plain_not_applicable_still_negated(self):
        rs = _rs("BP3", "BP3 is not applicable for this gene.")
        assert b._bp3_applicability(rs) == "not_applicable"


class TestBS1Strength:
    def test_very_strong_tier_when_no_strong(self):
        # GN023-style: BS1 applicable only at Very Strong + Supporting → the
        # chosen tier (first applicable) is Very Strong, not the Strong default.
        rs = {"criteriaCodes": [{
            "label": "BS1",
            "evidenceStrengths": [
                {"label": "Very Strong", "applicability": "Applicable",
                 "description": "MAF of >=0.003 (0.3%) for autosomal recessive."},
                {"label": "Supporting", "applicability": "Applicable",
                 "description": "MAF of >=0.0007 (0.07%) for autosomal recessive."},
            ],
        }]}
        assert b._bs1_strength(rs) == "VeryStrong"

    def test_strong_tier_default(self):
        rs = {"criteriaCodes": [{
            "label": "BS1",
            "evidenceStrengths": [
                {"label": "Strong", "applicability": "Applicable",
                 "description": "AF >= 0.001"},
            ],
        }]}
        assert b._bs1_strength(rs) == "Strong"

    def test_evaluator_emits_tier_strength(self, tmp_path):
        from acmg_classifier.criteria.allele_frequency import DiseaseThresholds
        from acmg_classifier.models.enums import CriterionStrength
        p = tmp_path / "disease_prevalence.tsv"
        p.write_text(
            "gene_symbol\tbs1_threshold\tbs1_strength\n"
            "MYO15A\t0.003\tVeryStrong\n",
            encoding="utf-8",
        )
        gt = DiseaseThresholds(p).get("MYO15A")
        assert gt.bs1_strength == CriterionStrength.VERY_STRONG


def test_bp1_truncating_includes_lof_classes():
    """#13: RASopathy GoF BP1 truncating target covers splice / start-loss /
    whole-gene deletion, not only nonsense/frameshift."""
    from acmg_classifier.criteria.benign.bp1 import _TRUNCATING
    from acmg_classifier.models.enums import ConsequenceType
    assert ConsequenceType.SPLICE_ACCEPTOR in _TRUNCATING
    assert ConsequenceType.SPLICE_DONOR in _TRUNCATING
    assert ConsequenceType.START_LOST in _TRUNCATING
    assert ConsequenceType.TRANSCRIPT_ABLATION in _TRUNCATING
