"""Unit tests for supplement TSV reader."""
import pytest
from pathlib import Path
from acmg_classifier.exceptions import SupplementParseError
from acmg_classifier.io.supplement_reader import parse_inline_evidence, read_supplement
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength

FIXTURE_TSV = Path(__file__).parent.parent / "fixtures" / "sample_supplement.tsv"


def test_read_supplement_parses_correctly():
    result = read_supplement(FIXTURE_TSV)
    assert "chr17:43044295:G:A" in result
    entries = result["chr17:43044295:G:A"]
    assert len(entries) == 2
    ps3 = next(e for e in entries if e.criterion == ACMGCriterion.PS3)
    assert ps3.strength == CriterionStrength.STRONG
    assert "PMID" in ps3.evidence


class TestParseInlineEvidence:
    def test_basic_pair(self):
        e = parse_inline_evidence(["PS3:strong"], "chr1:100:A:T")[0]
        assert e.criterion == ACMGCriterion.PS3
        assert e.strength == CriterionStrength.STRONG
        assert e.variant_id == "chr1:100:A:T"
        assert e.source == "explain --evidence"

    def test_with_note_containing_colon(self):
        e = parse_inline_evidence(["PM1:moderate:PMID:12345 hotspot"], "k")[0]
        assert e.criterion == ACMGCriterion.PM1
        assert e.strength == CriterionStrength.MODERATE
        assert e.evidence == "PMID:12345 hotspot"

    def test_case_insensitive_and_aliases(self):
        out = parse_inline_evidence(
            ["pvs1:VeryStrong", "ba1:stand_alone", "pm2:supp", "ps4:mod"], "k"
        )
        assert [o.criterion for o in out] == [
            ACMGCriterion.PVS1, ACMGCriterion.BA1, ACMGCriterion.PM2, ACMGCriterion.PS4]
        assert [o.strength for o in out] == [
            CriterionStrength.VERY_STRONG, CriterionStrength.VERY_STRONG,
            CriterionStrength.SUPPORTING, CriterionStrength.MODERATE]

    def test_invalid_format_raises(self):
        with pytest.raises(SupplementParseError):
            parse_inline_evidence(["PS3"], "k")

    def test_unknown_criterion_raises(self):
        with pytest.raises(SupplementParseError):
            parse_inline_evidence(["ZZ9:strong"], "k")

    def test_unknown_strength_raises(self):
        with pytest.raises(SupplementParseError):
            parse_inline_evidence(["PS3:bogus"], "k")
