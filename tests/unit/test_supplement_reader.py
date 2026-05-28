"""Unit tests for supplement TSV reader."""
from pathlib import Path
from acmg_classifier.io.supplement_reader import read_supplement
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
