"""VCEP classification lookup (Genome Medicine revision, Reviewer 3, comment 20)."""
from acmg_classifier.local_db.erepo_db import lookup, vcep_warning


def test_lookup_known_erepo_variant():
    rec = lookup("GRCh38", "chr10", 110964362, "A", "G")  # SHOC2 c.4A>G
    assert rec is not None
    assert rec["assertion"] == "Pathogenic"
    assert rec["expert_panel"] and rec["erepo_url"].startswith("https://erepo")
    assert rec["snapshot"] == "2026-05-28"
    # chrom prefix does not matter
    assert lookup("GRCh38", "10", 110964362, "A", "G") == rec


def test_unknown_variant_has_no_record_or_warning():
    assert lookup("GRCh38", "chr1", 1, "A", "G") is None
    assert vcep_warning("GRCh38", "chr1", 1, "A", "G") == []


def test_warning_text():
    w = vcep_warning("GRCh38", "chr10", 110964362, "A", "G")
    assert len(w) == 1 and w[0].startswith("VCEP_CLASSIFIED:")
    assert "Pathogenic" in w[0] and "precedence" in w[0]
