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


# ---------------------------------------------------------------------------
# VCEP classification in the batch (VCF) output
# ---------------------------------------------------------------------------

import csv as _csv

import pytest as _pytest


@_pytest.mark.parametrize("vcep,auto,exp", [
    ("Pathogenic", "Pathogenic", "exact"),
    ("Likely Pathogenic", "Pathogenic", "same_direction"),
    ("Likely Benign", "Benign", "same_direction"),
    ("Uncertain Significance", "VUS", "exact"),
    ("Uncertain Significance", "LikelyPathogenic", "discordant"),
    ("Pathogenic", "Benign", "discordant"),
    ("Likely Pathogenic|Pathogenic", "Pathogenic", "vcep_conflict"),
    ("", "Pathogenic", ""),
])
def test_concordance(vcep, auto, exp):
    from acmg_classifier.local_db.erepo_db import concordance
    assert concordance(vcep, auto) == exp


def test_vcep_record_and_tsv_columns(tmp_path):
    from acmg_classifier.io.tsv_writer import _HEADER, _VCEP_HEADER, write_tsv
    from acmg_classifier.local_db.erepo_db import vcep_record
    from acmg_classifier.models.classification import ClassificationResult
    from acmg_classifier.models.enums import Pathogenicity

    rec = vcep_record("GRCh38", "chr17", 43045709, "A", "G")   # ENIGMA: Likely Pathogenic
    assert rec is not None and rec.assertion == "Likely Pathogenic"
    assert "ENIGMA" in rec.expert_panel and rec.erepo_url.startswith("https://")
    assert vcep_record("GRCh38", "chr17", 43045709, "A", "T") is None

    def res(vid, vcep):
        return ClassificationResult(
            variant_id=vid, chrom="chr17", pos=43045709, ref="A", alt=vid[-1],
            criteria_results=[], classification_2015=Pathogenicity.PATHOGENIC,
            classification_2015_rules="", bayesian_score=6,
            classification_bayesian=Pathogenicity.LIKELY_PATHOGENIC, vcep=vcep)

    out = tmp_path / "out.tsv"
    write_tsv([res("chr17:43045709:A:G", rec), res("chr17:43045709:A:T", None)], out)
    rows = list(_csv.DictReader(out.open(encoding="utf-8"), delimiter="\t"))
    assert all(c in rows[0] for c in _VCEP_HEADER)
    assert len(rows[0]) == len(_HEADER)
    assert rows[0]["vcep_classification"] == "Likely Pathogenic"
    assert rows[0]["vcep_expert_panel"].startswith("ENIGMA")
    assert rows[0]["vcep_concordance_2015"] == "same_direction"
    assert rows[0]["vcep_concordance_bayesian"] == "exact"
    assert rows[1]["vcep_classification"] == "" and rows[1]["vcep_concordance_bayesian"] == ""
