"""PS1/PM5 comparators restricted to VCEP classifications (GM revision, R3 comment 3)."""
from acmg_classifier.criteria.vcep_comparators import (
    filter_comparators, load_erepo, load_rules, requires_vcep_comparator,
)
from acmg_classifier.models.annotation import ClinVarRecord


def _rec(chrom, pos, ref, alt, stars=1, sig="Pathogenic"):
    return ClinVarRecord(variation_id="1", clinical_significance=sig, review_status="x",
                         star_rating=stars, chrom=chrom, pos=pos, ref=ref, alt=alt)


def test_rules_loaded():
    rules = load_rules()
    assert requires_vcep_comparator("MLH1", "PS1").protein_level
    assert "PS1" in requires_vcep_comparator("CTLA4", "PS1").exclude_codes
    assert requires_vcep_comparator("BRCA1", "PS1") is None


def test_non_rule_gene_unchanged():
    hits = [_rec("17", 1, "A", "G")]
    kept, note = filter_comparators("BRCA1", "PS1", hits, "GRCh38")
    assert kept == hits and note is None


def test_vcep_gene_requires_expert_or_erepo():
    erepo = load_erepo("GRCh38")
    # pick an eRepo P/LP TP53 variant without PVS1/PM5 in its codes
    key = next(k for k, (a, c) in erepo.items() if "Pathogenic" in a and "PM5" not in c)
    chrom, pos, ref, alt = key.split(":")
    in_erepo = _rec(chrom, int(pos), ref, alt, stars=1)
    one_star = _rec("17", 999999999, "A", "G", stars=1)
    three_star = _rec("17", 999999998, "A", "G", stars=3)
    kept, note = filter_comparators("TP53", "PM5", [in_erepo, one_star, three_star], "GRCh38")
    assert in_erepo in kept and three_star in kept and one_star not in kept
    assert "2/3" in note


def test_exclude_codes_and_protein_level():
    erepo = load_erepo("GRCh38")
    pm5_key = next((k for k, (a, c) in erepo.items() if "Pathogenic" in a and "PM5" in c), None)
    if pm5_key:
        c, p, r, a = pm5_key.split(":")
        kept, _ = filter_comparators("PIK3R1", "PM5", [_rec(c, int(p), r, a)], "GRCh38")
        assert kept == []
    pvs1_key = next((k for k, (a, c) in erepo.items() if "Pathogenic" in a and "PVS1" in c), None)
    if pvs1_key:
        c, p, r, a = pvs1_key.split(":")
        kept, _ = filter_comparators("MLH1", "PS1", [_rec(c, int(p), r, a)], "GRCh38")
        assert kept == []
