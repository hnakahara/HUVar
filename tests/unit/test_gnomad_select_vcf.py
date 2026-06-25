"""_select_vcf_files: callset + per-chromosome filtering for the gnomAD builders."""
from acmg_classifier.setup.gnomad_builder import _select_vcf_files


def _make(d, names):
    for n in names:
        (d / n).write_text("x")


def _contigs(paths):
    return sorted(p.name.split(".sites.")[1] for p in paths)


def test_chromosome_filter_restricts_to_requested(tmp_path):
    _make(tmp_path, [
        "gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz",
        "gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz",
        "gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz",
        "gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz",
    ])
    out = _select_vcf_files(tmp_path, ("genomes",), ("chr13", "chr17"))
    assert _contigs(out) == ["chr13.vcf.bgz", "chr17.vcf.bgz"]


def test_chromosome_filter_no_prefix_false_match(tmp_path):
    # chr1 must NOT match chr13 / chr17 (dot-bounded token).
    _make(tmp_path, [
        "gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz",
        "gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz",
    ])
    out = _select_vcf_files(tmp_path, ("genomes",), ("chr1",))
    assert _contigs(out) == ["chr1.vcf.bgz"]


def test_callset_filter_excludes_other_callset(tmp_path):
    _make(tmp_path, [
        "gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz",
        "gnomad.exomes.v4.1.sites.chr13.vcf.bgz",
    ])
    out = _select_vcf_files(tmp_path, ("genomes",), ("chr13",))
    assert [p.name for p in out] == ["gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz"]


def test_no_chromosome_filter_returns_all_genomes(tmp_path):
    _make(tmp_path, [
        "gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz",
        "gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz",
    ])
    out = _select_vcf_files(tmp_path, ("genomes",))
    assert _contigs(out) == ["chr13.vcf.bgz", "chr17.vcf.bgz"]


def test_chromosome_filter_no_match_returns_empty(tmp_path):
    _make(tmp_path, ["gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz"])
    assert _select_vcf_files(tmp_path, ("genomes",), ("chr17",)) == []
