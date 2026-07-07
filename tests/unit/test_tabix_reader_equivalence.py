"""Persistent-handle tabix readers must equal the connection-per-call queries.

Phase 2 of the annotation speed-up replaces `pysam.TabixFile` open-per-variant
with a thread-local `TabixReader`. Each `query_*_reader` must therefore return
exactly what the legacy `query_*` returns. These tests build real bgzipped,
tabix-indexed fixtures and assert equality across present / absent / wrong-allele
cases and the chr-prefix fallback.

Requires pysam (skipped where it is not installed).
"""
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

from acmg_classifier.utils.tabix import TabixReader  # noqa: E402
from acmg_classifier.local_db.revel_db import query_revel, query_revel_reader  # noqa: E402
from acmg_classifier.local_db.bayesdel_db import (  # noqa: E402
    query_bayesdel, query_bayesdel_reader,
)
from acmg_classifier.local_db.cadd_db import query_cadd, query_cadd_reader  # noqa: E402
from acmg_classifier.local_db.alphamissense_db import (  # noqa: E402
    query_alphamissense, query_alphamissense_reader,
)
from acmg_classifier.local_db.clinvar_vcf import (  # noqa: E402
    query_clinvar_vcf, query_clinvar_vcf_reader,
)
from acmg_classifier.local_db.repeatmasker_db import (  # noqa: E402
    query_repeat, query_repeat_reader,
)


def _bgzip_index(path: Path, lines: list[str], *, preset=None,
                 seq_col=0, start_col=1, end_col=1) -> Path:
    """Write *lines* to a bgzipped, tabix-indexed file and return its .gz path."""
    raw = path
    raw.write_text(
        "".join(ln if ln.endswith("\n") else ln + "\n" for ln in lines),
        encoding="utf-8",
    )
    gz = Path(str(raw) + ".gz")
    pysam.tabix_compress(str(raw), str(gz), force=True)
    if preset:
        pysam.tabix_index(str(gz), preset=preset, force=True)
    else:
        pysam.tabix_index(str(gz), seq_col=seq_col, start_col=start_col,
                          end_col=end_col, force=True)
    return gz


# 5-column (chrom pos ref alt score) TSV shared by REVEL / BayesDel / CADD.
_SCORE_ROWS = [
    "chr1\t100\tA\tG\t0.750",
    "chr1\t100\tA\tT\t0.120",
    "chr1\t200\tC\tT\t0.900",
]
_SCORE_CASES = [
    ("chr1", 100, "A", "G"),   # present
    ("chr1", 100, "A", "T"),   # same pos, different alt
    ("chr1", 100, "A", "C"),   # uncovered alt → None
    ("chr1", 300, "G", "A"),   # absent → None
    ("1", 100, "A", "G"),      # bare chrom → chr-prefix fallback
]


@pytest.mark.parametrize("chrom,pos,ref,alt", _SCORE_CASES)
def test_revel_reader_matches_query(tmp_path, chrom, pos, ref, alt):
    gz = _bgzip_index(tmp_path / "revel.tsv", _SCORE_ROWS)
    assert query_revel_reader(TabixReader(gz), chrom, pos, ref, alt) == \
        query_revel(gz, chrom, pos, ref, alt)


@pytest.mark.parametrize("chrom,pos,ref,alt", _SCORE_CASES)
def test_bayesdel_reader_matches_query(tmp_path, chrom, pos, ref, alt):
    gz = _bgzip_index(tmp_path / "bayesdel.tsv", _SCORE_ROWS)
    assert query_bayesdel_reader(TabixReader(gz), chrom, pos, ref, alt) == \
        query_bayesdel(gz, chrom, pos, ref, alt)


@pytest.mark.parametrize("chrom,pos,ref,alt", _SCORE_CASES)
def test_cadd_reader_matches_query(tmp_path, chrom, pos, ref, alt):
    gz = _bgzip_index(tmp_path / "cadd.tsv", _SCORE_ROWS)
    assert query_cadd_reader(TabixReader(gz), chrom, pos, ref, alt) == \
        query_cadd(gz, chrom, pos, ref, alt)


_AM_ROWS = [
    "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript\tprotein_variant\tam_pathogenicity\tam_class",
    "chr1\t100\tA\tG\thg38\tP1\tENST1\tM1V\t0.980\tlikely_pathogenic",
    "chr1\t100\tA\tT\thg38\tP1\tENST1\tM1L\t0.100\tlikely_benign",
]


@pytest.mark.parametrize("chrom,pos,ref,alt", [
    ("chr1", 100, "A", "G"), ("chr1", 100, "A", "T"),
    ("chr1", 100, "A", "C"), ("1", 100, "A", "G"),
])
def test_alphamissense_reader_matches_query(tmp_path, chrom, pos, ref, alt):
    gz = _bgzip_index(tmp_path / "am.tsv", _AM_ROWS)
    assert query_alphamissense_reader(TabixReader(gz), chrom, pos, ref, alt) == \
        query_alphamissense(gz, chrom, pos, ref, alt)


_VCF_ROWS = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t100\t12345\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter",
    "chr1\t100\t.\tA\tT\t.\t.\tCLNSIG=Likely_benign;CLNREVSTAT=no_assertion_criteria_provided",
]


@pytest.mark.parametrize("chrom,pos,ref,alt", [
    ("chr1", 100, "A", "G"), ("chr1", 100, "A", "T"),
    ("chr1", 100, "A", "C"), ("1", 100, "A", "G"),
])
def test_clinvar_reader_matches_query(tmp_path, chrom, pos, ref, alt):
    gz = _bgzip_index(tmp_path / "clinvar.vcf", _VCF_ROWS, preset="vcf")
    assert query_clinvar_vcf_reader(TabixReader(gz), chrom, pos, ref, alt) == \
        query_clinvar_vcf(gz, chrom, pos, ref, alt)


_BED_ROWS = [
    "chr1\t50\t150\tAluY\t0\t+\tSINE\tAlu",
    "chr1\t400\t500\tL1\t0\t-\tLINE\tL1",
]


@pytest.mark.parametrize("chrom,pos", [
    ("chr1", 100), ("chr1", 450), ("chr1", 1000), ("1", 100),
])
def test_repeat_reader_matches_query(tmp_path, chrom, pos):
    gz = _bgzip_index(tmp_path / "repeat.bed", _BED_ROWS, preset="bed")
    assert query_repeat_reader(TabixReader(gz), chrom, pos) == \
        query_repeat(gz, chrom, pos)


def test_missing_file_readers_match_query(tmp_path):
    absent = tmp_path / "nope.tsv.gz"
    r = TabixReader(absent)
    assert query_revel_reader(r, "chr1", 100, "A", "G") is None
    assert query_bayesdel_reader(r, "chr1", 100, "A", "G") is None
    assert query_cadd_reader(r, "chr1", 100, "A", "G") is None
    assert query_alphamissense_reader(r, "chr1", 100, "A", "G") is None
    assert query_clinvar_vcf_reader(r, "chr1", 100, "A", "G") == []
    assert query_repeat_reader(r, "chr1", 100).in_repeat is False


def test_reader_reuses_handle_across_calls(tmp_path):
    gz = _bgzip_index(tmp_path / "revel.tsv", _SCORE_ROWS)
    reader = TabixReader(gz)
    a = query_revel_reader(reader, "chr1", 100, "A", "G")
    b = query_revel_reader(reader, "chr1", 200, "C", "T")
    assert a.score == 0.750 and b.score == 0.900
