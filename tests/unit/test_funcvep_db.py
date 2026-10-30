"""FuncVEP-CTI score reader: TSV and VCF layouts, gene disambiguation."""
from acmg_classifier.local_db.funcvep_db import _match_funcvep, is_vcf

TSV = [
    "17\t100\tG\tA\t0.9500\tBRCA1",
    "17\t100\tG\tA\t0.1000\tNBR2",      # overlapping gene, different score
    "17\t101\tC\tT\t0.2000\tBRCA1",
    "17\t102\tC\tT\t0.3000\tGENE_A",
    "17\t102\tC\tT\t0.3000\tGENE_B",    # same score in both genes
]

VCF = [
    "##fileformat=VCFv4.2",
    "chr17\t100\t.\tG\tA\t.\t.\tCTI=0.9500;GENE=BRCA1",
    "chr17\t100\t.\tG\tA\t.\t.\tCTI=0.1000;GENE=NBR2",
    "chr17\t101\t.\tC\tT\t.\t.\tCTI=0.2000;GENE=BRCA1",
    "chr17\t103\t.\tC\tT,G\t.\t.\tCTI=0.4000;GENE=BRCA1",   # multi-allelic
    "chr17\t104\t.\tC\tT\t.\t.\tGENE=BRCA1",                 # no CTI
    "chr17\t105\t.\tC\tT\t.\t.\tCTI=nan;GENE=BRCA1",
]


class TestIsVcf:
    def test_suffixes(self):
        assert is_vcf("/d/FuncVEP_CTI_hg38.vcf.gz")
        assert is_vcf("/d/x.vcf")
        assert not is_vcf("/d/funcvep_cti_grch38.tsv.gz")


class TestMatchTsv:
    def test_gene_match_wins_over_other_gene(self):
        d = _match_funcvep(TSV, 100, "G", "A", "BRCA1")
        assert (d.score, d.gene, d.match) == (0.95, "BRCA1", "gene")

    def test_conflicting_genes_without_match_abstains(self):
        d = _match_funcvep(TSV, 100, "G", "A", "SOMETHING_ELSE")
        assert d.score is None and d.match == "ambiguous"

    def test_single_row_without_gene(self):
        d = _match_funcvep(TSV, 101, "C", "T", None)
        assert (d.score, d.match) == (0.2, "single")

    def test_same_score_in_several_genes_is_used(self):
        d = _match_funcvep(TSV, 102, "C", "T", None)
        assert (d.score, d.match) == (0.3, "uniform")

    def test_uncovered_position(self):
        assert _match_funcvep(TSV, 999, "C", "T", None) is None

    def test_allele_must_match(self):
        assert _match_funcvep(TSV, 100, "G", "C", "BRCA1") is None


class TestMatchVcf:
    def test_gene_match(self):
        d = _match_funcvep(VCF, 100, "G", "A", "BRCA1", vcf=True)
        assert (d.score, d.gene, d.match) == (0.95, "BRCA1", "gene")

    def test_conflicting_genes_without_match_abstains(self):
        d = _match_funcvep(VCF, 100, "G", "A", "ZZZ", vcf=True)
        assert d.score is None and d.match == "ambiguous"

    def test_single_row(self):
        d = _match_funcvep(VCF, 101, "C", "T", None, vcf=True)
        assert (d.score, d.match) == (0.2, "single")

    def test_multiallelic_record_skipped(self):
        assert _match_funcvep(VCF, 103, "C", "T", None, vcf=True) is None

    def test_record_without_cti_skipped(self):
        assert _match_funcvep(VCF, 104, "C", "T", None, vcf=True) is None

    def test_nan_skipped(self):
        assert _match_funcvep(VCF, 105, "C", "T", None, vcf=True) is None

    def test_header_lines_ignored(self):
        # The '##' header must not be parsed as a record.
        assert _match_funcvep(VCF[:1], 100, "G", "A", None, vcf=True) is None


class TestIndexedFiles:
    """Both on-disk layouts must return the same thing through tabix, including
    the chr-prefix fallback (the VCF is chr-prefixed, the TSV is not)."""

    @staticmethod
    def _build(tmp_path):
        import pysam
        vcf = tmp_path / "FuncVEP_CTI_hg38.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            '##INFO=<ID=CTI,Number=1,Type=Float,Description="FuncVEP-CTI score">\n'
            '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr17\t43044295\t.\tG\tA\t.\t.\tCTI=0.9500;GENE=BRCA1\n"
            "chr17\t43044295\t.\tG\tA\t.\t.\tCTI=0.1000;GENE=NBR2\n"
            "chr17\t43045000\t.\tC\tT\t.\t.\tCTI=0.0150;GENE=BRCA1\n",
            encoding="utf-8",
        )
        pysam.tabix_compress(str(vcf), str(vcf) + ".gz", force=True)
        pysam.tabix_index(str(vcf) + ".gz", preset="vcf", force=True)

        tsv = tmp_path / "funcvep_cti_grch38.tsv"
        tsv.write_text(
            "17\t43044295\tG\tA\t0.9500\tBRCA1\n"
            "17\t43044295\tG\tA\t0.1000\tNBR2\n"
            "17\t43045000\tC\tT\t0.0150\tBRCA1\n",
            encoding="utf-8",
        )
        pysam.tabix_compress(str(tsv), str(tsv) + ".gz", force=True)
        pysam.tabix_index(str(tsv) + ".gz", seq_col=0, start_col=1, end_col=1, force=True)
        from pathlib import Path
        return Path(str(vcf) + ".gz"), Path(str(tsv) + ".gz")

    def test_vcf_and_tsv_agree(self, tmp_path):
        from acmg_classifier.local_db.funcvep_db import query_funcvep
        vcf_gz, tsv_gz = self._build(tmp_path)
        for path in (vcf_gz, tsv_gz):
            for chrom in ("chr17", "17"):   # chr-prefix fallback both ways
                d = query_funcvep(path, chrom, 43044295, "G", "A", gene="BRCA1")
                assert (d.score, d.gene, d.match) == (0.95, "BRCA1", "gene"), path
            d = query_funcvep(path, "17", 43045000, "C", "T", gene="BRCA1")
            assert d.score == 0.015
            assert query_funcvep(path, "17", 1, "G", "A", gene="BRCA1") is None

    def test_reader_path_matches(self, tmp_path):
        from acmg_classifier.local_db.funcvep_db import query_funcvep_reader
        from acmg_classifier.utils.tabix import TabixReader
        vcf_gz, tsv_gz = self._build(tmp_path)
        for path in (vcf_gz, tsv_gz):
            d = query_funcvep_reader(TabixReader(path), "17", 43044295, "G", "A", gene="NBR2")
            assert (d.score, d.gene) == (0.1, "NBR2"), path

    def test_missing_file_returns_none(self, tmp_path):
        from acmg_classifier.local_db.funcvep_db import query_funcvep
        assert query_funcvep(tmp_path / "nope.vcf.gz", "17", 1, "G", "A") is None
