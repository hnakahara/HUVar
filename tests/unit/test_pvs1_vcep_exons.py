"""Exon-aware VCEP PVS1 canonical-splice refinement (opt-in override table)."""
from pathlib import Path

from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType, CriterionStrength
from acmg_classifier.pvs1.vcep_pvs1 import evaluate_vcep_pvs1
from acmg_classifier.pvs1.vcep_pvs1_exons import SpliceExonOverrides, skipped_exon

S = CriterionStrength


def _pc(gene, consequence, intron=None):
    return ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=consequence, biotype="protein_coding", intron=intron,
    )


class TestSkippedExon:
    def test_donor_skips_same_intron_exon(self):
        assert skipped_exon(ConsequenceType.SPLICE_DONOR, "10/26") == 10

    def test_acceptor_skips_next_exon(self):
        assert skipped_exon(ConsequenceType.SPLICE_ACCEPTOR, "9/26") == 10

    def test_unparseable_returns_none(self):
        assert skipped_exon(ConsequenceType.SPLICE_DONOR, None) is None
        assert skipped_exon(ConsequenceType.SPLICE_DONOR, "x") is None


class TestOverridesLoader:
    def test_missing_file_is_empty(self):
        ov = SpliceExonOverrides(Path("does/not/exist.tsv"))
        assert not ov
        assert ov.lookup("DICER1", 10) is None

    def test_leading_comment_lines_skipped(self, tmp_path):
        # csv.DictReader must use the real header, not a leading '#' comment.
        p = tmp_path / "ov.tsv"
        p.write_text(
            "# a comment before the header\n"
            "# another comment\n"
            "gene\tskipped_exon\tstrength\n"
            "DICER1\t10\tstrong\n",
            encoding="utf-8",
        )
        ov = SpliceExonOverrides(p)
        assert ov.lookup("DICER1", 10) == S.STRONG

    def test_shipped_table_loads(self):
        # The committed override table should parse to a non-empty set with the
        # documented DICER1 entries.
        ov = SpliceExonOverrides(Path("data/shared/vcep_pvs1_splice_exons.tsv"))
        assert ov
        assert ov.lookup("DICER1", 10) == S.STRONG
        assert ov.lookup("DICER1", 5) == S.MODERATE
        assert ov.lookup("GAA", 20) == S.MODERATE

    def test_load_and_lookup(self, tmp_path):
        p = tmp_path / "ov.tsv"
        p.write_text(
            "gene\tskipped_exon\tstrength\tnote\n"
            "# DICER1\t99\tna\tcommented out\n"
            "DICER1\t10\tstrong\tin-frame >10%\n"
            "DICER1\t5\tmoderate\tin-frame <10%\n",
            encoding="utf-8",
        )
        ov = SpliceExonOverrides(p)
        assert ov
        assert ov.lookup("DICER1", 10) == S.STRONG
        assert ov.lookup("DICER1", 5) == S.MODERATE
        assert ov.lookup("DICER1", 99) is None  # commented row skipped
        assert ov.lookup("DICER1", 7) is None   # no entry


class TestSpliceRefinement:
    def _ov(self, tmp_path):
        p = tmp_path / "ov.tsv"
        p.write_text(
            "gene\tskipped_exon\tstrength\n"
            "DICER1\t10\tstrong\n"
            "CDKL5\t17\tmoderate\n",
            encoding="utf-8",
        )
        return SpliceExonOverrides(p)

    def test_override_downgrades_flat_default(self, tmp_path):
        # DICER1 flat splice default is Very Strong; exon-10 skip → Strong.
        s, ev = evaluate_vcep_pvs1(
            _pc("DICER1", ConsequenceType.SPLICE_DONOR, intron="10/26"),
            self._ov(tmp_path),
        )
        assert s == S.STRONG
        assert "exon-aware" in ev

    def test_acceptor_maps_to_next_exon(self, tmp_path):
        # Acceptor in intron 16 skips exon 17 → CDKL5 Moderate.
        s, _ = evaluate_vcep_pvs1(
            _pc("CDKL5", ConsequenceType.SPLICE_ACCEPTOR, intron="16/20"),
            self._ov(tmp_path),
        )
        assert s == S.MODERATE

    def test_no_matching_exon_keeps_flat_default(self, tmp_path):
        s, ev = evaluate_vcep_pvs1(
            _pc("DICER1", ConsequenceType.SPLICE_DONOR, intron="3/26"),
            self._ov(tmp_path),
        )
        assert s == S.VERY_STRONG
        assert "exon-aware" not in ev

    def test_no_overrides_keeps_flat_default(self):
        s, _ = evaluate_vcep_pvs1(
            _pc("DICER1", ConsequenceType.SPLICE_DONOR, intron="10/26"), None
        )
        assert s == S.VERY_STRONG
