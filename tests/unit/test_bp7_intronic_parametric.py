"""BP7 parametric intronic range (``bp7_intronic`` = "donor:N,acceptor:-M").

The Cardiomyopathy panel (MYH7, MYBPC3, TNNI3, TNNT2, TPM1, ACTC1, MYL2, MYL3)
and SLC6A8 extend BP7 to intronic variants "−4 and +7 outward" — more permissive
than the Walker −21 acceptor default but stricter than the noncanonical |dist|>=3
mode. These tests cover the cutoff parser, the build-script extraction (both
phrasings), the loader, and the evaluator's parametric eligibility.
"""
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.benign.bp7 import BP7Evaluator
from acmg_classifier.criteria.benign import bp7 as bp7_mod
from acmg_classifier.criteria.bp_genes import BPApplicability, _parse_intronic_cutoffs
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, SpliceScore
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


class TestParseCutoffs:
    def test_both_sides(self):
        assert _parse_intronic_cutoffs("donor:7,acceptor:-4") == (7, -4)

    def test_order_independent(self):
        assert _parse_intronic_cutoffs("acceptor:-4,donor:7") == (7, -4)

    def test_donor_only_defaults_acceptor(self):
        assert _parse_intronic_cutoffs("donor:9") == (9, -21)

    def test_garbage_is_none(self):
        assert _parse_intronic_cutoffs("noncanonical") is None


class TestBuildExtraction:
    def _bp7(self, desc):
        return {"criteriaCodes": [{"label": "BP7", "evidenceStrengths": [
            {"label": "Supporting", "applicability": "Applicable", "description": desc},
        ]}]}

    def test_cardiomyopathy_phrasing(self):
        rs = self._bp7("Also applicable to intronic variants outside the splice "
                       "consensus sequence (-4 and +7 outward) for which ...")
        assert bdt._bp7_intronic(rs) == "donor:7,acceptor:-4"

    def test_slc6a8_phrasing(self):
        rs = self._bp7("May also be applied ... for an intronic variant outside "
                       "the splice region (beyond -4bp or +7 bp)")
        assert bdt._bp7_intronic(rs) == "donor:7,acceptor:-4"

    def test_noncanonical_unaffected(self):
        rs = self._bp7("applies to any intronic variant except the canonical splice "
                       "dinucleotides")
        assert bdt._bp7_intronic(rs) == "noncanonical"

    def test_plain_is_blank(self):
        rs = self._bp7("A synonymous variant with no predicted splice impact.")
        assert bdt._bp7_intronic(rs) == ""


# ------------------------------ loader & evaluator ---------------------------

_TSV = (
    "gene_symbol\tbp7_phylop\tbp7_intronic\n"
    "MYH7\t\tdonor:7,acceptor:-4\n"
    "NOSPEC\t\t\n"
)


def _tsv(tmp_path):
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_TSV, encoding="utf-8")
    return p


def _cfg(tmp_path):
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _tsv(tmp_path)
    cfg.bp7_phylop_max = 2.0
    cfg.phylop_bigwig = None   # conservation gate skipped → isolate distance routing
    return cfg


def _snv():
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


def _ann(gene, consequence, dist):
    return AnnotationData(
        splice=SpliceScore(tool="openspliceai", is_available=True, max_delta=0.02),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=consequence, biotype="protein_coding",
            intron_distance_from_splice=dist,
        )],
    )


class TestLoader:
    def test_mode_and_cutoffs(self, tmp_path):
        s = BPApplicability(_tsv(tmp_path))
        assert s.bp7_intronic_mode("MYH7") == "parametric"
        assert s.bp7_intronic_cutoffs("MYH7") == (7, -4)

    def test_default_gene(self, tmp_path):
        s = BPApplicability(_tsv(tmp_path))
        assert s.bp7_intronic_mode("NOSPEC") == ""
        assert s.bp7_intronic_cutoffs("NOSPEC") is None


class TestParametricEvaluator:
    def _ev(self, tmp_path):
        return BP7Evaluator(_cfg(tmp_path))

    def test_acceptor_at_cutoff_fires(self, tmp_path):
        # -4 meets acceptor<=-4 (default mode would require -21).
        r = self._ev(tmp_path).evaluate(_snv(), _ann("MYH7", ConsequenceType.INTRON, -4))
        assert r.triggered

    def test_acceptor_inside_cutoff_blocked(self, tmp_path):
        # -3 is within the consensus for this gene (not <=-4, not >=+7).
        r = self._ev(tmp_path).evaluate(_snv(), _ann("MYH7", ConsequenceType.INTRON, -3))
        assert not r.triggered

    def test_acceptor_splice_region_routed(self, tmp_path):
        # -4 is annotated SPLICE_REGION; parametric mode pulls it in.
        r = self._ev(tmp_path).evaluate(_snv(), _ann("MYH7", ConsequenceType.SPLICE_REGION, -4))
        assert r.triggered

    def test_donor_at_cutoff_fires(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_snv(), _ann("MYH7", ConsequenceType.SPLICE_REGION, 7))
        assert r.triggered

    def test_donor_inside_cutoff_blocked(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_snv(), _ann("MYH7", ConsequenceType.INTRON, 6))
        assert not r.triggered

    def test_deep_acceptor_still_fires(self, tmp_path):
        r = self._ev(tmp_path).evaluate(_snv(), _ann("MYH7", ConsequenceType.INTRON, -30))
        assert r.triggered

    def test_default_gene_unaffected_at_minus4(self, tmp_path):
        # NOSPEC keeps the Walker default: -4 is NOT deep enough (needs -21).
        r = self._ev(tmp_path).evaluate(_snv(), _ann("NOSPEC", ConsequenceType.INTRON, -4))
        assert not r.triggered

    def test_intronic_eligible_unit(self):
        class _PC:
            intron_distance_from_splice = -4
        assert bp7_mod._intronic_eligible(_PC(), "parametric", (7, -4)) is True
        _PC.intron_distance_from_splice = -3
        assert bp7_mod._intronic_eligible(_PC(), "parametric", (7, -4)) is False
        _PC.intron_distance_from_splice = 7
        assert bp7_mod._intronic_eligible(_PC(), "parametric", (7, -4)) is True
