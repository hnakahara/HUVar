"""RUNX1 canonical +/-1,2 splice strengths from the Myeloid Malignancy VCEP
"Summary of splicing effects" table (per intron donor / acceptor)."""
import pytest

from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.enums import ConsequenceType as C, CriterionStrength as S
from acmg_classifier.pvs1.vcep_pvs1 import _canonical_splice_site, evaluate_vcep_pvs1


def _pc(hgvs_c, cons):
    return ConsequenceInfo(
        transcript_id="NM_001754.5", gene_id="ENSG", gene_symbol="RUNX1",
        consequence=cons, biotype="protein_coding", is_mane_select=True,
        hgvs_c=f"NM_001754.5:{hgvs_c}",
    )


@pytest.mark.parametrize("hgvs_c,cons,expected", [
    ("c.58+1G>A", C.SPLICE_DONOR, S.NOT_MET),
    ("c.97+1G>A", C.SPLICE_DONOR, S.NOT_MET),
    ("c.98-1G>A", C.SPLICE_ACCEPTOR, S.NOT_MET),
    ("c.351+1G>A", C.SPLICE_DONOR, S.VERY_STRONG),
    ("c.352-1G>A", C.SPLICE_ACCEPTOR, S.VERY_STRONG),
    ("c.508+2T>C", C.SPLICE_DONOR, S.VERY_STRONG),
    ("c.509-1G>A", C.SPLICE_ACCEPTOR, S.STRONG),
    ("c.613+1G>A", C.SPLICE_DONOR, S.STRONG),
    ("c.805+1G>T", C.SPLICE_DONOR, S.STRONG),
    ("c.806-1G>A", C.SPLICE_ACCEPTOR, S.STRONG),
    ("c.967+2_967+5del", C.SPLICE_DONOR, S.STRONG),
    ("c.968-2A>G", C.SPLICE_ACCEPTOR, S.STRONG),
])
def test_runx1_splice_table(hgvs_c, cons, expected):
    strength, ev = evaluate_vcep_pvs1(_pc(hgvs_c, cons))
    assert strength == expected
    assert "splicing-effects table" in ev


def test_exonic_side_deletion_uses_donor_site():
    # c.508_508+10del spans the exon end and the c.508 donor -> PVS1
    strength, _ = evaluate_vcep_pvs1(_pc("c.508_508+10del", C.SPLICE_DONOR))
    assert strength == S.VERY_STRONG


@pytest.mark.parametrize("h,site", [
    ("c.508+1G>A", ("donor", 508)), ("c.509-2A>G", ("acceptor", 509)),
    ("c.967+2_967+5del", ("donor", 967)), ("c.967+5G>A", None),
    ("c.98-12C>T", None), ("c.1234G>A", None),
])
def test_site_parser(h, site):
    assert _canonical_splice_site(h) == site
