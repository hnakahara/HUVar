"""APC-specific PVS1 decision tree (InSiGHT / Tayoun): codon-range gate +
allele-specific Lists A-E."""
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pvs1 import PVS1Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import (
    Assembly, ACMGCriterion, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.pvs1.apc import evaluate_apc_pvs1


def _pc(consequence, hgvs_c=None, protein_position=None, gene="APC"):
    return ConsequenceInfo(
        transcript_id="NM_000038.6", gene_id="ENSG", gene_symbol=gene,
        consequence=consequence, biotype="protein_coding",
        hgvs_c=hgvs_c, protein_position=protein_position,
    )


class TestApcListLookup:
    def test_list_a_pvs1(self):
        s, _ = evaluate_apc_pvs1(_pc(ConsequenceType.SPLICE_ACCEPTOR, "c.136-1G>A"))
        assert s == CriterionStrength.VERY_STRONG

    def test_allele_specific_downgrade(self):
        # c.835-1G>A is List A (PVS1); G>C/T is List C (Moderate).
        assert evaluate_apc_pvs1(_pc(ConsequenceType.SPLICE_ACCEPTOR, "c.835-1G>A"))[0] == CriterionStrength.VERY_STRONG
        assert evaluate_apc_pvs1(_pc(ConsequenceType.SPLICE_ACCEPTOR, "c.835-1G>C"))[0] == CriterionStrength.MODERATE

    def test_last_nt_exonic_strong(self):
        # "G to non-G last nucleotide" exonic change VEP calls missense → List B.
        s, _ = evaluate_apc_pvs1(_pc(ConsequenceType.MISSENSE, "c.422G>C"))
        assert s == CriterionStrength.STRONG

    def test_list_e_na(self):
        s, _ = evaluate_apc_pvs1(_pc(ConsequenceType.SPLICE_ACCEPTOR, "c.1959-1G>C"))
        assert s == CriterionStrength.NOT_MET

    def test_transcript_prefix_stripped(self):
        s, _ = evaluate_apc_pvs1(_pc(ConsequenceType.MISSENSE, "NM_000038.6:c.422G>T"))
        assert s == CriterionStrength.STRONG


class TestApcCodonGate:
    def test_in_range_pvs1(self):
        s, _ = evaluate_apc_pvs1(_pc(ConsequenceType.STOP_GAINED, "c.1500C>T", protein_position=500))
        assert s == CriterionStrength.VERY_STRONG

    def test_upstream_of_49_na(self):
        s, ev = evaluate_apc_pvs1(_pc(ConsequenceType.FRAMESHIFT, "c.100del", protein_position=34))
        assert s == CriterionStrength.NOT_MET
        assert "outside 49-2645" in ev

    def test_downstream_of_2645_na(self):
        s, _ = evaluate_apc_pvs1(_pc(ConsequenceType.STOP_GAINED, "c.8000C>T", protein_position=2700))
        assert s == CriterionStrength.NOT_MET

    def test_boundary_codons_inclusive(self):
        assert evaluate_apc_pvs1(_pc(ConsequenceType.STOP_GAINED, "c.x", protein_position=49))[0] == CriterionStrength.VERY_STRONG
        assert evaluate_apc_pvs1(_pc(ConsequenceType.STOP_GAINED, "c.x", protein_position=2645))[0] == CriterionStrength.VERY_STRONG

    def test_unknown_codon_defers(self):
        # No protein_position → defer to the generic tree (None).
        assert evaluate_apc_pvs1(_pc(ConsequenceType.FRAMESHIFT, "c.x", protein_position=None)) is None

    def test_non_truncating_non_listed_defers(self):
        # A plain missense not in any list → APC rules don't resolve it.
        assert evaluate_apc_pvs1(_pc(ConsequenceType.MISSENSE, "c.9999A>G")) is None


def _ev():
    cfg = MagicMock()
    cfg.disease_prevalence_tsv.exists.return_value = False  # no per-gene TSV
    return PVS1Evaluator(cfg)


def _ann(pc):
    return AnnotationData(consequences=[pc], gnomad=GnomADData(loeuf=0.10))


def _snv():
    return VariantRecord(chrom="chr5", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)


class TestApcEvaluatorRouting:
    def test_exonic_last_nt_fires_despite_missense(self):
        # APC c.422G>C is a missense by VEP but List B → PVS1_Strong. The generic
        # gate would have returned "Not a LoF consequence"; the APC path overrides.
        r = _ev().evaluate(_snv(), _ann(_pc(ConsequenceType.MISSENSE, "c.422G>C")))
        assert r.triggered
        assert r.strength == CriterionStrength.STRONG
        assert r.criterion == ACMGCriterion.PVS1

    def test_list_e_splice_withheld(self):
        r = _ev().evaluate(_snv(), _ann(_pc(ConsequenceType.SPLICE_ACCEPTOR, "c.1959-1G>T")))
        assert not r.triggered
        assert "N/A" in r.evidence

    def test_out_of_range_truncation_withheld(self):
        r = _ev().evaluate(_snv(), _ann(_pc(ConsequenceType.STOP_GAINED, "c.8000C>T", protein_position=2700)))
        assert not r.triggered

    def test_in_range_truncation_fires(self):
        r = _ev().evaluate(_snv(), _ann(_pc(ConsequenceType.STOP_GAINED, "c.1500C>T", protein_position=500)))
        assert r.triggered
        assert r.strength == CriterionStrength.VERY_STRONG

    def test_non_apc_gene_unaffected(self):
        # A BRCA1 missense still returns "Not a LoF consequence" (no APC path).
        r = _ev().evaluate(_snv(), _ann(_pc(ConsequenceType.MISSENSE, "c.422G>C", gene="BRCA1")))
        assert not r.triggered
        assert "Not a LoF" in r.evidence
