"""PM4 size-based Supporting downgrade (``pm4_supporting_max_aa`` column).

Many VCEPs apply PM4 at Supporting (not the default Moderate) for a small
in-frame indel — "single amino acid" (<=1 aa) or "< 3 amino acid residues"
(<=2 aa). The ``pm4_supporting_max_aa`` column encodes the size cutoff; the
evaluator measures the indel's net codon-length change and downgrades when it is
at or below the cutoff. Stop-loss is not size-scoped and keeps Moderate.
"""
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.pathogenic.pm4 import PM4Evaluator, _load_pm4_columns
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)

# GATM: single-aa -> Supporting (max 1). MECP2: <3 aa -> Supporting (max 2).
# MYH7: no size downgrade (default Moderate).
_TSV = (
    "gene_symbol\tpm4\tpm4_supporting_max_aa\n"
    "GATM\tapplicable\t1\n"
    "MECP2\tapplicable\t2\n"
    "MYH7\tapplicable\t\n"
    "BRCA1\tnot_applicable\t\n"
)


def _cfg(tmp_path):
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_TSV, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = p
    return cfg


def _ann(gene, consequence=ConsequenceType.INFRAME_DELETION):
    return AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=consequence, biotype="protein_coding",
    )])


def _var(ref, alt):
    return VariantRecord(chrom="chr1", pos=100, ref=ref, alt=alt, assembly=Assembly.GRCH38)


# 3-bp (1 aa) and 6-bp (2 aa) in-frame deletions.
_DEL_1AA = _var("ACGT", "A")       # diff 3 -> 1 aa
_DEL_2AA = _var("ACGTACG", "A")    # diff 6 -> 2 aa


class TestLoader:
    def test_reads_max_aa_and_declined(self, tmp_path):
        declined, max_aa = _load_pm4_columns(_cfg(tmp_path).disease_prevalence_tsv)
        assert "BRCA1" in declined
        assert max_aa == {"GATM": 1, "MECP2": 2}


class TestSizeDowngrade:
    def test_single_aa_downgrades_to_supporting(self, tmp_path):
        r = PM4Evaluator(_cfg(tmp_path)).evaluate(_DEL_1AA, _ann("GATM"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_two_aa_above_cutoff_stays_moderate(self, tmp_path):
        # GATM cutoff is 1 aa; a 2-aa indel keeps the default Moderate.
        r = PM4Evaluator(_cfg(tmp_path)).evaluate(_DEL_2AA, _ann("GATM"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_two_aa_within_cutoff_supporting(self, tmp_path):
        # MECP2 cutoff is 2 aa; a 2-aa indel downgrades to Supporting.
        r = PM4Evaluator(_cfg(tmp_path)).evaluate(_DEL_2AA, _ann("MECP2"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_stop_loss_not_size_scoped(self, tmp_path):
        # Stop-loss keeps Moderate even for a size-downgrade gene.
        r = PM4Evaluator(_cfg(tmp_path)).evaluate(
            _var("A", "T"), _ann("GATM", ConsequenceType.STOP_LOST)
        )
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_gene_without_cutoff_stays_moderate(self, tmp_path):
        r = PM4Evaluator(_cfg(tmp_path)).evaluate(_DEL_1AA, _ann("MYH7"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE


class TestParser:
    def _pm4(self, *strengths):
        return {"criteriaCodes": [{"label": "PM4", "evidenceStrengths": [
            {"label": lbl, "applicability": "Applicable", "description": desc}
            for lbl, desc in strengths
        ]}]}

    def test_single_amino_acid_is_1(self):
        rs = self._pm4(("Supporting", "Downgrade to PM4_Supporting for an in-frame "
                                      "deletion/insertion of a single amino acid."))
        assert bdt._pm4_supporting_max_aa(rs) == "1"

    def test_one_amino_acid_is_1(self):
        rs = self._pm4(("Supporting", "In frame deletion/insertions of one amino acid."))
        assert bdt._pm4_supporting_max_aa(rs) == "1"

    def test_less_than_three_is_2(self):
        rs = self._pm4(("Supporting", "Smaller in-frame events (< 3 amino acid residues)."))
        assert bdt._pm4_supporting_max_aa(rs) == "2"

    def test_no_size_supporting_is_blank(self):
        rs = self._pm4(("Moderate", "Protein length change outside a repeat region."))
        assert bdt._pm4_supporting_max_aa(rs) == ""
