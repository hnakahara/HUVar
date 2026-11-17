"""PS3/BS3 from the variant-level functional tables the VCEPs attach to their CSpec
specifications (resources/shared/vcep_functional_codes.tsv): TP53 Table S3, OTC
yeast complementation (Lo 2023, VCEP thresholds), GALT activity table, PAH
(current eRepo VCEP classifications, else the pilot assignments) and IDUA
Appendix 3."""
from pathlib import Path

import pytest

from acmg_classifier.config import Config
from acmg_classifier.criteria.benign.manual import ManualBenignEvaluator
from acmg_classifier.criteria.functional_assays import (
    SRC_VCEP_TABLE, FunctionalAssays, _seq_matches,
)
from acmg_classifier.criteria.pathogenic.ps3 import PS3Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, SpliceScore
from acmg_classifier.models.enums import ACMGCriterion as C, Assembly, ConsequenceType as CT
from acmg_classifier.models.enums import CriterionStrength as S
from acmg_classifier.models.variant import VariantRecord


def _cfg(**kw):
    return Config(data_dir=Path("./nonexistent-data"), **kw)


def _ann(gene, tx, hgvs_p=None, hgvs_c=None, codons=None, csq=CT.MISSENSE, splice=0.01):
    c = ConsequenceInfo(transcript_id=tx, gene_id="g", gene_symbol=gene, consequence=csq,
                        biotype="protein_coding", is_mane_select=True, hgvs_p=hgvs_p,
                        hgvs_c=hgvs_c, codons=codons)
    sp = None if splice is None else SpliceScore(tool="openspliceai", max_delta=splice,
                                                 is_available=True)
    return AnnotationData(consequences=[c], splice=sp)


_V = VariantRecord(chrom="chr17", pos=7675088, ref="C", alt="T", assembly=Assembly.GRCH38)


def _call(ann):
    return FunctionalAssays(_cfg()).call(_V, ann)


# --- TP53 Table S3 -----------------------------------------------------------------

def test_tp53_table_code():
    fa = _call(_ann("TP53", "NM_000546.6", "NP_000537.3:p.Arg175His"))
    assert fa.criterion == C.PS3 and fa.strength == S.STRONG
    assert fa.evidence.startswith(SRC_VCEP_TABLE) and "Table S3" in fa.evidence


def test_tp53_no_evidence_is_assessed_without_mavedb_fallback():
    # Ala138Gly: Kato non-functional but LOF not by the majority -> "No evidence".
    fa = _call(_ann("TP53", "NM_000546.6", "NP_000537.3:p.Ala138Gly"))
    assert fa.criterion is None and fa.evidence.startswith(SRC_VCEP_TABLE)


def test_tp53_matched_from_vep_amino_acids_without_hgvs_p():
    c = ConsequenceInfo(transcript_id="NM_000546.6", gene_id="g", gene_symbol="TP53",
                        consequence=CT.MISSENSE, biotype="protein_coding", is_mane_select=True,
                        amino_acids="R/H", protein_position=175)
    fa = _call(AnnotationData(consequences=[c]))
    assert fa.criterion == C.PS3


@pytest.mark.parametrize("codons, hgvs_c, expected", [
    ("caC/caA", None, (C.PS3, S.STRONG)),         # His178Gln C>A
    ("caC/caG", None, (None, None)),              # His178Gln C>G: no evidence
])
def test_tp53_codon_specific_codes(codons, hgvs_c, expected):
    fa = _call(_ann("TP53", "NM_000546.6", "NP_000537.3:p.His178Gln", hgvs_c, codons))
    assert (fa.criterion, fa.strength) == expected


def test_tp53_codon_specific_unresolved_without_codons():
    fa = _call(_ann("TP53", "NM_000546.6", "NP_000537.3:p.His178Gln"))
    assert fa.criterion is None and "not resolved" in fa.evidence


def test_tp53_c_dot_specific_row():
    fa = _call(_ann("TP53", "NM_000546.6", "NP_000537.3:p.Gly187Asp",
                    "NM_000546.6:c.560G>A", "gGt/gAt"))
    assert fa.criterion == C.BS3 and fa.strength == S.SUPPORTING


def test_tp53_codon_change_rows():
    fa = _call(_ann("TP53", "NM_000546.6", "NP_000537.3:p.Arg282Val", None, "CGG/GTC"))
    assert fa.criterion == C.BS3
    fa = _call(_ann("TP53", "NM_000546.6", "NP_000537.3:p.Arg282Val", None, "CGG/GTG"))
    assert fa.criterion == C.PS3 and fa.strength == S.SUPPORTING


@pytest.mark.parametrize("splice", [None, 0.5])
def test_bs3_withheld_without_low_splice_score(splice):
    fa = _call(_ann("TP53", "NM_000546.6", "NP_000537.3:p.Ala129Ser", splice=splice))
    assert fa.criterion is None and "BS3 withheld" in fa.evidence


def test_tp53_table_before_nci_flowchart():
    ann = _ann("TP53", "NM_000546.6", "NP_000537.3:p.Arg175His")
    r = PS3Evaluator(_cfg(use_tp53_functional=True)).evaluate(_V, ann)
    assert r.triggered and r.evidence.startswith(SRC_VCEP_TABLE)


def test_seq_matches():
    assert _seq_matches("CC>TT", "", "CCg/TTg")
    assert _seq_matches("CGG>CTA", "", "cgg/cta")
    assert not _seq_matches("G>T", "", "Gtg/Ctg")
    assert _seq_matches("c.560G>A", "c.560G>A", None)
    assert _seq_matches("G>T", "", None) is None


# --- OTC yeast complementation ------------------------------------------------------

@pytest.mark.parametrize("p, expected", [
    ("Arg141Gln", (C.PS3, S.MODERATE)),     # growth -0.05
    ("Thr150Ile", (C.BS3, S.MODERATE)),     # growth 0.94
    ("Leu166Phe", (C.BS3, S.SUPPORTING)),   # growth 0.72
    ("Gln270Arg", (None, None)),            # SMG motif
])
def test_otc_thresholds(p, expected):
    fa = _call(_ann("OTC", "NM_000531.6", f"NP_000522.3:p.{p}"))
    assert (fa.criterion, fa.strength) == expected
    assert "Lo et al. 2023" in fa.evidence


def test_otc_not_for_predicted_splice_impact():
    fa = _call(_ann("OTC", "NM_000531.6", "NP_000522.3:p.Arg141Gln", splice=0.5))
    assert fa.criterion is None and "splicing impact" in fa.evidence


def test_otc_bs3_through_manual_benign_evaluator():
    ann = _ann("OTC", "NM_000531.6", "NP_000522.3:p.Thr150Ile")
    bs3 = next(r for r in ManualBenignEvaluator(_cfg()).evaluate(_V, ann)
               if r.criterion == C.BS3)
    assert bs3.triggered and bs3.strength == S.MODERATE


# --- GALT / PAH ----------------------------------------------------------------------

def test_galt_codes():
    fa = _call(_ann("GALT", "NM_000155.4", "NP_000146.2:p.Gln188Arg"))
    assert fa.criterion == C.PS3 and fa.strength == S.SUPPORTING
    fa = _call(_ann("GALT", "NM_000155.4", "NP_000146.2:p.Asn314Asp"))
    assert fa.criterion == C.BS3 and fa.strength == S.SUPPORTING
    fa = _call(_ann("GALT", "NM_000155.4", "NP_000146.2:p.Thr350Ala"))
    assert fa.criterion is None and "discordant" in fa.evidence


def test_pah_current_vcep_classification_wins():
    # c.1222C>T (p.Arg408Trp): eRepo VCEP classification uses PS3 (Strong).
    fa = _call(_ann("PAH", "NM_000277.3", "NP_000268.1:p.Arg408Trp", "NM_000277.3:c.1222C>T"))
    assert fa.criterion == C.PS3 and fa.strength == S.STRONG and "eRepo" in fa.evidence
    # c.922C>T: pilot PS3_Supporting, current classification without PS3.
    fa = _call(_ann("PAH", "NM_000277.3", "NP_000268.1:p.Arg308Cys", "NM_000277.3:c.922C>T"))
    assert fa.criterion is None and "pilot code PS3_supporting" in fa.evidence
    # c.208T>C: pilot PS3, current PS3_Supporting.
    fa = _call(_ann("PAH", "NM_000277.3", "NP_000268.1:p.Ser70Pro", "NM_000277.3:c.208T>C"))
    assert fa.criterion == C.PS3 and fa.strength == S.SUPPORTING


def test_pah_pilot_only_ps3_stored_as_moderate():
    fa = _call(_ann("PAH", "NM_000277.3", "NP_000268.1:p.Leu245Gln", "NM_000277.3:c.734T>A"))
    assert fa.criterion == C.PS3 and fa.strength == S.MODERATE and "Kanavy" in fa.evidence


@pytest.mark.parametrize("c, p, expected", [
    ("c.623G>A", "Gly208Asp", (C.PS3, S.MODERATE)),     # Matte 0.35% + Western blot
    ("c.719A>G", "His240Arg", (C.PS3, S.SUPPORTING)),   # Yu 0.4-0.9% and Beesley ~1.5%
    ("c.236C>T", "Ala79Val", (C.PS3, S.SUPPORTING)),    # Yogalingam <1%
    ("c.965T>A", "Val322Glu", (C.BS3, S.SUPPORTING)),   # Yu ~10% (pseudodeficiency range)
    ("c.246C>G", "His82Gln", (C.BS3, S.SUPPORTING)),    # Yu 26%, Yogalingam 25%
    ("c.793G>C", "Gly265Arg", (None, None)),            # Yogalingam 5%: between thresholds
])
def test_idua_appendix3(c, p, expected):
    fa = _call(_ann("IDUA", "NM_000203.5", f"NP_000194.2:p.{p}", f"NM_000203.5:{c}"))
    assert (fa.criterion, fa.strength) == expected
    assert "IDUA VCEP Appendix 3" in fa.evidence


def test_null_variants_and_other_transcripts_are_skipped():
    assert _call(_ann("GALT", "NM_000155.4", "NP_000146.2:p.Arg204Ter",
                      csq=CT.STOP_GAINED)) is None
    assert FunctionalAssays(_cfg()).vcep_call(
        _V, _ann("OTC", "NM_999999.1", "NP_1:p.Arg141Gln")) is None
