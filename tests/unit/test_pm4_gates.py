"""PM4 conservation gate, deletion-content gate, and mutual-exclusion pass.

- Conservation (RPE65/CTLA4/PIK3R1): an in-frame indel earns PM4 only at a
  conserved position (phyloP > cutoff); the gate is skipped when phyloP is
  unavailable.
- Deletion-content (SCID panel): an in-frame deletion earns PM4 only if its
  deleted span contains a ClinVar P/LP (Moderate) or VUS (Supporting) variant.
- Mutual exclusion: the registry suppresses PM4 when PVS1 (FBN1) or PVS1/PP3
  (KCNQ1/CTLA4/PIK3R1) also fired.
"""
import sqlite3
from unittest.mock import MagicMock

from acmg_classifier.criteria.pm4_regions import PM4Regions
from acmg_classifier.criteria.pathogenic.pm4 import PM4Evaluator
from acmg_classifier.criteria.registry import CriteriaRegistry
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    Assembly, ACMGCriterion, ConsequenceType, CriterionStrength,
)
from acmg_classifier.models.variant import VariantRecord

_REGIONS = (
    "gene_symbol\tstrength\tregions\tresidues\n"
    "RPE65\tconserved_phylop\t2.0\t\n"
    "CTLA4\tconserved_phylop\t2.0\t\n"
    "CTLA4\texcludes\tPVS1,PP3\t\n"
    "FBN1\texcludes\tPVS1\t\n"
    "DCLRE1C\tdeletion_content\tyes\t\n"
    "ABCA4\tnt_phylop\t7.367\t\n"
)
_DP = (
    "gene_symbol\tpm4\tpm4_supporting_max_aa\n"
    "RPE65\tapplicable\t1\n"
)


def _paths(tmp_path):
    reg = tmp_path / "pm4_regions.tsv"
    reg.write_text(_REGIONS, encoding="utf-8")
    dp = tmp_path / "dp.tsv"
    dp.write_text(_DP, encoding="utf-8")
    return reg, dp


def _phylop_stub(score):
    class _Stub:
        def is_available(self):
            return score is not None

        def value(self, chrom, pos):
            return score
    return _Stub()


def _cfg(tmp_path, clinvar=None):
    reg, dp = _paths(tmp_path)
    cfg = MagicMock()
    cfg.pm4_regions_tsv = reg
    cfg.disease_prevalence_tsv = dp
    cfg.clinvar_sqlite = clinvar or (tmp_path / "absent.sqlite")
    return cfg


def _ann(gene, consequence, pos=None):
    return AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=consequence, biotype="protein_coding", protein_position=pos,
    )])


def _del(ref="ACGTACG", alt="A", pos=100):  # 2 aa deletion by default
    return VariantRecord(chrom="chr1", pos=pos, ref=ref, alt=alt, assembly=Assembly.GRCH38)


class TestConservationGate:
    def test_conserved_fires(self, tmp_path):
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(3.0)   # > 2.0 → conserved
        r = ev.evaluate(_del(), _ann("RPE65", ConsequenceType.INFRAME_DELETION, 200))
        assert r.triggered and r.strength == CriterionStrength.MODERATE  # 2 aa

    def test_conserved_single_aa_supporting(self, tmp_path):
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(3.0)
        r = ev.evaluate(_del("ACGT", "A"), _ann("RPE65", ConsequenceType.INFRAME_DELETION, 200))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING  # 1 aa

    def test_not_conserved_withheld(self, tmp_path):
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(0.5)   # <= 2.0 → not conserved
        r = ev.evaluate(_del(), _ann("RPE65", ConsequenceType.INFRAME_DELETION, 200))
        assert not r.triggered and "conserved" in r.evidence

    def test_phylop_unavailable_skips_gate(self, tmp_path):
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(None)  # unavailable → gate skipped → fires
        r = ev.evaluate(_del(), _ann("RPE65", ConsequenceType.INFRAME_DELETION, 200))
        assert r.triggered and r.strength == CriterionStrength.MODERATE


# --------------------------- deletion-content --------------------------------

_SCHEMA = """
CREATE TABLE variants (
    variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER
);
"""


def _clinvar(tmp_path, rows):
    p = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, chrom, pos, gene_symbol, "
        "clinical_significance, star_rating) VALUES (?,?,?,?,?,?)",
        rows,
    )
    con.commit()
    con.close()
    return p


class TestDeletionContent:
    # The deletion spans genomic [100, 106] (ref length 7).
    def test_plp_in_region_moderate(self, tmp_path):
        db = _clinvar(tmp_path, [("1", "chr1", 103, "DCLRE1C", "Pathogenic", 2)])
        ev = PM4Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_del(), _ann("DCLRE1C", ConsequenceType.INFRAME_DELETION, 50))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_vus_in_region_supporting(self, tmp_path):
        db = _clinvar(tmp_path, [("1", "chr1", 103, "DCLRE1C", "Uncertain significance", 2)])
        ev = PM4Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_del(), _ann("DCLRE1C", ConsequenceType.INFRAME_DELETION, 50))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_empty_region_not_met(self, tmp_path):
        db = _clinvar(tmp_path, [("1", "chr1", 500, "DCLRE1C", "Pathogenic", 2)])  # outside span
        ev = PM4Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_del(), _ann("DCLRE1C", ConsequenceType.INFRAME_DELETION, 50))
        assert not r.triggered and "no ClinVar" in r.evidence

    def test_insertion_not_content_gated(self, tmp_path):
        # The content rule is delete-only; an insertion uses the flat default.
        db = _clinvar(tmp_path, [])
        ev = PM4Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_del("A", "AGCT"), _ann("DCLRE1C", ConsequenceType.INFRAME_INSERTION, 50))
        assert r.triggered and r.strength == CriterionStrength.MODERATE


# --------------------------- mutual exclusion --------------------------------

class TestPM4Exclusions:
    def _registry(self, tmp_path):
        reg_tsv, _ = _paths(tmp_path)
        reg = object.__new__(CriteriaRegistry)
        reg._pm4_regions = PM4Regions(reg_tsv)
        return reg

    def test_fbn1_pm4_suppressed_by_pvs1(self, tmp_path):
        reg = self._registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PM4),
            CriteriaResult.met(ACMGCriterion.PVS1),
        ]
        reg._apply_pm4_exclusions(results, _ann("FBN1", ConsequenceType.INFRAME_DELETION))
        pm4 = next(r for r in results if r.criterion == ACMGCriterion.PM4)
        assert pm4.suppressed and "not with PVS1" in pm4.evidence

    def test_ctla4_pm4_suppressed_by_pp3(self, tmp_path):
        reg = self._registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PM4),
            CriteriaResult.met(ACMGCriterion.PP3),
        ]
        reg._apply_pm4_exclusions(results, _ann("CTLA4", ConsequenceType.INFRAME_DELETION))
        assert results[0].suppressed

    def test_pm4_kept_without_clash(self, tmp_path):
        reg = self._registry(tmp_path)
        results = [CriteriaResult.met(ACMGCriterion.PM4)]
        reg._apply_pm4_exclusions(results, _ann("FBN1", ConsequenceType.INFRAME_DELETION))
        assert not results[0].suppressed

    def test_non_listed_gene_untouched(self, tmp_path):
        reg = self._registry(tmp_path)
        results = [
            CriteriaResult.met(ACMGCriterion.PM4),
            CriteriaResult.met(ACMGCriterion.PVS1),
        ]
        reg._apply_pm4_exclusions(results, _ann("RPE65", ConsequenceType.INFRAME_DELETION))
        assert not results[0].suppressed


def _snv(ref="A", alt="G", pos=100):
    return VariantRecord(chrom="chr1", pos=pos, ref=ref, alt=alt, assembly=Assembly.GRCH38)


class TestNtPhylop:
    """ABCA4: a synonymous/missense variant affecting >1 nucleotide at a highly
    conserved position (phyloP >= 7.367) earns PM4_Moderate. A single-nucleotide
    change (SNV) does NOT qualify; <7.367, non-syn/mis, or phyloP off → not met."""

    def test_single_nt_synonymous_not_met(self, tmp_path):
        # A single-nucleotide synonymous SNV does not qualify (>1 nt required).
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(9.25)
        r = ev.evaluate(_snv(), _ann("ABCA4", ConsequenceType.SYNONYMOUS))
        assert not r.triggered

    def test_single_nt_missense_not_met(self, tmp_path):
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(7.367)
        r = ev.evaluate(_snv(), _ann("ABCA4", ConsequenceType.MISSENSE))
        assert not r.triggered

    def test_multi_nt_moderate(self, tmp_path):
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(8.0)
        r = ev.evaluate(_snv("AC", "GT"), _ann("ABCA4", ConsequenceType.MISSENSE))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_multi_nt_below_cutoff_not_met(self, tmp_path):
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(7.0)
        r = ev.evaluate(_snv("AC", "GT"), _ann("ABCA4", ConsequenceType.SYNONYMOUS))
        assert not r.triggered

    def test_phylop_unavailable_not_met(self, tmp_path):
        ev = PM4Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(None)
        r = ev.evaluate(_snv("AC", "GT"), _ann("ABCA4", ConsequenceType.SYNONYMOUS))
        assert not r.triggered
