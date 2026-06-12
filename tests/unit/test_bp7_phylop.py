"""Per-gene BP7 phyloP "highly conserved" cutoff from ClinGen VCEP specs.

Most VCEPs define "not highly conserved" (BP7-eligible) as phyloP below a cutoff
that differs from the global default (phyloP100way 2.0): the neurodevelopmental
and coagulation panels use 0.1, VHL 0.2, the platelet GP genes 1.5, and RPGR 0
(only accelerated positions qualify). These tests cover the cspec extraction, the
TSV loader, the committed TSV, and the evaluator's per-gene gate.
"""
import importlib.util
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.benign.bp7 import BP7Evaluator
from acmg_classifier.criteria.bp_genes import BPApplicability
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, SpliceScore,
)
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

# --- load the build script as a module (it lives under scripts/, not the pkg) -
_BDT_PATH = Path(__file__).resolve().parents[2] / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT_PATH)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


def _bp7_code(desc, applicability="Applicable"):
    return {"label": "BP7", "evidenceStrengths": [
        {"label": "Supporting", "applicability": applicability, "description": desc},
    ]}


def _rs(*codes):
    return {"criteriaCodes": list(codes)}


class TestBp7PhylopExtraction:
    def test_eligible_side_lt(self):
        # "not highly conserved = PhyloP < 0.1" → cutoff 0.1.
        rs = _rs(_bp7_code("Defined 'not highly conserved' as PhyloP score <0.1 "
                           "and/or PhastCons <1."))
        assert bdt._bp7_phylop(rs) == "0.1"

    def test_conserved_side_gt(self):
        # "conservation = PhyloP > 0.1" → same boundary 0.1.
        rs = _rs(_bp7_code("Evolutionary conservation is defined as a PhyloP > 0.1."))
        assert bdt._bp7_phylop(rs) == "0.1"

    def test_le_value(self):
        rs = _rs(_bp7_code("PhyloP score is <=0.2 for conservation."))
        assert bdt._bp7_phylop(rs) == "0.2"

    def test_score_or_less(self):
        # GP1BA wording: "reference PhyloP (score = 1.5 or less)".
        rs = _rs(_bp7_code("SpliceAI (score = 0.2 or less) and reference PhyloP "
                           "(score = 1.5 or less) to assess conservation."))
        assert bdt._bp7_phylop(rs) == "1.5"

    def test_zero_cutoff(self):
        # RPGR: only accelerated (phyloP < 0) positions are BP7-eligible.
        rs = _rs(_bp7_code("SpliceAI <= 0.2 AND PhyloP < 0 for conservation."))
        assert bdt._bp7_phylop(rs) == "0"

    def test_phylop100way_default(self):
        rs = _rs(_bp7_code("the nucleotide is not highly conserved (phyloP100 way < 2.0)."))
        assert bdt._bp7_phylop(rs) == "2.0"

    def test_no_phylop_number(self):
        rs = _rs(_bp7_code("the nucleotide is not highly conserved."))
        assert bdt._bp7_phylop(rs) == ""

    def test_not_applicable_strength_ignored(self):
        rs = _rs(_bp7_code("PhyloP < 0.1", applicability="Not applicable"))
        assert bdt._bp7_phylop(rs) == ""


class TestBp7ConservationNotApplicable:
    """VCEPs that declare conservation non-informative → sentinel 'na'."""

    def test_not_required(self):
        # SCID T-/B-cell genes (ADA, IL7R, ...).
        rs = _rs(_bp7_code("Given poor conservation among vertebrates, nucleotide "
                           "conservation is not required in order to apply BP7."))
        assert bdt._bp7_phylop(rs) == "na"

    def test_does_not_have_to_be_considered(self):
        rs = _rs(_bp7_code("Conservation does not have to be considered for this "
                           "code to apply."))
        assert bdt._bp7_phylop(rs) == "na"

    def test_no_requirement_to_assess(self):
        # TP53 wording.
        rs = _rs(_bp7_code("No requirement to assess for nucleotide conservation "
                           "for rule application as per Walker et al., 2023."))
        assert bdt._bp7_phylop(rs) == "na"

    def test_no_conservation_requirement(self):
        # GALT wording.
        rs = _rs(_bp7_code("Use 'as is', but no conservation requirement."))
        assert bdt._bp7_phylop(rs) == "na"

    def test_not_considered_informative(self):
        # The user's reported wording.
        rs = _rs(_bp7_code("BP4 and BP7 can be added unless variant is in an "
                           "excluded region. Evolutionary conservation is not "
                           "considered informative for application of this code."))
        assert bdt._bp7_phylop(rs) == "na"

    def test_numeric_cutoff_wins_over_na(self):
        # A spec stating a number is a numeric gate, not 'na'.
        rs = _rs(_bp7_code("PhyloP < 0.1; conservation is required."))
        assert bdt._bp7_phylop(rs) == "0.1"


_TSV = (
    "gene_symbol\tbp7_phylop\n"
    "FOXG1\t0.1\n"
    "VHL\t0.2\n"
    "RPGR\t0\n"
    "GP1BA\t1.5\n"
    "TP53\tna\n"
    "NOSPEC\t\n"
)


def _tsv(tmp_path: Path) -> Path:
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_TSV, encoding="utf-8")
    return p


class TestBp7PhylopLoader:
    def test_reads_values(self, tmp_path):
        spec = BPApplicability(_tsv(tmp_path))
        assert spec.bp7_phylop("FOXG1") == 0.1
        assert spec.bp7_phylop("VHL") == 0.2
        assert spec.bp7_phylop("GP1BA") == 1.5

    def test_zero_is_a_cutoff_not_missing(self, tmp_path):
        # "0" must load as 0.0 (RPGR), distinct from "" → None.
        spec = BPApplicability(_tsv(tmp_path))
        assert spec.bp7_phylop("RPGR") == 0.0

    def test_blank_is_none(self, tmp_path):
        spec = BPApplicability(_tsv(tmp_path))
        assert spec.bp7_phylop("NOSPEC") is None

    def test_unknown_gene_none(self, tmp_path):
        spec = BPApplicability(_tsv(tmp_path))
        assert spec.bp7_phylop("UNKNOWN") is None

    def test_na_sentinel(self, tmp_path):
        # 'na' is not a numeric cutoff; it flags the conservation gate as off.
        spec = BPApplicability(_tsv(tmp_path))
        assert spec.bp7_phylop("TP53") is None
        assert spec.bp7_conservation_na("TP53") is True
        assert spec.bp7_conservation_na("FOXG1") is False
        assert spec.bp7_conservation_na("UNKNOWN") is False


class TestCommittedTsv:
    """Spot-check the committed TSV reflects the cspec extraction."""

    def _rows(self):
        import csv
        tsv = Path(__file__).resolve().parents[2] / "resources" / "clingen" / "disease_prevalence.tsv"
        with tsv.open(encoding="utf-8") as f:
            return {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}

    def test_tightened_genes(self):
        rows = self._rows()
        for gene in ("FOXG1", "MECP2", "CDKL5", "SLC9A6", "TCF4", "UBE3A",
                     "AKT3", "MTOR", "PIK3CA", "PIK3R2", "F8", "F9",
                     "SERPINC1", "RS1"):
            assert rows[gene]["bp7_phylop"] == "0.1", gene
        assert rows["VHL"]["bp7_phylop"] == "0.2"
        assert rows["RPGR"]["bp7_phylop"] == "0"
        for gene in ("GP1BA", "GP1BB", "GP9"):
            assert rows[gene]["bp7_phylop"] == "1.5", gene

    def test_phylop100way_default_genes(self):
        rows = self._rows()
        for gene in ("GCK", "HNF1A", "HNF4A", "KCNQ1"):
            assert rows[gene]["bp7_phylop"] == "2.0", gene

    def test_conservation_non_informative_genes(self):
        # Auto-detected from in-text wording.
        rows = self._rows()
        for gene in ("TP53", "ABCD1", "GALT", "ADA", "DCLRE1C", "IL7R",
                     "JAK3", "RAG1", "RAG2", "IL2RG"):
            assert rows[gene]["bp7_phylop"] == "na", gene

    def test_curated_lca_genes(self):
        # Leber Congenital Amaurosis VCEP (RPE65, GUCY2D, AIPL1): the
        # conservation-non-informative policy is in the published spec prose but
        # absent from the JSON-LD criteriaCode descriptions, so the build forces
        # 'na' via the curated _BP7_CONSERVATION_NA set.
        rows = self._rows()
        for gene in ("RPE65", "GUCY2D", "AIPL1"):
            assert rows[gene]["bp7_phylop"] == "na", gene


# --- evaluator per-gene gate -------------------------------------------------

def _phylop_stub(score):
    class _Stub:
        def is_available(self):
            return True

        def value(self, chrom, pos):
            return score
    return _Stub()


def _cfg(tmp_path: Path) -> MagicMock:
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = _tsv(tmp_path)
    cfg.bp7_phylop_max = 2.0
    cfg.phylop_bigwig = None
    return cfg


def _ann(gene: str) -> AnnotationData:
    return AnnotationData(
        splice=SpliceScore(tool="openspliceai", is_available=True, max_delta=0.02),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.SYNONYMOUS, biotype="protein_coding",
        )],
    )


def _snv() -> VariantRecord:
    return VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


class TestBp7PerGeneGate:
    def test_pergene_cutoff_blocks_below_default(self, tmp_path):
        # phyloP=0.5 passes the global 2.0 default but exceeds FOXG1's 0.1
        # cutoff → BP7 blocked for FOXG1.
        ev = BP7Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(0.5)
        r = ev.evaluate(_snv(), _ann("FOXG1"))
        assert not r.triggered
        assert "conserved" in r.evidence.lower()

    def test_default_gene_fires_at_same_score(self, tmp_path):
        # A gene with no VCEP cutoff uses the 2.0 default; phyloP=0.5 < 2.0 → fires.
        ev = BP7Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(0.5)
        r = ev.evaluate(_snv(), _ann("NOSPEC"))
        assert r.triggered

    def test_pergene_cutoff_fires_below_cutoff(self, tmp_path):
        # phyloP=0.05 is below FOXG1's 0.1 cutoff → not highly conserved → fires.
        ev = BP7Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(0.05)
        r = ev.evaluate(_snv(), _ann("FOXG1"))
        assert r.triggered

    def test_rpgr_zero_cutoff_blocks_neutral(self, tmp_path):
        # RPGR cutoff 0: a neutral position (phyloP=0.0 >= 0) is "conserved" →
        # blocked; only accelerated (negative) positions are BP7-eligible.
        ev = BP7Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(0.0)
        r = ev.evaluate(_snv(), _ann("RPGR"))
        assert not r.triggered

    def test_na_gene_skips_conservation_gate(self, tmp_path):
        # TP53 ('na'): conservation non-informative → even a highly conserved
        # position (phyloP=7.5) does NOT block BP7; it fires on splice alone.
        ev = BP7Evaluator(_cfg(tmp_path))
        ev._phylop = _phylop_stub(7.5)
        r = ev.evaluate(_snv(), _ann("TP53"))
        assert r.triggered
