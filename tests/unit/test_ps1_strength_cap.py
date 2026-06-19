"""PS1 per-gene strength cap (``ps1_max`` column).

Some VCEPs downgrade PS1 below its Strong default — RMRP's spec states
"Downgraded to PS1_Supporting". The ``ps1_max`` column carries that ceiling and
the evaluator clamps the comparator-derived strength to it. These tests cover
the loader, the evaluator clamp, the build-script extraction, and the committed
TSV value for RMRP.
"""
import csv
import importlib.util
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.ps1_genes import PS1Spec
from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_ROOT = Path(__file__).resolve().parents[2]
_BDT = _ROOT / "scripts" / "build_disease_thresholds.py"
_spec = importlib.util.spec_from_file_location("build_disease_thresholds", _BDT)
bdt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdt)


# ------------------------------- loader --------------------------------------

class TestLoader:
    def _tsv(self, tmp_path):
        p = tmp_path / "dp.tsv"
        p.write_text(
            "gene_symbol\tps1\tps1_splice\tps1_max\n"
            "RMRP\tapplicable\t\tSupporting\n"
            "GENE_MOD\tapplicable\t\tModerate\n"
            "BRCA1\tapplicable\t\t\n",
            encoding="utf-8",
        )
        return p

    def test_reads_cap(self, tmp_path):
        s = PS1Spec(self._tsv(tmp_path))
        assert s.max_strength("RMRP") == CriterionStrength.SUPPORTING
        assert s.max_strength("GENE_MOD") == CriterionStrength.MODERATE
        assert s.max_strength("BRCA1") is None
        assert s.max_strength("UNSEEN") is None


# ------------------------------ evaluator ------------------------------------

_SCHEMA = """
CREATE TABLE variants (
    variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT,
    gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER
);
"""


def _db(tmp_path, rows):
    p = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany(
        "INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, "
        "hgvs_p, amino_acid_change, codon_position, clinical_significance, "
        "review_status, star_rating) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
        rows,
    )
    con.commit()
    con.close()
    return p


def _cfg(tmp_path, clinvar):
    tsv = tmp_path / "dp.tsv"
    tsv.write_text(
        "gene_symbol\tps1\tps1_splice\tps1_max\n"
        "GENECAP\tapplicable\t\tSupporting\n"
        "GENEFULL\tapplicable\t\t\n",
        encoding="utf-8",
    )
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = clinvar
    return cfg


def _ann(gene):
    return AnnotationData(
        gnomad=GnomADData(),
        consequences=[ConsequenceInfo(
            transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
            is_mane_select=True, protein_position=175, hgvs_p="NP_0:p.Gly175Arg",
            amino_acids="G/R",
        )],
    )


def _snv():
    # Different nucleotide than the comparator below (so it is not self-excluded).
    return VariantRecord(chrom="chr1", pos=100, ref="C", alt="T", assembly=Assembly.GRCH38)


class TestEvaluatorCap:
    def _pathogenic_sibling(self, tmp_path, gene):
        # Same AA change (G175R) via a DIFFERENT nucleotide, classified Pathogenic
        # → comparator strength would be Strong without a cap.
        return _db(tmp_path, [(
            "1", "chr1", 100, "C", "G", gene, "NP_0:p.Gly175Arg", "G175R", 175,
            "Pathogenic", "criteria provided, single submitter", 2,
        )])

    def test_capped_gene_clamped_to_supporting(self, tmp_path):
        db = self._pathogenic_sibling(tmp_path, "GENECAP")
        r = PS1Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("GENECAP"))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_uncapped_gene_stays_strong(self, tmp_path):
        db = self._pathogenic_sibling(tmp_path, "GENEFULL")
        r = PS1Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("GENEFULL"))
        assert r.triggered and r.strength == CriterionStrength.STRONG


# --------------------------- build-script parser -----------------------------

class TestPs1MaxExtraction:
    def _ps1(self, *strengths):
        return {"criteriaCodes": [{"label": "PS1", "evidenceStrengths": [
            {"label": lbl, "applicability": "Applicable", "description": desc}
            for lbl, desc in strengths
        ]}]}

    def test_supporting_only_caps(self):
        rs = self._ps1(("Supporting", "Downgraded to PS1_Supporting. Same nucleotide position rule."))
        assert bdt._ps1_max(rs) == "Supporting"

    def test_strong_applicable_no_cap(self):
        rs = self._ps1(("Strong", "Same amino acid as established pathogenic variant."))
        assert bdt._ps1_max(rs) == ""

    def test_moderate_top_caps_moderate(self):
        rs = self._ps1(("Moderate", "..."), ("Supporting", "..."))
        assert bdt._ps1_max(rs) == "Moderate"


# ----------------------------- committed TSV ---------------------------------

def test_committed_rmrp_capped_supporting():
    tsv = _ROOT / "resources" / "shared" / "disease_prevalence.tsv"
    with tsv.open(encoding="utf-8") as f:
        rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
    assert rows["RMRP"]["ps1_max"] == "Supporting"
