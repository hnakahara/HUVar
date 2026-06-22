"""PS1 paralogue / analogous-residue route (ps1_paralog_group / _strength).

RASopathy VCEP (GN004) grants PS1 from the same amino-acid change at the
analogous (same-numbered) residue of a sibling gene — HRAS/NRAS/KRAS,
MAP2K1/MAP2K2, SOS1/SOS2 — at the full comparator-derived strength. HBA2 (GN173)
grants only PS1_Moderate from its paralogue HBA1. The paralogue route fires only
when the same-gene rule did not.
"""
import csv
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.ps1_genes import PS1Spec
from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_ROOT = Path(__file__).resolve().parents[2]

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
        "gene_symbol\tps1\tps1_splice\tps1_max\tps1_paralog_group\tps1_paralog_strength\n"
        "KRAS\tapplicable\t\t\tHRAS,NRAS\t\n"
        "HBA2\tapplicable\t\t\tHBA1\tModerate\n"
        "BRCA1\tapplicable\t\t\t\t\n",
        encoding="utf-8",
    )
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = clinvar
    return cfg


def _ann(gene, codon=12, change=("Gly", "Val")):
    return AnnotationData(
        gnomad=GnomADData(),
        consequences=[ConsequenceInfo(
            transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
            is_mane_select=True, protein_position=codon,
            hgvs_p=f"NP_0:p.{change[0]}{codon}{change[1]}",
            amino_acids="G/V",
        )],
    )


def _snv():
    return VariantRecord(chrom="chr12", pos=100, ref="C", alt="A", assembly=Assembly.GRCH38)


class TestLoader:
    def test_group_and_strength(self, tmp_path):
        s = PS1Spec(_cfg(tmp_path, None).disease_prevalence_tsv)
        assert s.paralog_group("KRAS") == ("HRAS", "NRAS")
        assert s.paralog_strength("KRAS") is None
        assert s.paralog_group("HBA2") == ("HBA1",)
        assert s.paralog_strength("HBA2") == CriterionStrength.MODERATE
        assert s.paralog_group("BRCA1") == ()


class TestParalogEvaluator:
    def test_rasopathy_paralog_strong(self, tmp_path):
        # No KRAS G12V in ClinVar, but HRAS G12V Pathogenic exists → PS1 Strong
        # via the analogous-residue route.
        db = _db(tmp_path, [(
            "1", "chr11", 534289, "C", "A", "HRAS", "NP_h:p.Gly12Val", "G12V", 12,
            "Pathogenic", "criteria provided, single submitter", 2,
        )])
        r = PS1Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("KRAS"))
        assert r.triggered and r.strength == CriterionStrength.STRONG
        assert "paralogue" in r.evidence

    def test_hba2_paralog_fixed_moderate(self, tmp_path):
        # HBA1 paralogue hit (even Pathogenic) → fixed Moderate for HBA2.
        db = _db(tmp_path, [(
            "1", "chr16", 176680, "C", "A", "HBA1", "NP_a:p.Gly12Val", "G12V", 12,
            "Pathogenic", "criteria provided, single submitter", 2,
        )])
        r = PS1Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("HBA2"))
        assert r.triggered and r.strength == CriterionStrength.MODERATE

    def test_same_gene_takes_precedence(self, tmp_path):
        # A same-gene KRAS hit fires the normal route; the paralog route is not used.
        db = _db(tmp_path, [
            ("self", "chr12", 100, "C", "A", "KRAS", "NP_k:p.Gly12Val", "G12V", 12,
             "Pathogenic", "criteria provided", 2),  # this is the variant itself → excluded
            ("sib", "chr12", 101, "C", "G", "KRAS", "NP_k:p.Gly12Val", "G12V", 12,
             "Pathogenic", "criteria provided", 2),  # same-gene, same codon (within 2bp)
        ])
        r = PS1Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("KRAS"))
        assert r.triggered and r.strength == CriterionStrength.STRONG
        assert "paralogue" not in r.evidence   # same-gene route, not paralog

    def test_no_hit_anywhere_not_met(self, tmp_path):
        db = _db(tmp_path, [(
            "1", "chr11", 534289, "C", "A", "HRAS", "NP_h:p.Gly13Asp", "G13D", 13,
            "Pathogenic", "criteria provided", 2,
        )])  # different residue → no match for G12V
        r = PS1Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("KRAS"))
        assert not r.triggered

    def test_non_paralog_gene_unaffected(self, tmp_path):
        db = _db(tmp_path, [(
            "1", "chr11", 534289, "C", "A", "HRAS", "NP_h:p.Gly12Val", "G12V", 12,
            "Pathogenic", "criteria provided", 2,
        )])
        # BRCA1 has no paralog group → the HRAS hit must not leak in.
        r = PS1Evaluator(_cfg(tmp_path, db)).evaluate(_snv(), _ann("BRCA1"))
        assert not r.triggered


def test_committed_tsv_paralog_groups():
    tsv = _ROOT / "resources" / "shared" / "disease_prevalence.tsv"
    with tsv.open(encoding="utf-8") as f:
        rows = {r["gene_symbol"]: r for r in csv.DictReader(f, delimiter="\t")}
    assert rows["HRAS"]["ps1_paralog_group"] == "NRAS,KRAS"
    assert rows["HBA2"]["ps1_paralog_group"] == "HBA1"
    assert rows["HBA2"]["ps1_paralog_strength"] == "Moderate"
    assert rows["SOS1"]["ps1_paralog_group"] == "SOS2"


# --------------------- SCN paralogue alignment route -------------------------

from acmg_classifier.criteria.ps1_paralog import PS1ParalogMap  # noqa: E402

_ROOT = Path(__file__).resolve().parents[2]
_MAP_TSV = _ROOT / "resources" / "shared" / "ps1_paralog_map.tsv"


class TestParalogMapLoader:
    def test_committed_map_analogs(self):
        m = PS1ParalogMap(_MAP_TSV)
        assert m.has_gene("SCN1A") and m.has_gene("SCN8A")
        assert not m.has_gene("BRCA1")
        # Alignment row "40  40  40  44": SCN1A40 ↔ SCN2A40/SCN3A40/SCN8A44.
        assert m.analogs("SCN1A", 40) == {"SCN2A": 40, "SCN3A": 40, "SCN8A": 44}
        assert m.analogs("SCN8A", 44)["SCN1A"] == 40

    def test_missing_file(self, tmp_path):
        m = PS1ParalogMap(tmp_path / "nope.tsv")
        assert not m.has_gene("SCN1A")
        assert m.analogs("SCN1A", 40) == {}


def _scn_cfg(tmp_path, clinvar, map_tsv):
    dp = tmp_path / "dp.tsv"
    dp.write_text(
        "gene_symbol\tps1\tps1_splice\tps1_max\tps1_paralog_group\tps1_paralog_strength\n"
        "SCN2A\tapplicable\t\t\t\t\n",
        encoding="utf-8",
    )
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = dp
    cfg.clinvar_sqlite = clinvar
    cfg.ps1_paralog_map_tsv = map_tsv
    return cfg


def _map(tmp_path):
    p = tmp_path / "map.tsv"
    p.write_text("SCN1A\tSCN2A\tSCN3A\tSCN8A\n40\t40\t40\t44\n", encoding="utf-8")
    return p


def _scn_ann(gene, codon, change="R/W"):
    return AnnotationData(
        gnomad=GnomADData(),
        consequences=[ConsequenceInfo(
            transcript_id="NM_1", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
            is_mane_select=True, protein_position=codon,
            hgvs_p=f"NP:p.Arg{codon}Trp", amino_acids=change,
        )],
    )


class TestScnParalogEvaluator:
    def test_two_paralog_pathogenic_strong(self, tmp_path):
        # SCN2A R40W: analogous SCN1A R40W (Path) + SCN8A R44W (Path) → 2 → Strong.
        db = _db(tmp_path, [
            ("1", "chr2", 1, "C", "A", "SCN1A", "NP:p.Arg40Trp", "R40W", 40,
             "Pathogenic", "criteria provided", 2),
            ("2", "chr12", 1, "C", "A", "SCN8A", "NP:p.Arg44Trp", "R44W", 44,
             "Pathogenic", "criteria provided", 2),
        ])
        r = PS1Evaluator(_scn_cfg(tmp_path, db, _map(tmp_path))).evaluate(
            _snv(), _scn_ann("SCN2A", 40))
        assert r.triggered and r.strength == CriterionStrength.STRONG
        assert "analogous paralogue" in r.evidence

    def test_single_paralog_supporting(self, tmp_path):
        # Only one paralogue hit (LP) → Supporting.
        db = _db(tmp_path, [
            ("1", "chr2", 1, "C", "A", "SCN1A", "NP:p.Arg40Trp", "R40W", 40,
             "Likely pathogenic", "criteria provided", 2),
        ])
        r = PS1Evaluator(_scn_cfg(tmp_path, db, _map(tmp_path))).evaluate(
            _snv(), _scn_ann("SCN2A", 40))
        assert r.triggered and r.strength == CriterionStrength.SUPPORTING

    def test_no_paralog_hit_not_met(self, tmp_path):
        db = _db(tmp_path, [
            ("1", "chr2", 1, "C", "A", "SCN1A", "NP:p.Arg40Gln", "R40Q", 40,
             "Pathogenic", "criteria provided", 2),  # different alt → no match for R40W
        ])
        r = PS1Evaluator(_scn_cfg(tmp_path, db, _map(tmp_path))).evaluate(
            _snv(), _scn_ann("SCN2A", 40))
        assert not r.triggered

    def test_same_gene_takes_precedence(self, tmp_path):
        # A same-gene SCN2A hit (within codon window) fires first; paralog unused.
        db = _db(tmp_path, [
            ("self", "chr12", 101, "C", "G", "SCN2A", "NP:p.Arg40Trp", "R40W", 40,
             "Pathogenic", "criteria provided", 2),
        ])
        ev = PS1Evaluator(_scn_cfg(tmp_path, db, _map(tmp_path)))
        # candidate at chr12:100 so the sibling at 101 is a different-nucleotide
        # same-codon same-gene hit.
        r = ev.evaluate(_snv(), _scn_ann("SCN2A", 40))
        assert r.triggered and "paralogue" not in r.evidence


class TestKCNQ1Paralog:
    """KCNQ1 (GN112): a Pathogenic KCNQ2 same-AA variant at the analogous residue
    → PS1_Moderate (fixed); LP-only does not qualify."""

    def _cfg(self, tmp_path, clinvar):
        dp = tmp_path / "dp.tsv"
        dp.write_text(
            "gene_symbol\tps1\tps1_splice\tps1_max\tps1_paralog_group\tps1_paralog_strength\n"
            "KCNQ1\tapplicable\t\t\t\tModerate\n",
            encoding="utf-8",
        )
        m = tmp_path / "map.tsv"
        m.write_text("SCN1A\tSCN2A\tSCN3A\tSCN8A\tKCNQ1\tKCNQ2\n"
                     "\t\t\t\t144\t114\n", encoding="utf-8")
        cfg = MagicMock()
        cfg.disease_prevalence_tsv = dp
        cfg.clinvar_sqlite = clinvar
        cfg.ps1_paralog_map_tsv = m
        return cfg

    def _ann(self, codon=144, change="T/A"):
        return AnnotationData(
            gnomad=GnomADData(),
            consequences=[ConsequenceInfo(
                transcript_id="NM_1", gene_id="ENSG", gene_symbol="KCNQ1",
                consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
                is_mane_select=True, protein_position=codon,
                hgvs_p=f"NP:p.Thr{codon}Ala", amino_acids=change,
            )],
        )

    def test_kcnq2_pathogenic_moderate(self, tmp_path):
        # KCNQ1 p.Thr144Ala → analogous KCNQ2 p.Thr114Ala Pathogenic → Moderate.
        db = _db(tmp_path, [
            ("1", "chr20", 1, "A", "G", "KCNQ2", "NP:p.Thr114Ala", "T114A", 114,
             "Pathogenic", "criteria provided", 2),
        ])
        r = PS1Evaluator(self._cfg(tmp_path, db)).evaluate(_snv(), self._ann())
        assert r.triggered and r.strength == CriterionStrength.MODERATE
        assert "analogous paralogue" in r.evidence

    def test_kcnq2_lp_only_not_met(self, tmp_path):
        # The KCNQ1 rule requires KCNQ2 *Pathogenic*; LP-only does not qualify.
        db = _db(tmp_path, [
            ("1", "chr20", 1, "A", "G", "KCNQ2", "NP:p.Thr114Ala", "T114A", 114,
             "Likely pathogenic", "criteria provided", 2),
        ])
        r = PS1Evaluator(self._cfg(tmp_path, db)).evaluate(_snv(), self._ann())
        assert not r.triggered

    def test_no_kcnq2_hit_not_met(self, tmp_path):
        db = _db(tmp_path, [])
        r = PS1Evaluator(self._cfg(tmp_path, db)).evaluate(_snv(), self._ann())
        assert not r.triggered

    def test_committed_map_has_kcnq(self):
        m = PS1ParalogMap(_MAP_TSV)
        assert m.has_gene("KCNQ1") and m.has_gene("KCNQ2")
        assert m.analogs("KCNQ1", 144) == {"KCNQ2": 114}
