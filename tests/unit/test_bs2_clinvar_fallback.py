"""BS2 ClinVar expert-panel fallback for genes whose VCEP bars gnomAD data.

CDH1 / TP53 / SERPINC1 force gnomAD-based BS2 to ``not_applicable`` because their
BS2 needs an internal cohort gnomAD cannot supply. When a >=3-star ClinVar review
has already applied BS2 to the specific variant, the evaluator harvests that
expert-panel judgement instead of withholding BS2 entirely.
"""
import sqlite3
from pathlib import Path
from unittest.mock import MagicMock

from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
from acmg_classifier.local_db import clinvar_sqlite
from acmg_classifier.local_db.clinvar_sqlite import query_bs2_benign_evidence
from acmg_classifier.models.annotation import AnnotationData, GnomADData, ConsequenceInfo
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_TSV = (
    "gene_symbol\tbs2\tinheritance\tbs2_count\tbs2_female_only\tbs2_hom_only\n"
    "CDH1\tnot_applicable\tAD\t\t\t\n"      # gnomAD BS2 barred → ClinVar fallback
    "GENE0\tapplicable\tAD\t\t\t\n"          # ordinary gnomAD-based BS2
)

# Schema must match clinvar_builder._CREATE_TABLE (19 columns).
_COLS = (
    "variation_id", "chrom", "pos", "ref", "alt", "gene_symbol", "hgvs_c",
    "hgvs_p", "amino_acid_change", "codon_position", "clinical_significance",
    "review_status", "star_rating", "last_evaluated", "affected_cases",
    "functional_evidence", "segregation_evidence", "bs2_evidence", "bs2_strength",
)


def _build_db(tmp_path: Path, rows: list[dict]) -> Path:
    db = tmp_path / "clinvar.sqlite"
    con = sqlite3.connect(str(db))
    con.execute(
        "CREATE TABLE variants ("
        "variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, "
        "gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT, "
        "codon_position INTEGER, clinical_significance TEXT, review_status TEXT, "
        "star_rating INTEGER, last_evaluated TEXT, affected_cases INTEGER, "
        "functional_evidence INTEGER, segregation_evidence INTEGER, "
        "bs2_evidence INTEGER DEFAULT 0, bs2_strength TEXT)"
    )
    for r in rows:
        vals = [r.get(c) for c in _COLS]
        con.execute(
            f"INSERT INTO variants VALUES ({','.join('?' * len(_COLS))})", vals
        )
    con.commit()
    con.close()
    # The query layer caches one immutable connection per path; drop it so a
    # freshly-built per-test DB is reopened rather than served from a stale cache.
    clinvar_sqlite._CONN_CACHE.pop(str(db), None)
    return db


def _bs2_row(strength: str | None = "Strong", stars: int = 3, **over) -> dict:
    base = dict(
        variation_id="9999", chrom="16", pos=68812210, ref="T", alt="C",
        gene_symbol="CDH1", hgvs_c=None, hgvs_p="NP_004351.1:p.Trp409Arg",
        amino_acid_change="W409R", codon_position=409,
        clinical_significance="Likely benign",
        review_status="reviewed by expert panel", star_rating=stars,
        last_evaluated="2024-01-01", affected_cases=0, functional_evidence=0,
        segregation_evidence=0, bs2_evidence=1, bs2_strength=strength,
    )
    base.update(over)
    return base


def _cfg(tmp_path: Path, db: Path) -> MagicMock:
    cfg = MagicMock()
    tsv = tmp_path / "disease_prevalence.tsv"
    tsv.write_text(_TSV, encoding="utf-8")
    cfg.disease_prevalence_tsv = tsv
    cfg.clinvar_sqlite = db
    cfg.bs2_min_homalt = 2
    cfg.bs2_min_hemi = 2
    cfg.bs2_min_het = 3
    return cfg


def _variant() -> VariantRecord:
    return VariantRecord(
        chrom="chr16", pos=68812210, ref="T", alt="C", assembly=Assembly.GRCH38
    )


def _ann(gene: str, gnomad: GnomADData | None = None) -> AnnotationData:
    return AnnotationData(
        gnomad=gnomad,
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
        )],
    )


class TestQueryBs2BenignEvidence:
    def test_expert_panel_bs2_found(self, tmp_path):
        db = _build_db(tmp_path, [_bs2_row("Strong")])
        assert query_bs2_benign_evidence(db, "chr16", 68812210, "T", "C") == (
            True, "Strong")

    def test_supporting_strength_returned(self, tmp_path):
        db = _build_db(tmp_path, [_bs2_row("Supporting")])
        assert query_bs2_benign_evidence(db, "16", 68812210, "T", "C") == (
            True, "Supporting")

    def test_strongest_across_conditions(self, tmp_path):
        # Two RCVs (different conditions) for the same variant cite different
        # strengths — the strongest wins.
        db = _build_db(tmp_path, [
            _bs2_row("Supporting", variation_id="1"),
            _bs2_row("Moderate", variation_id="2"),
        ])
        assert query_bs2_benign_evidence(db, "16", 68812210, "T", "C") == (
            True, "Moderate")

    def test_single_submitter_not_returned(self, tmp_path):
        # star_rating < 3 is excluded by the query.
        db = _build_db(tmp_path, [_bs2_row("Strong", stars=1)])
        assert query_bs2_benign_evidence(db, "16", 68812210, "T", "C") == (
            False, None)

    def test_no_bs2_row(self, tmp_path):
        db = _build_db(tmp_path, [_bs2_row("Strong", bs2_evidence=0)])
        assert query_bs2_benign_evidence(db, "16", 68812210, "T", "C") == (
            False, None)


class TestEvaluatorFallback:
    def test_not_applicable_gene_met_via_clinvar(self, tmp_path):
        db = _build_db(tmp_path, [_bs2_row("Strong")])
        ev = BS2Evaluator(_cfg(tmp_path, db))
        # gnomAD present but VCEP bars population BS2; ClinVar 3★ rescues it.
        r = ev.evaluate(_variant(), _ann("CDH1", GnomADData(filter_pass=True, ac=2)))
        assert r.triggered
        assert r.strength == CriterionStrength.STRONG
        assert "expert-panel" in r.evidence

    def test_fallback_without_gnomad_record(self, tmp_path):
        # Internal-cohort variants are often absent from gnomAD; the fallback
        # must still fire (the old early-return on missing gnomAD blocked it).
        db = _build_db(tmp_path, [_bs2_row("Strong")])
        ev = BS2Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_variant(), _ann("CDH1", gnomad=None))
        assert r.triggered

    def test_supporting_strength_applied(self, tmp_path):
        db = _build_db(tmp_path, [_bs2_row("Supporting")])
        ev = BS2Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_variant(), _ann("CDH1", gnomad=None))
        assert r.triggered
        assert r.strength == CriterionStrength.SUPPORTING

    def test_not_applicable_gene_without_clinvar_not_met(self, tmp_path):
        db = _build_db(tmp_path, [])  # empty ClinVar
        ev = BS2Evaluator(_cfg(tmp_path, db))
        r = ev.evaluate(_variant(), _ann("CDH1", GnomADData(filter_pass=True, ac=99)))
        assert not r.triggered
        assert "no clinvar" in r.evidence.lower()
