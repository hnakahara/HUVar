"""SCID VCEP PM5 point table for nonsense variants (criteria/pm5_ptc.py, group C)."""
import sqlite3
from unittest.mock import MagicMock

import pytest

from acmg_classifier.criteria.pm5_ptc import SCID_EVIDENCE_PREFIX, evaluate_scid_pm5
from acmg_classifier.criteria.registry import _apply_scid_pm5
from acmg_classifier.models.annotation import ConsequenceInfo
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion as A, Assembly, ConsequenceType, CriterionStrength as S,
)
from acmg_classifier.models.variant import VariantRecord

# ADA NM_000022: CDS starts at cDNA 93; exon 4 = codons ~73-120 (NMD), exon 12
# (final) = codons 360-363.
_SCHEMA = """CREATE TABLE variants (variation_id TEXT, chrom TEXT, pos INTEGER, ref TEXT,
    alt TEXT, gene_symbol TEXT, hgvs_c TEXT, hgvs_p TEXT, amino_acid_change TEXT,
    codon_position INTEGER, clinical_significance TEXT, review_status TEXT,
    star_rating INTEGER)"""


def _db(tmp_path, rows):
    p = tmp_path / "cv.sqlite"
    con = sqlite3.connect(p)
    con.execute(_SCHEMA)
    con.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit()
    con.close()
    return p


def _row(vid, codon, sig, fs=None):
    hp = (f"NP_000013.2:p.Leu{codon}fsTer{fs}" if fs else f"NP_000013.2:p.Arg{codon}Ter")
    return (vid, "20", 1000 + int(vid), "C", "T", "ADA", f"NM_000022.4:c.{codon * 3 - 2}C>T",
            hp, None, codon, sig, "reviewed by expert panel", 3)


@pytest.fixture(autouse=True)
def _accept_all(monkeypatch):
    from acmg_classifier.criteria import vcep_comparators
    monkeypatch.setattr(vcep_comparators, "filter_comparators",
                        lambda gene, crit, hits, asm: (hits, ""))


def _run(tmp_path, rows, codon):
    cfg = MagicMock()
    cfg.clinvar_sqlite = _db(tmp_path, rows)
    cfg.pm5_min_stars = 1
    pc = ConsequenceInfo(transcript_id="NM_000022.4", gene_id="ENSG", gene_symbol="ADA",
                         consequence=ConsequenceType.STOP_GAINED, biotype="protein_coding",
                         protein_position=codon, hgvs_p=f"NP_000013.2:p.Gln{codon}Ter")
    v = VariantRecord(chrom="chr20", pos=1, ref="C", alt="T", assembly=Assembly.GRCH38)
    return evaluate_scid_pm5(pc, cfg, v)


def test_nmd_exon_points(tmp_path):
    st, ev = _run(tmp_path, [_row("1", 90, "Pathogenic"), _row("2", 110, "Likely pathogenic"),
                             _row("3", 200, "Pathogenic")], 100)       # codon 200: other exon
    assert st == S.MODERATE and "2 point(s)" in ev and ev.startswith(SCID_EVIDENCE_PREFIX)
    (tmp_path / "x").mkdir()
    st, ev = _run(tmp_path / "x", [_row("1", 90, "Pathogenic"), _row("4", 115, "Benign")], 100)
    assert st == S.NOT_MET and "-1 point(s)" in ev and "reconsider PVS1" in ev


def test_final_exon_points(tmp_path):
    # Final exon: P/LP with a PTC downstream +1; B/LB with a PTC upstream -2;
    # a P/LP upstream does not count.
    st, ev = _run(tmp_path, [_row("1", 363, "Pathogenic"), _row("2", 360, "Pathogenic")], 362)
    assert st == S.SUPPORTING and "1 point(s)" in ev and "final exon" in ev


def test_four_points_strong(tmp_path):
    rows = [_row(str(i), c, "Pathogenic") for i, c in enumerate((80, 85, 95, 105), 1)]
    assert _run(tmp_path, rows, 100)[0] == S.STRONG


def _res(pm5_strength, pvs1=None, pm2=True):
    out = [CriteriaResult.met(A.PM5, pm5_strength, f"{SCID_EVIDENCE_PREFIX} ADA point table): 4")]
    out.append(CriteriaResult.met(A.PM2, S.SUPPORTING, "rare") if pm2
               else CriteriaResult.not_met(A.PM2, "common"))
    out.append(CriteriaResult.met(A.PVS1, pvs1, "null") if pvs1
               else CriteriaResult.not_met(A.PVS1, "-"))
    return out


def test_registry_pm2_and_pvs1():
    r = _res(S.STRONG, pm2=False)
    _apply_scid_pm5(r)
    assert r[0].suppressed
    r = _res(S.STRONG, pvs1=S.STRONG)
    _apply_scid_pm5(r)
    assert r[0].strength == S.MODERATE and "downgraded" in r[0].evidence
    r = _res(S.MODERATE, pvs1=S.VERY_STRONG)
    _apply_scid_pm5(r)
    assert r[0].strength == S.SUPPORTING
    r = _res(S.MODERATE, pvs1=S.STRONG)
    _apply_scid_pm5(r)
    assert r[0].strength == S.MODERATE
