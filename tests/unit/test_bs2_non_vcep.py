"""BS2 for genes without a VCEP specification (Genome Medicine revision,
Reviewer 3): no gnomAD-based BS2 by any route (homozygotes, hemizygotes or
heterozygous carriers); BS2 only from curated evidence or an expert-panel
(>=3-star) assertion in ClinVar."""
from pathlib import Path
from unittest.mock import MagicMock, patch

from acmg_classifier.criteria.benign.bs2 import BS2Evaluator
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import Assembly, ConsequenceType
from acmg_classifier.models.variant import VariantRecord

_DP = "gene_symbol\tbs2\tinheritance\tbs2_count\tbs2_female_only\tbs2_hom_only\nVCEPGENE\tapplicable\tAR\t\t\t\n"


def _cfg(tmp_path: Path) -> MagicMock:
    dp = tmp_path / "disease_prevalence.tsv"
    dp.write_text(_DP, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = dp
    cfg.bs2_min_homalt, cfg.bs2_min_hemi, cfg.bs2_min_het = 2, 2, 3
    cfg.exclude_self_expert_panel = False
    cfg.curated_evidence = True
    cfg.clinvar_sqlite = tmp_path / "absent.sqlite"
    return cfg


def _ann(gene: str, **gd) -> AnnotationData:
    return AnnotationData(
        gnomad=GnomADData(filter_pass=True, **gd),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding")],
    )


_V = VariantRecord(chrom="chr1", pos=100, ref="A", alt="G", assembly=Assembly.GRCH38)


def test_no_vcep_no_gnomad_route(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path))
    for gd in (dict(ac=50, nhomalt=5), dict(ac=50, nhemi=5), dict(ac=50, nhomalt=0, loeuf=0.1)):
        r = ev.evaluate(_V, _ann("NOVCEP", **gd))
        assert not r.triggered
        assert "no VCEP specification" in r.evidence and "curated evidence only" in r.evidence


def test_no_vcep_expert_panel_import(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path))
    with patch("acmg_classifier.local_db.clinvar_sqlite.query_bs2_benign_evidence",
               return_value=(True, "Supporting")):
        r = ev.evaluate(_V, _ann("NOVCEP", ac=0))
    assert r.triggered and "expert-panel" in r.evidence


def test_vcep_gene_still_uses_gnomad(tmp_path):
    ev = BS2Evaluator(_cfg(tmp_path))
    assert ev.evaluate(_V, _ann("VCEPGENE", ac=50, nhomalt=2)).triggered
    assert not ev.evaluate(_V, _ann("VCEPGENE", ac=50, nhomalt=1)).triggered


def _cfg_rows(tmp_path: Path, rows: str) -> MagicMock:
    cfg = _cfg(tmp_path)
    cfg.disease_prevalence_tsv.write_text(
        "gene_symbol\tbs2\tinheritance\tbs2_count\tbs2_strength\tbs2_female_only\tbs2_hom_only\tbs2_xl_hom\n"
        + rows, encoding="utf-8")
    return cfg


def test_xl_female_homozygotes(tmp_path):
    # MTM1: >=1 hemizygous male or >=1 homozygous female
    ev = BS2Evaluator(_cfg_rows(tmp_path, "MTM1\tapplicable\tXL\t1\t\t\t\t1\n"
                                         "OTHERXL\tapplicable\tXL\t1\t\t\t\t\n"))
    assert ev.evaluate(_V, _ann("MTM1", ac=2, nhemi=1)).triggered
    assert ev.evaluate(_V, _ann("MTM1", ac=2, nhomalt=1)).triggered
    assert not ev.evaluate(_V, _ann("MTM1", ac=2)).triggered
    # without the flag, female homozygotes do not count
    assert not ev.evaluate(_V, _ann("OTHERXL", ac=2, nhomalt=1)).triggered


def test_hom_only_count_one_and_tiers(tmp_path):
    ev = BS2Evaluator(_cfg_rows(tmp_path, "ETHE1\tapplicable\t\t1\t\t\t1\t\n"
                                         "SCN1B\tapplicable\tAD,AR\t1\t\t\t1\t\n"
                                         "ADA\tapplicable\tAR\t\tStrong:3,Supporting:1\t\t1\t\n"))
    assert ev.evaluate(_V, _ann("ETHE1", ac=2, nhomalt=1)).triggered
    assert not ev.evaluate(_V, _ann("SCN1B", ac=50, nhomalt=0)).triggered   # carriers ignored
    assert ev.evaluate(_V, _ann("SCN1B", ac=50, nhomalt=1)).triggered
    from acmg_classifier.models.enums import CriterionStrength
    assert ev.evaluate(_V, _ann("ADA", ac=4, nhomalt=1)).strength == CriterionStrength.SUPPORTING
    assert ev.evaluate(_V, _ann("ADA", ac=8, nhomalt=3)).strength == CriterionStrength.STRONG
