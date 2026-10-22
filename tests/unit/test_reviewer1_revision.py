"""Reviewer 1 (Genome Medicine revision): unsupported alleles, large indels,
gene-disease validity cap exemption for VCEP genes, PS4 vs PM3, PS1 star threshold."""
from pathlib import Path
from types import SimpleNamespace

import pytest

from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import (
    ACMGCriterion as C, Assembly, ConsequenceType, CriterionStrength as S, Pathogenicity,
)
from acmg_classifier.models.variant import (
    LARGE_INDEL_BP, VariantRecord, large_indel_warning, unsupported_allele_reason,
)


# --- R1-1: variant types -------------------------------------------------------

@pytest.mark.parametrize("ref,alt,ok", [
    ("A", "G", True), ("AT", "A", True), ("A", "ATTT", True), ("AC", "GT", True),
    ("A", "<DEL>", False), ("A", "<DUP:TANDEM>", False), ("A", "<CNV>", False),
    ("A", "A[chr2:100[", False), ("A", "]chr2:100]A", False), ("A", "*", False),
    ("A", ".", False), ("A", "", False),
])
def test_unsupported_allele_reason(ref, alt, ok):
    assert (unsupported_allele_reason(ref, alt) is None) is ok


def test_large_indel_warning():
    assert large_indel_warning("A", "A" + "T" * LARGE_INDEL_BP) is not None
    assert large_indel_warning("A", "A" + "T" * (LARGE_INDEL_BP - 1)) is None
    assert large_indel_warning("A" * 60, "T" * 60) is None  # same-length MNV-like: no flag


def test_write_skipped_reasons(tmp_path):
    from acmg_classifier.io.tsv_writer import write_skipped
    recs = [VariantRecord(chrom="1", pos=10, ref="A", alt="<DEL>", assembly=Assembly.GRCH38),
            VariantRecord(chrom="1", pos=11, ref="A", alt=".", assembly=Assembly.GRCH38)]
    out = tmp_path / "x.not_annotated.tsv"
    write_skipped(recs, out, reasons={r.key: unsupported_allele_reason(r.ref, r.alt) for r in recs})
    text = out.read_text()
    assert "symbolic ALT allele <DEL>" in text and "ALT='.'" in text


def test_annotate_one_rejects_symbolic():
    from acmg_classifier.pipeline.pipeline import annotate_one
    with pytest.raises(ValueError, match="symbolic"):
        annotate_one("chr1", 100, "A", "<DEL>", SimpleNamespace(assembly=Assembly.GRCH38))


# --- R1-3: gene-disease validity cap, VCEP genes exempt ---------------------------

def _gdv(tmp_path):
    p = tmp_path / "gdv.csv"
    p.write_text('"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP",'
                 '"CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n'
                 '"VCEPGENE","HGNC:1","other disease","MONDO:1","AD","SOP9","Refuted","x","2020","g"\n'
                 '"WEAK1","HGNC:2","d","MONDO:2","AD","SOP9","Limited","x","2020","g"\n',
                 encoding="utf-8")
    return p


def _dp(tmp_path):
    p = tmp_path / "dp.tsv"
    p.write_text("gene_symbol\tinheritance\tsource_vcep\nVCEPGENE\tAD\tSome Expert Panel\n", encoding="utf-8")
    return p


def test_gdv_cap_exempts_vcep_genes(tmp_path):
    from acmg_classifier.classification.gene_validity import apply_gene_validity, _vcep_genes
    from acmg_classifier.local_db.clingen_gene_db import load_gene_validity
    load_gene_validity.cache_clear()
    _vcep_genes.cache_clear()
    cfg = SimpleNamespace(clingen_gene_validity_csv=_gdv(tmp_path), disease_prevalence_tsv=_dp(tmp_path))
    a, b, w = apply_gene_validity("VCEPGENE", Pathogenicity.PATHOGENIC, Pathogenicity.LIKELY_PATHOGENIC, cfg)
    assert a == Pathogenicity.PATHOGENIC and b == Pathogenicity.LIKELY_PATHOGENIC
    assert w and w[0].startswith("GENE_VALIDITY_NOTE") and "Some Expert Panel" in w[0]
    a, b, w = apply_gene_validity("WEAK1", Pathogenicity.PATHOGENIC, Pathogenicity.PATHOGENIC, cfg)
    assert a == b == Pathogenicity.VUS and "capped at VUS" in w[0]


# --- R1-5: PS4 vs PM3 --------------------------------------------------------------

def _ann(gene):
    return AnnotationData(consequences=[ConsequenceInfo(
        transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
        consequence=ConsequenceType.MISSENSE, biotype="protein_coding")])


def test_ps4_suppressed_when_pm3_met_non_vcep(tmp_path):
    from acmg_classifier.classification.gene_validity import _vcep_genes
    from acmg_classifier.criteria.registry import _apply_ps4_pm3_exclusion
    _vcep_genes.cache_clear()
    cfg = SimpleNamespace(disease_prevalence_tsv=_dp(tmp_path))
    res = [CriteriaResult.met(C.PS4, S.STRONG, "[manual] cases"), CriteriaResult.met(C.PM3, S.MODERATE, "[manual] trans")]
    _apply_ps4_pm3_exclusion(res, _ann("OTHERGENE"), cfg)
    assert res[0].suppressed and "because PM3 is met" in res[0].evidence and not res[1].suppressed


def test_ps4_kept_for_vcep_gene_or_without_pm3(tmp_path):
    from acmg_classifier.classification.gene_validity import _vcep_genes
    from acmg_classifier.criteria.registry import _apply_ps4_pm3_exclusion
    _vcep_genes.cache_clear()
    cfg = SimpleNamespace(disease_prevalence_tsv=_dp(tmp_path))
    res = [CriteriaResult.met(C.PS4, S.STRONG, "x"), CriteriaResult.met(C.PM3, S.MODERATE, "y")]
    _apply_ps4_pm3_exclusion(res, _ann("VCEPGENE"), cfg)
    assert not res[0].suppressed
    res = [CriteriaResult.met(C.PS4, S.STRONG, "x"), CriteriaResult.not_met(C.PM3, "none")]
    _apply_ps4_pm3_exclusion(res, _ann("OTHERGENE"), cfg)
    assert not res[0].suppressed


# --- R1-7: PS1 star threshold --------------------------------------------------------

def test_ps1_min_stars_config():
    from acmg_classifier.config import Config
    from acmg_classifier.criteria.pathogenic.ps1 import PS1Evaluator
    assert Config(data_dir=Path("./data")).ps1_min_stars == 2
    ev = PS1Evaluator.__new__(PS1Evaluator)
    ev._cfg = SimpleNamespace(ps1_min_stars=2)
    assert ev._min_stars == 2
    ev._cfg = SimpleNamespace()
    assert ev._min_stars == 2
    ev._cfg = SimpleNamespace(ps1_min_stars=1)
    assert ev._min_stars == 1
