"""PM2 on the overall gnomAD AF for the InSiGHT MMR genes (pm2_basis=global), and
UTF-8 tabix lines (ClinVar disease names such as "Muir-Torré_syndrome")."""
from pathlib import Path
from unittest.mock import MagicMock

import pysam

from acmg_classifier.criteria.pathogenic.pm2 import PM2Evaluator
from acmg_classifier.criteria.pm2_genes import PM2Spec
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord

_HEADER = "gene_symbol\tpm2_threshold\tpm2_strength\tpm2_basis\tpm2_subpop\tpm2_zygosity\n"
_ROWS = "MLH1\t0.00002\t\tglobal\t\t\nMSH2\t0.00002\t\t\t\t\n"


def _cfg(tmp_path: Path) -> MagicMock:
    p = tmp_path / "disease_prevalence.tsv"
    p.write_text(_HEADER + _ROWS, encoding="utf-8")
    cfg = MagicMock()
    cfg.disease_prevalence_tsv = p
    cfg.gene_inheritance_tsv = tmp_path / "missing_inheritance.tsv"
    return cfg


def _ann(gene: str, **gnomad) -> AnnotationData:
    return AnnotationData(
        gnomad=GnomADData(filter_pass=True, **gnomad),
        consequences=[ConsequenceInfo(
            transcript_id="NM_x", gene_id="ENSG", gene_symbol=gene,
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding")],
    )


def _snv() -> VariantRecord:
    return VariantRecord(chrom="chr3", pos=37048973, ref="C", alt="T", assembly=Assembly.GRCH38)


def test_global_basis_loaded(tmp_path):
    spec = PM2Spec(tmp_path / "x.tsv")
    assert spec.get("MLH1") is None
    spec = PM2Spec(_cfg(tmp_path).disease_prevalence_tsv)
    assert spec.get("MLH1").use_global and not spec.get("MLH1").use_faf
    assert not spec.get("MSH2").use_global


def test_global_basis_uses_overall_af(tmp_path):
    # MLH1 c.2059C>T: overall AF 7.4e-6, GrpMax point AF 7.7e-5, FAF95 3.6e-5.
    ev = PM2Evaluator(_cfg(tmp_path))
    r = ev.evaluate(_snv(), _ann("MLH1", ac=12, af=7.4e-6, popmax_af=7.7e-5,
                                 faf95_popmax=3.6e-5))
    assert r.triggered and r.strength == CriterionStrength.SUPPORTING
    assert "global AF" in r.evidence


def test_default_basis_still_grpmax(tmp_path):
    ev = PM2Evaluator(_cfg(tmp_path))
    r = ev.evaluate(_snv(), _ann("MSH2", ac=12, af=7.4e-6, popmax_af=7.7e-5,
                                 faf95_popmax=3.6e-5))
    assert not r.triggered


def test_clinvar_lines_with_utf8_are_read(tmp_path):
    from acmg_classifier.local_db.clinvar_vcf import query_clinvar_vcf, query_clinvar_vcf_reader
    from acmg_classifier.utils.tabix import TabixReader
    vcf = tmp_path / "clinvar.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "3\t37048948\t3779867\tAGAATGCGCTATGTTCTATTCCATCCG\tA\t.\t.\t"
        "CLNDN=Muir-Torré_syndrome;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic\n"
        "3\t37048973\t90014\tC\tT\t.\t.\t"
        "CLNDN=Lynch_syndrome;CLNREVSTAT=reviewed_by_expert_panel;CLNSIG=Pathogenic\n",
        encoding="utf-8")
    gz = Path(pysam.tabix_index(str(vcf), preset="vcf", force=True))
    for chrom in ("chr3", "3"):
        recs = query_clinvar_vcf(gz, chrom, 37048973, "C", "T")
        assert [(r.variation_id, r.star_rating) for r in recs] == [("90014", 3)]
        recs = query_clinvar_vcf_reader(TabixReader(gz), chrom, 37048973, "C", "T")
        assert [r.variation_id for r in recs] == ["90014"]
    assert query_clinvar_vcf(gz, "chr7", 100, "C", "T") == []
