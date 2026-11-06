"""VCF sample columns (sample/depth/allele_count/VAF, one row per carrier) and the inference
value columns (scores and cutoffs the criteria actually used)."""
from pathlib import Path

import pytest

from acmg_classifier.config import Config
from acmg_classifier.criteria import insilico_values as iv
from acmg_classifier.criteria.benign import bp4
from acmg_classifier.criteria.pathogenic import pp3
from acmg_classifier.criteria.registry import CriteriaRegistry
from acmg_classifier.io.tsv_writer import _HEADER, _INFERENCE_HEADER, _build_row, write_tsv
from acmg_classifier.io.vcf_reader import read_vcf
from acmg_classifier.models.annotation import (
    AnnotationData, ConsequenceInfo, ESM1bData, FuncVEPData, GnomADData,
)
from acmg_classifier.models.classification import ClassificationResult
from acmg_classifier.models.enums import (
    Assembly, ConsequenceType, InSilicoTool, Pathogenicity,
)
from acmg_classifier.models.variant import SampleCall, VariantRecord

_VCF = (
    "##fileformat=VCFv4.2\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">\n'
    '##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations">\n'
    "##contig=<ID=chr1,length=1000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tblood\ttumor\tmother\n"
    # biallelic: blood het, tumor het (no DP -> sum of AD), mother hom-ref
    "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP\t0/1:6,4:10\t0|1:30,10:.\t0/0:20,0:20\n"
    # multi-allelic: blood 1/2, tumor carries only ALT 2
    "chr1\t200\t.\tC\tT,G\t50\tPASS\t.\tGT:AD:DP\t1/2:2,5,3:10\t0/2:8,0,2:10\t0/0:9,0,0:9\n"
    # AD missing -> FORMAT/AF; missing genotype
    "chr1\t300\t.\tT\tC\t50\tPASS\t.\tGT:AF\t0/1:0.25\t./.:.\t0/0:0\n"
    # nobody carries it
    "chr1\t400\t.\tG\tA\t50\tPASS\t.\tGT\t0/0\t0/0\t./.\n"
    # freebayes style: no AD, AO + DP -> count AO, VAF AO / DP
    "chr1\t500\t.\tA\tT\t50\tPASS\t.\tGT:DP:AO\t0/1:20:5\t0/0:10:0\t./.:.:.\n"
)


@pytest.fixture()
def vcf(tmp_path):
    p = tmp_path / "s.vcf"
    p.write_text(_VCF)
    return p


def test_sample_calls_from_format_fields(vcf):
    recs = {r.key: r for r in read_vcf(vcf, Assembly.GRCH38)}
    calls = {k: [(c.sample, c.depth, c.allele_count, c.vaf) for c in r.sample_calls]
             for k, r in recs.items()}
    assert calls["chr1:100:A:G"] == [("blood", 10, 4, 0.4), ("tumor", 40, 10, 0.25)]
    assert calls["chr1:200:C:T"] == [("blood", 10, 5, 0.5)]
    assert calls["chr1:200:C:G"] == [("blood", 10, 3, 0.3), ("tumor", 10, 2, 0.2)]
    assert calls["chr1:300:T:C"] == [("blood", None, None, 0.25)]
    assert calls["chr1:400:G:A"] == []
    assert calls["chr1:500:A:T"] == [("blood", 20, 5, 0.25)]


def test_one_tsv_row_per_carrier_sample(vcf, tmp_path):
    results = []
    for r in read_vcf(vcf, Assembly.GRCH38):
        results.append(ClassificationResult(
            variant_id=r.key, chrom=r.chrom, pos=r.pos, ref=r.ref, alt=r.alt,
            sample_calls=r.sample_calls, criteria_results=[],
            classification_2015=Pathogenicity.VUS, classification_2015_rules="",
            bayesian_score=0, classification_bayesian=Pathogenicity.VUS))
    out = tmp_path / "o.tsv"
    write_tsv(results, out)
    lines = out.read_text().splitlines()
    head = lines[0].split("\t")
    rows = [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]
    assert "gt" not in head and "dp" not in head
    got = [(r["variant_id"], r["sample"], r["depth"], r["allele_count"], r["vaf"]) for r in rows]
    assert got == [
        ("chr1:100:A:G", "blood", "10", "4", "0.4"),
        ("chr1:100:A:G", "tumor", "40", "10", "0.25"),
        ("chr1:200:C:T", "blood", "10", "5", "0.5"),
        ("chr1:200:C:G", "blood", "10", "3", "0.3"),
        ("chr1:200:C:G", "tumor", "10", "2", "0.2"),
        ("chr1:300:T:C", "blood", "", "", "0.25"),
        ("chr1:400:G:A", "", "", "", ""),            # nobody carries it: one row
        ("chr1:500:A:T", "blood", "20", "5", "0.25"),
    ]
    assert head.index("sample") == head.index("filter") + 1


def test_sites_only_vcf(tmp_path):
    p = tmp_path / "sites.vcf"
    p.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "chr1\t100\t.\tA\tG\t50\tPASS\t.\n")
    (rec,) = read_vcf(p, Assembly.GRCH38)
    assert rec.sample_calls == []


# --- inference values ----------------------------------------------------------

def _missense_ann(**kw):
    return AnnotationData(
        gnomad=kw.pop("gnomad", GnomADData(af=0.0, ac=0, an=100000)),
        consequences=[ConsequenceInfo(
            transcript_id="NM_000059.4", gene_id="ENSG", gene_symbol="BRCA2",
            consequence=ConsequenceType.MISSENSE, biotype="protein_coding",
            is_mane_select=True, protein_position=10, amino_acids="R/C",
            hgvs_c="NM_000059.4:c.28C>T", hgvs_p="NP_000050.3:p.Arg10Cys")],
        **kw)


_V = VariantRecord(chrom="chr13", pos=32316488, ref="C", alt="T", assembly=Assembly.GRCH38)


def _values(cfg, ann):
    res = CriteriaRegistry(cfg).evaluate_all(_V, ann)
    from acmg_classifier.io.tsv_writer import inference_cells
    return dict(zip(_INFERENCE_HEADER, inference_cells(res)))


def test_funcvep_score_and_cutoffs_are_reported():
    cfg = Config(data_dir=Path("./data"), insilico_tool=InSilicoTool.FUNCVEP)
    vals = _values(cfg, _missense_ann(funcvep=FuncVEPData(score=0.93, gene="BRCA2")))
    assert vals["insilico_tool"] == "FuncVEP-CTI"
    assert vals["insilico_score"] == "0.9300"
    assert vals["pp3_insilico_cutoffs"].startswith("Supporting>=0.6903")
    assert "Strong>=0.9109" in vals["pp3_insilico_cutoffs"]
    assert vals["bp4_insilico_cutoffs"].startswith("Supporting<=0.1352")


def test_esm1b_and_frequency_values():
    cfg = Config(data_dir=Path("./data"), insilico_tool=InSilicoTool.ESM1B)
    gd = GnomADData(af=0.02, ac=2000, an=100000, faf95_popmax=0.015, popmax_af=0.02,
                    filter_pass=True)
    vals = _values(cfg, _missense_ann(esm1b=ESM1bData(llr=-3.0), gnomad=gd))
    assert (vals["insilico_tool"], vals["insilico_score"]) == ("ESM1b", "-3.000")
    assert vals["pp3_insilico_cutoffs"] == ("Supporting<=-10.7; Moderate<=-12.2; "
                                            "ThreePoint<=-14; Strong<=-24")
    assert vals["af_metric"].startswith("gnomAD")
    assert vals["af_value"] == "0.015"
    assert vals["ba1_cutoff"] and vals["bs1_cutoff"]
    assert vals["pm2_value"] and vals["pm2_cutoff"]


def test_inference_columns_follow_the_annotation_columns():
    assert _HEADER.index("insilico_tool") == _HEADER.index("repeat_class") + 1


@pytest.mark.parametrize("fn,table,op", [
    (pp3._alphamissense_pp3, iv.AM_PP3, ">="), (pp3._esm1b_pp3, iv.ESM1B_PP3, "<="),
    (bp4._alphamissense_bp4, iv.AM_BP4, "<="), (bp4._esm1b_bp4, iv.ESM1B_BP4, ">="),
])
def test_reported_tiers_match_the_strength_functions(fn, table, op):
    for strength, cutoff in table.items():
        assert fn(cutoff) == strength
        # just outside the cutoff the tier is no longer reached
        eps = 1e-6 if op == "<=" else -1e-6
        assert fn(cutoff + eps) != strength


def test_row_without_call_has_empty_sample_cells():
    r = ClassificationResult(
        variant_id="chr1:1:A:G", chrom="chr1", pos=1, ref="A", alt="G", criteria_results=[],
        classification_2015=Pathogenicity.VUS, classification_2015_rules="",
        bayesian_score=0, classification_bayesian=Pathogenicity.VUS)
    row = dict(zip(_HEADER, _build_row(r)))
    assert (row["sample"], row["depth"], row["allele_count"], row["vaf"]) == ("", "", "", "")
    row = dict(zip(_HEADER, _build_row(
        r, SampleCall(sample="s1", depth=12, allele_count=6, vaf=0.5))))
    cells = (row["sample"], row["depth"], row["allele_count"], row["vaf"])
    assert cells == ("s1", "12", "6", "0.5")
