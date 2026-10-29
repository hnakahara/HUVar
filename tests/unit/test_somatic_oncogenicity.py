"""Somatic oncogenicity (ClinGen/CGC/VICC SOP, Horak et al. 2022) with the bundled
NCG gene roles and cancerhotspots.org v2 tables."""
import io
from pathlib import Path

import pytest

from acmg_classifier.config import Config
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.somatic.evaluator import SomaticEvaluator
from acmg_classifier.somatic.models import (
    ManualOncEvidence, OncCriterion as C, Oncogenicity, classify_points,
)
from acmg_classifier.somatic.pipeline import (
    SomaticClassifier, parse_inline_evidence, read_somatic_supplement,
)
from acmg_classifier.somatic.report import HEADER, build_row
from acmg_classifier.somatic.resources import short_protein_change

CT = ConsequenceType
A = ACMGCriterion
S = CriterionStrength


def _cfg(**kw):
    return Config(data_dir=Path("./nonexistent-data"), **kw)


def _ann(gene, csq, pos=None, aa=None, hgvs_p=None, tx="NM_TEST.1", exon=None, gnomad=None):
    c = ConsequenceInfo(transcript_id=tx, gene_id="g", gene_symbol=gene, consequence=csq,
                        biotype="protein_coding", is_mane_select=True, protein_position=pos,
                        amino_acids=aa, hgvs_p=hgvs_p, exon=exon)
    return AnnotationData(consequences=[c],
                          gnomad=gnomad if gnomad is not None else GnomADData(af=0.0, ac=0, an=0))


def _v(chrom="chr7", pos=140753336, ref="A", alt="T"):
    return VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, assembly=Assembly.GRCH38)


def _germ(*met):
    """Germline HUVar results: the given (criterion, strength, evidence) as met."""
    out = {c: CriteriaResult.not_met(c) for c in (A.PP3, A.BP4, A.BP7, A.PM1, A.PS3, A.BS3)}
    for c, st, ev in met:
        out[c] = CriteriaResult.met(c, st, ev)
    return list(out.values())


def _run(ann, germ=(), manual=None, variant=None):
    crit, ctx = SomaticEvaluator(_cfg()).evaluate(variant or _v(), ann, _germ(*germ), manual)
    by = {r.criterion: r for r in crit}
    score = sum(r.points for r in crit)
    return by, score, ctx


def _active(by):
    return {c for c, r in by.items() if r.active}


@pytest.mark.parametrize("score,exp", [
    (10, Oncogenicity.ONCOGENIC), (9, Oncogenicity.LIKELY_ONCOGENIC), (6, Oncogenicity.LIKELY_ONCOGENIC),
    (5, Oncogenicity.VUS), (0, Oncogenicity.VUS), (-1, Oncogenicity.LIKELY_BENIGN),
    (-6, Oncogenicity.LIKELY_BENIGN), (-7, Oncogenicity.BENIGN),
])
def test_sop_categories(score, exp):
    assert classify_points(score) == exp


def test_braf_v600e_hotspot_os3():
    ann = _ann("BRAF", CT.MISSENSE, 600, "V/E", "NP_004324.2:p.Val600Glu")
    by, score, ctx = _run(ann, [(A.PP3, S.MODERATE, "ESM1b")])
    assert _active(by) == {C.OS3, C.OP1, C.OP4}
    assert "897 samples" in by[C.OS3].evidence and "833" in by[C.OS3].evidence
    assert score == 6 and classify_points(score) == Oncogenicity.LIKELY_ONCOGENIC


def test_hotspot_tiers_om3_op3_and_ref_mismatch():
    by, _, _ = _run(_ann("ALK", CT.MISSENSE, 1174, "F/L"))
    assert by[C.OM3].active and not by[C.OS3].met            # 14 samples at residue, 11 F1174L
    by, _, _ = _run(_ann("BRAF", CT.MISSENSE, 600, "V/R"))
    assert by[C.OP3].active and "4 samples" in by[C.OP3].evidence
    by, _, _ = _run(_ann("BRAF", CT.MISSENSE, 600, "A/E"))   # numbering mismatch
    assert not (by[C.OS3].met or by[C.OM3].met or by[C.OP3].met)
    assert "different transcript numbering" in by[C.OS3].evidence


def test_ovs1_tumour_suppressor_only():
    tp53 = _ann("TP53", CT.STOP_GAINED, 196, "R/*", "NP_000537.3:p.Arg196Ter", exon="6/11")
    by, score, _ = _run(tp53, [(A.PP3, S.SUPPORTING, "splice")], variant=_v("chr17"))
    assert by[C.OVS1].active and by[C.OP4].active
    assert by[C.OP1].suppressed                                # null variant: not counted twice
    assert score == 9
    kras = _ann("KRAS", CT.STOP_GAINED, 100, "Q/*", "NP_004976.2:p.Gln100Ter", exon="3/5")
    by, _, _ = _run(kras)
    assert not by[C.OVS1].met and "oncogene" in by[C.OVS1].evidence
    notch1 = _ann("NOTCH1", CT.FRAMESHIFT, 2400, None, "NP_060087.3:p.Pro2400fs", exon="34/34")
    by, _, _ = _run(notch1)
    assert not by[C.OVS1].met and "dual" in by[C.OVS1].evidence


def test_om2_inframe_and_indel_hotspot():
    ann = _ann("TP53", CT.INFRAME_DELETION, 255, None, "NP_000537.3:p.Ile255del")
    by, _, ctx = _run(ann, variant=_v("chr17"))
    assert by[C.OM2].active
    assert by[C.OP3].active and "I255del" in by[C.OP3].evidence   # region 229-292, 9 samples
    assert short_protein_change("NP_000537.3:p.Glu286_Glu287del") == "E286_E287del"
    by, _, _ = _run(_ann("NOTAGENE", CT.INFRAME_DELETION, 10, None, "NP_1.1:p.Lys10del"))
    assert not by[C.OM2].met


def test_frequency_criteria_exclusive():
    common = GnomADData(af=0.08, ac=9000, an=100000, popmax_af=0.2, popmax_pop="afr")
    by, score, _ = _run(_ann("TP53", CT.MISSENSE, 72, "P/R", gnomad=common))
    assert _active(by) == {C.SBVS1} and score == -8
    mid = GnomADData(af=0.01, ac=1000, an=100000, popmax_af=0.03, popmax_pop="eas")
    by, _, _ = _run(_ann("TP53", CT.MISSENSE, 72, "P/R", gnomad=mid))
    assert _active(by) == {C.SBS1}
    rare = GnomADData(af=2e-6, ac=3, an=1500000, popmax_af=4e-6, popmax_pop="nfe")
    by, _, _ = _run(_ann("TP53", CT.MISSENSE, 72, "P/R", gnomad=rare))
    assert by[C.OP4].active and "extremely low" in by[C.OP4].evidence


def test_functional_assays_only_in_tsg_and_strong():
    ps3 = (A.PS3, S.STRONG, "[MaveDB] urn:mavedb:x abnormal")
    by, _, _ = _run(_ann("TP53", CT.MISSENSE, 273, "R/C"), [ps3])
    assert by[C.OS2].active
    by, _, _ = _run(_ann("TP53", CT.MISSENSE, 273, "R/C"), [(A.PS3, S.MODERATE, "[MaveDB] x")])
    assert not by[C.OS2].met and "Strong" in by[C.OS2].evidence
    by, _, _ = _run(_ann("KRAS", CT.MISSENSE, 12, "G/D"), [ps3])
    assert not by[C.OS2].met and "tumour suppressor" in by[C.OS2].evidence
    by, _, _ = _run(_ann("TP53", CT.MISSENSE, 273, "R/C"), [(A.PS3, S.STRONG, "[Text mining] 3 reports")])
    assert not by[C.OS2].met
    bs3 = (A.BS3, S.STRONG, "[MaveDB] normal")
    by, _, _ = _run(_ann("TP53", CT.MISSENSE, 72, "P/R"), [bs3])
    assert by[C.SBS2].active


def test_om1_from_pm1_and_exclusions():
    by, _, _ = _run(_ann("BRAF", CT.MISSENSE, 600, "V/E"), [(A.PM1, S.MODERATE, "hotspot")])
    assert by[C.OM1].suppressed and by[C.OS3].active        # OM1 not with OS3
    by, _, _ = _run(_ann("TP53", CT.MISSENSE, 300, "P/L"), [(A.PM1, S.MODERATE, "DBD")])
    assert by[C.OM1].active


def test_in_silico_benign_and_synonymous():
    by, _, _ = _run(_ann("TP53", CT.MISSENSE, 72, "P/R"), [(A.BP4, S.SUPPORTING, "ESM1b")])
    assert by[C.SBP1].active
    by, _, _ = _run(_ann("TP53", CT.SYNONYMOUS, 72, "P"), [(A.BP7, S.SUPPORTING, "no splice")])
    assert by[C.SBP2].active


def test_manual_evidence_op2_and_override():
    ann = _ann("RB1", CT.MISSENSE, 600, "L/P")
    manual = parse_inline_evidence(["OP2:retinoblastoma", "OP4:not_met:CH suspected"], "k")
    by, _, _ = _run(ann, manual=manual)
    assert by[C.OP2].active and by[C.OP2].manual and "retinoblastoma" in by[C.OP2].evidence
    assert not by[C.OP4].met and by[C.OP4].manual
    by, _, _ = _run(ann)
    assert not by[C.OP2].met
    with pytest.raises(ValueError):
        parse_inline_evidence(["PS3:strong"], "k")


def test_supplement_tsv(tmp_path):
    p = tmp_path / "s.tsv"
    p.write_text("variant_id\tcriterion\tmet\tevidence\nchr1:1:A:G\tOS2\tyes\tPMID 1\n"
                 "chr1:1:A:G\tOP4\tno\tCH\n", encoding="utf-8")
    e = read_somatic_supplement(p)["chr1:1:A:G"]
    assert [(x.criterion, x.met) for x in e] == [(C.OS2, True), (C.OP4, False)]


def test_combination_rules_follow_sop_table2():
    ev = SomaticEvaluator(_cfg())
    from acmg_classifier.somatic.models import OncCriterionResult
    res = {c: OncCriterionResult(criterion=c) for c in C}
    for c in (C.OS1, C.OS3, C.OM4, C.OM1):
        res[c].met = True
    ev._combination_rules(res)
    # OS1 and OS3 both count (SOP); OM1 / OM4 are not used with them.
    assert res[C.OS1].active and res[C.OS3].active
    assert res[C.OM4].suppressed and res[C.OM1].suppressed


def test_classifier_end_to_end_and_tsv_row():
    ann = _ann("BRAF", CT.MISSENSE, 600, "V/E", "NP_004324.2:p.Val600Glu")
    r = SomaticClassifier(_cfg()).classify(_v(), ann, [ManualOncEvidence(
        variant_id="x", criterion=C.OS2, evidence="PMID 12068308")])
    assert C.OS3.value in r.met_criteria and C.OS2.value in r.met_criteria
    assert r.gene_role == "oncogene" and r.classification in (Oncogenicity.LIKELY_ONCOGENIC,
                                                             Oncogenicity.ONCOGENIC)
    row = build_row(r)
    assert len(row) == len(HEADER)
    assert row[HEADER.index("OS3")] == "+4"
    from acmg_classifier.somatic.report import write_somatic_tsv
    import sys
    buf = io.StringIO()
    old, sys.stdout = sys.stdout, buf
    try:
        write_somatic_tsv([r], None)
    finally:
        sys.stdout = old
    assert buf.getvalue().count("\n") == 2


def test_non_cancer_gene_warning():
    r = SomaticClassifier(_cfg()).classify(_v(), _ann("MYH7", CT.STOP_GAINED, 10, "Q/*"))
    assert r.gene_role == "not_cancer_gene" and any("GENE_ROLE" in w for w in r.warnings)


def test_ovs1_nmd_escape_by_truncated_fraction():
    # APC last exon (NMD escape): p.Arg1450Ter removes ~49% of the protein -> OVS1;
    # a truncation in the last 10% is left for manual review.
    early = _ann("APC", CT.STOP_GAINED, 1450, "R/*", "NP_000029.2:p.Arg1450Ter",
                 tx="NM_000038.6", exon="16/16")
    by, _, _ = _run(early, variant=_v("chr5"))
    assert by[C.OVS1].active and "NMD escape" in by[C.OVS1].evidence and "/2843" in by[C.OVS1].evidence
    late = _ann("APC", CT.STOP_GAINED, 2800, "S/*", "NP_000029.2:p.Ser2800Ter",
                tx="NM_000038.6", exon="16/16")
    by, _, _ = _run(late, variant=_v("chr5"))
    assert not by[C.OVS1].met and "review manually" in by[C.OVS1].evidence


# --- ClinVar somatic oncogenicity comparators ------------------------------------

_V600E = ("13961", "7", 140753336, "A", "T", "BRAF", "NP_004324.2:p.Val600Glu", "V600E", 600)
_V600L = ("1964069", "7", 140753335, "C", "A", "BRAF", "NP_004324.2:p.Val600Leu", "V600L", 600)
_MULTI = "criteria provided, multiple submitters"
_SINGLE = "criteria provided, single submitter"


def _clinvar_data(tmp_path, sqlite_rows=(), vcf_rows=(), old_schema=False):
    """ClinVar SQLite (rows: base tuple + (oncogenicity, review status, stars)) and,
    optionally, a tabix ClinVar VCF with ONC fields, under a data dir."""
    import sqlite3

    import pysam

    from acmg_classifier.local_db.clinvar_vcf import onc_star_rating
    from acmg_classifier.setup.clinvar_builder import _CREATE_TABLE

    cfg = Config(data_dir=tmp_path)
    cfg.clinvar_sqlite.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(cfg.clinvar_sqlite)
    ddl = _CREATE_TABLE
    if old_schema:
        ddl = ddl.replace(",\n    oncogenicity TEXT,\n    onc_review_status TEXT,\n"
                          "    onc_stars INTEGER DEFAULT 0", "")
    con.executescript(ddl)
    for base, onc, rs in sqlite_rows:
        cols = ("variation_id, chrom, pos, ref, alt, gene_symbol, hgvs_p, amino_acid_change, "
                "codon_position, clinical_significance, review_status, star_rating")
        vals = list(base) + ["", "", 0]
        if not old_schema:
            cols += ", oncogenicity, onc_review_status, onc_stars"
            vals += [onc, rs, onc_star_rating(rs)]
        con.execute(f"INSERT INTO variants ({cols}) VALUES ({','.join('?' * len(vals))})", vals)
    con.commit()
    con.close()
    if vcf_rows:
        plain = cfg.clinvar_vcf.with_suffix("")
        lines = ["##fileformat=VCFv4.1", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
        lines += ["\t".join(map(str, r)) for r in sorted(vcf_rows, key=lambda r: r[1])]
        plain.write_text("\n".join(lines) + "\n", encoding="utf-8")
        pysam.tabix_compress(str(plain), str(cfg.clinvar_vcf), force=True)
        pysam.tabix_index(str(cfg.clinvar_vcf), preset="vcf", force=True)
    return cfg


def _braf(cfg, germ=()):
    ann = _ann("BRAF", CT.MISSENSE, 600, "V/E", "NP_004324.2:p.Val600Glu",
               gnomad=GnomADData(af=1.24e-6, ac=2, an=1600000, popmax_af=1.1e-5, popmax_pop="sas"))
    crit, _ = SomaticEvaluator(cfg).evaluate(_v(), ann, _germ(*germ))
    return {r.criterion: r for r in crit}, sum(r.points for r in crit)


def test_own_clinvar_oncogenicity_is_shown_but_not_counted(tmp_path):
    cfg = _clinvar_data(tmp_path, [(_V600E, "Oncogenic", _MULTI)])
    by, score = _braf(cfg, [(A.PP3, S.STRONG, "ESM1b")])
    assert not by[C.OS1].met
    assert "this variant is ClinVar Oncogenic" in by[C.OS1].evidence and "not counted" in by[C.OS1].evidence
    assert by[C.OS3].active
    assert by[C.OP4].active and "outlier" in by[C.OP4].evidence      # known oncogenic itself
    assert score == 6 and classify_points(score) == Oncogenicity.LIKELY_ONCOGENIC


def test_os1_from_a_different_nucleotide_change(tmp_path):
    # Same amino-acid change V600E from c.1799_1800delinsAA (another allele).
    other = ("376069", "7", 140753335, "CA", "TT", "BRAF", "NP_004324.2:p.Val600Glu", "V600E", 600)
    cfg = _clinvar_data(tmp_path, [(_V600E, "Oncogenic", _MULTI), (other, "Oncogenic", _MULTI)])
    by, score = _braf(cfg)
    assert by[C.OS1].active and "ClinVar record 376069" in by[C.OS1].evidence
    assert by[C.OS3].active                                        # OS1 and OS3 both count


def test_onc_below_threshold_and_om4(tmp_path):
    cfg = _clinvar_data(tmp_path, [(_V600L, "Likely oncogenic", _MULTI),
                                   (("9", "7", 140753334, "A", "G", "BRAF", "NP_004324.2:p.Val600Ala",
                                     "V600A", 600), "Oncogenic", _SINGLE)])
    by, _ = _braf(cfg)
    assert by[C.OM4].met and "p.Val600Leu" in by[C.OM4].evidence
    assert by[C.OM4].suppressed                                    # OS3 (BRAF V600 hotspot) met
    assert not by[C.OS1].met
    # V600A has a single-submitter record only: listed under OM4 as below the threshold
    # only when OM4 is not met (here V600L meets it).


def test_vcf_fallback_for_old_clinvar_db(tmp_path):
    cfg = _clinvar_data(tmp_path, [(_V600L, None, "")], old_schema=True,
                        vcf_rows=[("7", 140753335, "1964069", "C", "A", ".", ".",
                                   "ONC=Oncogenic;ONCREVSTAT=criteria_provided,_multiple_submitters"),
                                  ("7", 140753336, "13961", "A", "T", ".", ".",
                                   "ONC=Oncogenic;ONCREVSTAT=criteria_provided,_multiple_submitters")])
    by, _ = _braf(cfg)
    assert by[C.OM4].met and "ClinVar record 1964069" in by[C.OM4].evidence
    assert not by[C.OS1].met and "this variant is ClinVar Oncogenic" in by[C.OS1].evidence


def test_builder_reads_rcv_oncogenicity():
    import xml.etree.ElementTree as ET

    from acmg_classifier.setup.clinvar_builder import _rcv_oncogenicity
    rcv = ET.fromstring(
        "<ReferenceClinVarAssertion><Classifications><OncogenicityClassification>"
        "<ReviewStatus>criteria provided, multiple submitters</ReviewStatus>"
        "<Description>Oncogenic</Description></OncogenicityClassification>"
        "</Classifications></ReferenceClinVarAssertion>")
    assert _rcv_oncogenicity(rcv) == ("Oncogenic", "criteria provided, multiple submitters", 2)
    assert _rcv_oncogenicity(ET.fromstring("<ReferenceClinVarAssertion/>")) == ("", "", 0)


@pytest.mark.parametrize("review,met", [(_SINGLE, True), ("no assertion criteria provided", False)])
def test_os1_onc_one_star_is_enough(tmp_path, review, met):
    other = ("376069", "7", 140753335, "CA", "TT", "BRAF", "NP_004324.2:p.Val600Glu", "V600E", 600)
    cfg = _clinvar_data(tmp_path, [(other, "Oncogenic", review)])
    by, _ = _braf(cfg)
    assert by[C.OS1].met is met
    if not met:
        assert "below the star threshold" in by[C.OS1].evidence


def test_om4_onc_still_needs_two_stars(tmp_path):
    cfg = _clinvar_data(tmp_path, [(_V600L, "Oncogenic", _SINGLE)])
    by, _ = _braf(cfg)
    assert not by[C.OM4].met and "below the star threshold" in by[C.OM4].evidence


def test_germline_plp_not_used_for_os1_om4(tmp_path):
    import sqlite3
    cfg = _clinvar_data(tmp_path)
    con = sqlite3.connect(cfg.clinvar_sqlite)
    for vid, pos, ref, alt, hp, aa in (("1", 140753335, "CA", "TT", "p.Val600Glu", "V600E"),
                                       ("2", 140753335, "C", "A", "p.Val600Leu", "V600L")):
        con.execute("INSERT INTO variants (variation_id, chrom, pos, ref, alt, gene_symbol, hgvs_p, "
                    "amino_acid_change, codon_position, clinical_significance, review_status, "
                    "star_rating) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    (vid, "7", pos, ref, alt, "BRAF", "NP_004324.2:" + hp, aa, 600, "Pathogenic",
                     "reviewed by expert panel", 3))
    con.commit()
    con.close()
    by, _ = _braf(cfg)
    assert not by[C.OS1].met and not by[C.OM4].met
