"""AMP/ASCO/CAP 2017 tiers: OncoTree -> DOID mapping, CIViC matching and grading."""
import csv
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from acmg_classifier.config import Config
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.enums import Assembly, ConsequenceType as CT
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.somatic.amp import assign, clinvar_calls
from acmg_classifier.somatic.civic_clinical import (
    HEADER, VariantFacts, classify_variant, load_clinical, match_records,
)
from acmg_classifier.somatic.pipeline import SomaticClassifier
from acmg_classifier.somatic.tumor import disease_ontology, relation, tumor_type

DO = disease_ontology()


def test_oncotree_mapping_and_relation():
    skcm, luad, idc = tumor_type("SKCM"), tumor_type("luad"), tumor_type("IDC")
    assert skcm.doids == ("DOID:8923",) and luad.doids == ("DOID:3910",)
    assert idc.method == "parent:BREAST" and idc.doids == ("DOID:1612",)
    assert relation(skcm, "1909", DO) == "same"          # melanoma: more general
    assert relation(skcm, "162", DO) == "same"           # cancer (tumour-agnostic)
    assert relation(luad, "3908", DO) == "same"          # NSCLC
    assert relation(tumor_type("NSCLC"), "3910", DO) == "subtype"
    assert relation(luad, "9256", DO) == "other"         # colorectal cancer
    assert relation(skcm, "", DO) == "unknown"
    with pytest.raises(ValueError):
        tumor_type("NOPE")


@pytest.mark.parametrize("name,hgvs,exp", [
    ("V600E", "NP_004324.2:p.Val600Glu", [("exact", "V600E")]),
    ("Y772_A775DUP", "", [("exact", "Y772_A775dup")]),
    ("V600", "", [("codon", "V600")]),
    ("G12/G13", "", [("codon", "G12"), ("codon", "G13")]),
    ("Q157P/R", "", [("exact", "Q157P"), ("exact", "Q157R")]),
    ("Exon 19 Deletion", "", [("exon", "19:deletion")]),
    ("Exon 14 Skipping Mutation", "", [("exon", "14:skipping")]),
    ("Rare Exon 18-19 Mutation", "", [("exon", "18:mutation"), ("exon", "19:mutation")]),
    ("Loss-of-function", "", [("gene", "lof")]),
    ("Mutation", "", [("gene", "mutation")]),
    ("Amplification", "", []),
])
def test_classify_variant(name, hgvs, exp):
    assert classify_variant(name, hgvs) == exp


def _write(tmp_path, rows):
    p = tmp_path / "shared" / "civic_clinical.tsv"
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, HEADER, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({**{h: "" for h in HEADER}, **r})
    load_clinical.cache_clear()
    return Config(data_dir=tmp_path)


def _ev(match, key, level, doid, disease, rid, etype="Predictive", direction="Supports", **kw):
    return {"gene": kw.pop("gene", "BRAF"), "match_type": match, "key": key, "record": "evidence",
            "record_id": rid, "evidence_type": etype, "direction": direction,
            "significance": "Sensitivity/Response", "level": level, "rating": "4", "doid": doid,
            "disease": disease, "therapies": "Dabrafenib", "citation_id": "1", **kw}


def _ann(gene="BRAF", pos=600, aa="V/E", hp="NP_004324.2:p.Val600Glu", csq=CT.MISSENSE,
         exon="15/18", af=0.0, popmax=None):
    c = ConsequenceInfo(transcript_id="NM_X.1", gene_id="g", gene_symbol=gene, consequence=csq,
                        biotype="protein_coding", is_mane_select=True, protein_position=pos,
                        amino_acids=aa, hgvs_p=hp, exon=exon)
    return AnnotationData(consequences=[c], gnomad=GnomADData(af=af, ac=int(af * 1e5), an=100000,
                                                              popmax_af=popmax))


_V = VariantRecord(chrom="chr7", pos=140753336, ref="A", alt="T", assembly=Assembly.GRCH38)


def test_tiers_by_tumour_type(tmp_path):
    cfg = _write(tmp_path, [_ev("exact", "V600E", "A", "1909", "Melanoma", "1"),
                            _ev("codon", "V600", "B", "3908", "Lung Non-small Cell Carcinoma", "2"),
                            _ev("exact", "V600E", "D", "9256", "Colorectal Cancer", "3"),
                            _ev("exact", "V600E", "A", "1909", "Melanoma", "4",
                                direction="Does Not Support")])
    r = SomaticClassifier(cfg, "SKCM").classify(_V, _ann())
    assert (r.amp_tier, r.amp_level) == ("Tier I", "A") and "EID1" in r.amp_evidence
    r = SomaticClassifier(cfg, "LUAD").classify(_V, _ann())
    assert (r.amp_tier, r.amp_level) == ("Tier I", "B") and "residue V600" in r.amp_evidence
    r = SomaticClassifier(cfg, "PAAD").classify(_V, _ann())
    assert (r.amp_tier, r.amp_level) == ("Tier II", "C")      # level A/B in other tumours
    r = SomaticClassifier(cfg, None).classify(_V, _ann())
    assert r.amp_tier == "" and r.amp_candidates and "tier not assigned" in r.amp_note


def test_assertion_amp_category_and_cap(tmp_path):
    cfg = _write(tmp_path, [{"gene": "BRAF", "match_type": "exact", "key": "V600E",
                             "record": "assertion", "record_id": "7", "evidence_type": "Predictive",
                             "direction": "Supports", "significance": "Sensitivity/Response",
                             "doid": "1909", "disease": "Melanoma", "amp_category": "Tier I - Level A",
                             "nccn_guideline": "Melanoma", "regulatory_approval": "true"}])
    r = SomaticClassifier(cfg, "MEL").classify(_V, _ann())
    assert (r.amp_tier, r.amp_level) == ("Tier I", "A") and "AID7" in r.amp_evidence
    assert "NCCN Melanoma" in r.amp_evidence
    r = SomaticClassifier(cfg, "LUAD").classify(_V, _ann())
    assert (r.amp_tier, r.amp_level) == ("Tier II", "C")


def test_gene_and_exon_level_need_oncogenic(tmp_path):
    cfg = _write(tmp_path, [_ev("gene", "mutation", "A", "3908", "NSCLC", "9", gene="EGFR"),
                            _ev("exon", "19:deletion", "A", "3908", "NSCLC", "10", gene="EGFR")])
    db = load_clinical(cfg.civic_clinical_tsv)
    vus = VariantFacts("EGFR", "A123T", "A", 123, "missense_variant", 3, None, False, False)
    assert match_records(db, vus) == []
    onc = VariantFacts("EGFR", "A123T", "A", 123, "missense_variant", 3, None, False, True)
    assert [h[1] for h in match_records(db, onc)] == ["gene-level mutation (variant classified oncogenic)"]
    dl = VariantFacts("EGFR", "E746_A750del", "E", 746, "inframe_deletion", 19, None, False, False)
    assert [h[1] for h in match_records(db, dl)] == ["exon 19 deletion"]


def test_tier_iii_and_iv(tmp_path):
    cfg = _write(tmp_path, [])
    r = SomaticClassifier(cfg, "SKCM").classify(_V, _ann(gene="TP53", pos=72, aa="P/R", hp="p.Pro72Arg"))
    assert r.amp_tier == "Tier III"
    r = SomaticClassifier(cfg, "SKCM").classify(
        _V, _ann(gene="TP53", pos=72, aa="P/R", hp="p.Pro72Arg", af=0.3, popmax=0.6))
    assert r.amp_tier == "Tier IV"


def test_clinvar_somatic_impact_calls():
    tt = tumor_type("SKCM")
    rec = [{"variation_id": "13961", "sci": "Tier I - Strong", "sci_review_status":
            "criteria provided, single submitter", "sci_stars": 1,
            "sci_type": "therapeutic; sensitivity/response", "conditions": "Melanoma"}]
    c = clinvar_calls(rec, tt, DO)
    assert c[0].tier == "I" and c[0].category == "therapeutic"
    c = clinvar_calls(rec, tumor_type("LUAD"), DO)
    assert (c[0].tier, c[0].level) == ("II", "C")
    res = assign(clinvar_calls(rec, tt, DO), tt, DO, benign=False, common=False)
    assert res.tier == "Tier I" and "ClinVar" in res.note


def test_builder_reads_somatic_clinical_impact():
    from acmg_classifier.setup.clinvar_builder import _rcv_conditions, _rcv_somatic_impact
    rcv = ET.fromstring(
        "<ReferenceClinVarAssertion><Classifications><SomaticClinicalImpact>"
        "<ReviewStatus>criteria provided, single submitter</ReviewStatus>"
        "<Description ClinicalImpactAssertionType=\"therapeutic\" "
        "ClinicalImpactClinicalSignificance=\"sensitivity/response\">Tier I - Strong</Description>"
        "</SomaticClinicalImpact></Classifications><TraitSet><Trait><Name>"
        "<ElementValue Type=\"Preferred\">Melanoma</ElementValue></Name></Trait></TraitSet>"
        "</ReferenceClinVarAssertion>")
    assert _rcv_somatic_impact(rcv) == ("Tier I - Strong", "criteria provided, single submitter", 1,
                                        "therapeutic; sensitivity/response")
    assert _rcv_conditions(rcv) == "Melanoma"


def test_packaged_tables_present():
    root = Path(__file__).resolve().parents[2] / "resources" / "shared"
    assert (root / "oncotree_doid.tsv").exists() and (root / "doid_hierarchy.tsv").exists()


def test_report_keeps_bracketed_text_and_orders_evidence(tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("COLUMNS", "400")
    from acmg_classifier.somatic.report import print_somatic_report
    cfg = _write(tmp_path, [_ev("codon", "V600", "B", "1909", "Melanoma", "2", rating="5"),
                            _ev("exact", "V600E", "A", "1909", "Melanoma", "1", rating="3")])
    r = SomaticClassifier(cfg, "SKCM").classify(_V, _ann())
    assert r.amp_evidence.index("EID1 ") < r.amp_evidence.index("EID2 ")   # level A first
    print_somatic_report(r)
    out = capsys.readouterr().out
    assert "[residue V600]" in out and "[exact]" in out
