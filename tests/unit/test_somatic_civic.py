"""OS2 / SBS2 from CIViC functional evidence and oncogenic assertions."""
from pathlib import Path

from acmg_classifier.config import Config
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo, GnomADData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, Assembly, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.somatic.civic import build_civic_table, civic_call, load_civic
from acmg_classifier.somatic.evaluator import SomaticEvaluator
from acmg_classifier.somatic.models import OncCriterion as C, Oncogenicity, classify_points

A, S, CT = ACMGCriterion, CriterionStrength, ConsequenceType

_EV_HDR = ["molecular_profile", "molecular_profile_id", "disease", "doid", "phenotypes", "therapies",
           "therapy_interaction_type", "evidence_type", "evidence_direction", "evidence_level",
           "significance", "evidence_statement", "citation_id", "source_type", "asco_abstract_id",
           "citation", "nct_ids", "rating", "evidence_status", "evidence_id", "variant_origin",
           "last_review_date", "evidence_civic_url", "molecular_profile_civic_url", "is_flagged"]
_VAR_HDR = ["variant_id", "variant_civic_url", "feature_type", "feature_id", "feature_name",
            "feature_civic_url", "variant", "variant_aliases", "is_flagged", "variant_groups",
            "variant_types", "single_variant_molecular_profile_id", "last_review_date", "gene",
            "entrez_id", "chromosome", "start", "stop", "reference_bases", "variant_bases",
            "representative_transcript", "ensembl_version", "reference_build", "hgvs_descriptions",
            "allele_registry_id", "clinvar_ids"]
_AS_HDR = ["molecular_profile", "molecular_profile_id", "disease", "doid", "phenotypes", "therapies",
           "assertion_type", "assertion_direction", "significance", "acmg_codes", "amp_category",
           "clingen_codes", "nccn_guideline", "nccn_guideline_version", "regulatory_approval",
           "fda_companion_test", "assertion_summary", "assertion_description", "assertion_id",
           "evidence_item_ids", "last_review_date", "assertion_civic_url",
           "evidence_items_civic_url", "molecular_profile_civic_url", "is_flagged"]


def _tsv(path: Path, header, rows):
    lines = ["\t".join(header)] + ["\t".join(str(r.get(h, "")) for h in header) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _ev(mp, sig, pmid, rating, eid, etype="Functional", direction="Supports"):
    return {"molecular_profile_id": mp, "evidence_type": etype, "evidence_direction": direction,
            "evidence_level": "D", "significance": sig, "citation_id": pmid, "rating": rating,
            "evidence_status": "accepted", "evidence_id": eid, "variant_origin": "Somatic"}


def _civic(tmp_path, evidence, assertions=()):
    variants = [
        {"variant_id": 12, "feature_type": "Gene", "variant": "V600E", "gene": "BRAF",
         "single_variant_molecular_profile_id": 12},
        {"variant_id": 30, "feature_type": "Gene", "variant": "R273C", "gene": "TP53",
         "single_variant_molecular_profile_id": 30},
        {"variant_id": 40, "feature_type": "Gene", "variant": "EXON 19 DELETION", "gene": "EGFR",
         "single_variant_molecular_profile_id": 40,
         "hgvs_descriptions": "NP_005219.2:p.Glu746_Ala750del"},
    ]
    cfg = Config(data_dir=tmp_path)
    info = build_civic_table(_tsv(tmp_path / "ev.tsv", _EV_HDR, evidence),
                             _tsv(tmp_path / "var.tsv", _VAR_HDR, variants),
                             _tsv(tmp_path / "as.tsv", _AS_HDR, assertions),
                             cfg.civic_somatic_tsv, release="01-Oct-2026")
    load_civic.cache_clear()
    return cfg, info


def _ann(gene, pos, aa, hgvs_p, csq=CT.MISSENSE, tx="NM_TEST.1"):
    c = ConsequenceInfo(transcript_id=tx, gene_id="g", gene_symbol=gene, consequence=csq,
                        biotype="protein_coding", is_mane_select=True, protein_position=pos,
                        amino_acids=aa, hgvs_p=hgvs_p)
    return AnnotationData(consequences=[c], gnomad=GnomADData(af=1.24e-6, ac=2, an=1600000,
                                                              popmax_af=1.1e-5, popmax_pop="sas"))


def _run(cfg, ann, germ=(), chrom="chr7", pos=140753336):
    g = {c: CriteriaResult.not_met(c) for c in (A.PP3, A.BP4, A.BP7, A.PM1, A.PS3, A.BS3)}
    for c, st, ev in germ:
        g[c] = CriteriaResult.met(c, st, ev)
    v = VariantRecord(chrom=chrom, pos=pos, ref="A", alt="T", assembly=Assembly.GRCH38)
    crit, _ = SomaticEvaluator(cfg).evaluate(v, ann, list(g.values()))
    return {r.criterion: r for r in crit}, sum(r.points for r in crit)


_BRAF = ("BRAF", 600, "V/E", "NP_004324.2:p.Val600Glu")


def test_build_keeps_functional_items_and_coded_assertions(tmp_path):
    _, info = _civic(tmp_path, [
        _ev(12, "Gain of Function", "12068308", 4, 1),
        _ev(12, "Sensitivity/Response", "1", 5, 2, etype="Predictive"),     # not functional
        _ev(99, "Gain of Function", "3", 5, 3),                              # complex profile
        _ev(40, "Gain of Function", "4", 4, 4),
    ], [{"molecular_profile_id": 30, "assertion_type": "Oncogenic", "assertion_direction": "Supports",
         "significance": "Oncogenic", "clingen_codes": "OS2,OS3,OP4", "assertion_id": 7}])
    assert info["evidence_items"] == 2 and info["assertions"] == 1


def test_braf_v600e_gain_of_function_gives_os2(tmp_path):
    cfg, _ = _civic(tmp_path, [_ev(12, "Gain of Function", "12068308", 4, 1),
                               _ev(12, "Gain of Function", "15035987", 3, 2)])
    by, score = _run(cfg, _ann(*_BRAF), [(A.PP3, S.STRONG, "ESM1b")])
    assert by[C.OS2].active and "EID1" in by[C.OS2].evidence and "2 publications" in by[C.OS2].evidence
    assert score == 10 and classify_points(score) == Oncogenicity.ONCOGENIC   # OS2+OS3+OP1+OP4


def test_threshold_single_source_or_low_rating(tmp_path):
    cfg, _ = _civic(tmp_path, [_ev(12, "Gain of Function", "12068308", 4, 1),
                               _ev(12, "Gain of Function", "12068308", 5, 2),   # same paper
                               _ev(12, "Gain of Function", "999", 2, 3)])        # rating 2
    by, _ = _run(cfg, _ann(*_BRAF))
    assert not by[C.OS2].met and "below the threshold" in by[C.OS2].evidence


def test_role_direction_and_conflict():
    ev = lambda sig, pmid, r=4: {"record": "evidence", "record_id": pmid, "significance": sig,  # noqa: E731
                                 "direction": "Supports", "rating": r, "citation_id": pmid, "codes": ""}
    lof = [ev("Loss of Function", "1"), ev("Loss of Function", "2")]
    assert civic_call(lof, "TSG").criterion == "OS2"
    assert civic_call(lof, "oncogene").criterion is None          # LoF does not fit an oncogene
    neutral = [ev("Unaltered Function", "3"), ev("Unaltered Function", "4")]
    assert civic_call(neutral, "oncogene").criterion == "SBS2"
    gof = [ev("Gain of Function", "5"), ev("Gain of Function", "6")]
    assert civic_call(gof + neutral, "oncogene").criterion is None   # conflict
    assert civic_call(gof, None).criterion is None                   # unknown gene role


def test_assertion_coded_os2(tmp_path):
    cfg, _ = _civic(tmp_path, [], [{"molecular_profile_id": 30, "assertion_type": "Oncogenic",
                                    "assertion_direction": "Supports", "significance": "Oncogenic",
                                    "clingen_codes": "OS2,OM1", "assertion_id": 7}])
    by, _ = _run(cfg, _ann("TP53", 273, "R/C", "NP_000537.3:p.Arg273Cys"), chrom="chr17", pos=1)
    assert by[C.OS2].active and "AID7" in by[C.OS2].evidence


def test_inframe_deletion_key(tmp_path):
    cfg, _ = _civic(tmp_path, [_ev(40, "Gain of Function", "1", 4, 1), _ev(40, "Gain of Function", "2", 4, 2)])
    by, _ = _run(cfg, _ann("EGFR", 746, None, "NP_005219.2:p.Glu746_Ala750del", CT.INFRAME_DELETION),
                 chrom="chr7", pos=55174772)
    assert by[C.OS2].active


def test_conflict_with_calibrated_assay(tmp_path):
    cfg, _ = _civic(tmp_path, [_ev(30, "Unaltered Function", "1", 4, 1),
                               _ev(30, "Unaltered Function", "2", 4, 2)])
    by, _ = _run(cfg, _ann("TP53", 273, "R/C", "NP_000537.3:p.Arg273Cys"),
                 [(A.PS3, S.STRONG, "[MaveDB] urn:x abnormal")], chrom="chr17", pos=1)
    assert not by[C.OS2].met and not by[C.SBS2].met and "conflicts" in by[C.OS2].evidence


def test_table_not_installed_note(tmp_path):
    by, _ = _run(Config(data_dir=tmp_path), _ann(*_BRAF))
    assert not by[C.OS2].met and "CIViC table not installed" in by[C.OS2].evidence
