"""Expert-panel criterion import: clause-level negation and expert-panel SCVs only
(Genome Medicine revision; phrasing taken from ClinGen VCEP ClinVar summaries)."""
import xml.etree.ElementTree as ET

from acmg_classifier.setup.clinvar_builder import _ep_scv_texts, _mine_bs2, _mine_ep_criteria


def test_applied_codes_are_imported():
    t = ("This variant was identified in 4 unrelated individuals with hyperglycemia "
         "(PS4_Moderate; PMID: 1). Functional studies in HEK293 cells showed increased "
         "sensitivity, PS3_Moderate (PMID: 2).")
    assert _mine_ep_criteria(t) == "PS3:Moderate;PS4:Moderate"


def test_criteria_list_after_semicolon():
    t = ("In summary, this variant meets the criteria to be classified as pathogenic based on "
         "the ACMG/AMP criteria applied, as specified by the ClinGen Hearing Loss VCEP; "
         "PM2_Supporting, PP3, PM3_Strong, PP1, PP4.")
    assert _mine_ep_criteria(t) == "PP1:Supporting"


def test_negated_phrasings_are_not_imported():
    for t in (
        "however, PS4_Moderate cannot be applied because this number is below the ClinGen MDEP threshold.",
        "Given that it was only identified in one case, this is below the MDEP VCEP threshold to apply PS4_Moderate.",
        "this data is currently insufficient to establish co-segregation and apply PP1.",
        "PP1, PS2, and PM6 were not assessed due to absence of co-segregation.",
        "these studies do not meet the criteria set forth by the MDEP for the application of PS3 or BS3.",
        "did not meet the OddsPath threshold for PS3_Supporting (> 2.1).",
        "Functional studies have not been conducted for this variant (PS3 not assessed).",
        "with the variant present in the homozygous state (potential PP1_Strong; PMID: 23105016).",
        "No PS4 or PP4 points have been applied due to the variant meeting BA1.",
        "The PP3 and PS3_Supporting codes were not met as PVS1 has been applied instead.",
    ):
        assert _mine_ep_criteria(t) is None, t


def test_negation_of_one_code_does_not_suppress_another():
    t = ("Identified in 4 probands (PS4_Moderate), but segregation data were insufficient "
         "to apply PP1.")
    assert _mine_ep_criteria(t) == "PS4:Moderate"


def test_bs2_clause_level():
    assert _mine_bs2("Observed in 3 homozygous unaffected adults (BS2). BA1 not applicable.") == (1, "Strong")
    assert _mine_bs2("BS2 not met.") == (0, None)
    assert _mine_bs2("BS2_Supporting cannot be applied.") == (0, None)


def test_only_expert_panel_scvs_are_read():
    xml = """<ClinVarSet>
      <ClinVarAssertion><ReviewStatus>criteria provided, single submitter</ReviewStatus>
        <Comment>Lab criteria: PS4_Supporting, PM2.</Comment></ClinVarAssertion>
      <ClinVarAssertion><ReviewStatus>reviewed by expert panel</ReviewStatus>
        <Comment>Applied: BA1.</Comment></ClinVarAssertion>
    </ClinVarSet>"""
    texts = _ep_scv_texts(ET.fromstring(xml))
    assert texts == ["Applied: BA1."]
    assert _mine_ep_criteria(" ".join(texts)) is None
