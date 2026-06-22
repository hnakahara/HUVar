"""PS3 SCV free-text mining patterns (clinvar_builder).

An SCV is counted toward PS3 when the combined comment/description text matches
a damaging-functional-assay indicator (_FUNCTIONAL_POS), does NOT match a
negation (_FUNCTIONAL_NEG), AND — when the claim is quantitative
(_QUANT_FUNCTIONAL: activity / levels / % / protein stability) — cites a PMID
(_PMID_RE). "loss of function" was removed as a positive trigger and "assay"
added.
"""
from acmg_classifier.setup.clinvar_builder import (
    _FUNCTIONAL_POS, _FUNCTIONAL_NEG, _QUANT_FUNCTIONAL, _PMID_RE,
)


def _is_functional(text: str) -> bool:
    """Replicates the build-time counting rule (clinvar_builder, SCV loop)."""
    if not (_FUNCTIONAL_POS.search(text) and not _FUNCTIONAL_NEG.search(text)):
        return False
    if _QUANT_FUNCTIONAL.search(text) and not _PMID_RE.search(text):
        return False
    return True


class TestPositive:
    def test_assay_now_matches(self):
        # Qualitative (splicing) assay — no PMID required.
        assert _is_functional("A minigene assay demonstrated aberrant splicing.")

    def test_quantitative_study_with_pmid_matches(self):
        assert _is_functional(
            "Functional studies showed reduced enzyme activity (PMID: 12345678)."
        )

    def test_abolish_with_pmid_matches(self):
        assert _is_functional(
            "The substitution abolishes catalytic activity in vitro (PMID 9990021)."
        )


class TestQuantitativeRequiresPmid:
    def test_quantitative_without_pmid_excluded(self):
        # Reduced enzyme activity but cited only by author/year → not counted.
        assert not _is_functional(
            "Functional studies demonstrated significantly reduced enzymatic activity "
            "(Takahashi et al. 2006; Li et al. 2014)."
        )

    def test_abolished_activity_without_pmid_excluded(self):
        assert not _is_functional(
            "In vitro assays confirmed that the mutation abolished RPE65 enzymatic activity."
        )

    def test_percent_activity_without_pmid_excluded(self):
        assert not _is_functional("Enzyme activity was reduced to <10% of wild type.")

    def test_protein_stability_without_pmid_excluded(self):
        assert not _is_functional(
            "Western blot showed decreased protein stability for the variant."
        )

    def test_percent_activity_with_pmids_matches(self):
        assert _is_functional(
            "Enzyme activity was reduced to <10% of wild type (PMIDs: 11786058, 34492281)."
        )


class TestLossOfFunctionRemoved:
    def test_bare_loss_of_function_no_longer_fires(self):
        assert not _is_functional(
            "This variant is predicted to cause loss of function of the protein."
        )

    def test_loss_of_function_qualitative_still_fires(self):
        # Qualitative LoF statement from a functional study (no activity numbers) —
        # no PMID required.
        assert _is_functional(
            "Functional studies demonstrate loss of function of the gene product."
        )


class TestNegationGuard:
    def test_assay_with_normal_result_excluded(self):
        assert not _is_functional("The assay showed normal protein function.")

    def test_no_functional_assay_excluded(self):
        assert not _is_functional("There is no functional assay data for this variant.")

    def test_functional_assays_not_reported_excluded(self):
        # ACADVL-style: the assay was never done → not PS3 evidence.
        assert not _is_functional(
            "To our knowledge, functional assays have not been reported for this variant."
        )

    def test_not_reported_in_functional_studies_excluded(self):
        # DCLRE1C-style: "... has not been reported ... or in functional studies."
        assert not _is_functional(
            "To our knowledge, this variant has not been reported in the literature in "
            "individuals affected with SCID/DCLRE1C-related conditions or in functional studies."
        )

    def test_experimental_studies_not_available_excluded(self):
        # GAA-style: "results of experimental studies are not available."
        assert not _is_functional(
            "To our knowledge, this variant has not been reported in the literature in "
            "individuals with Pompe disease, and results of experimental studies are not available."
        )

    def test_qualitative_assay_positive_still_fires(self):
        # Qualitative splicing/recombination assay with a PMID — fires.
        assert _is_functional(
            "An in vitro V(D)J recombination assay showed an abnormal result (PMID: 30622176)."
        )

    def test_not_reported_then_quantitative_clause_with_pmid_fires(self):
        # Negation about individuals, separate clause with a real (cited) result.
        assert _is_functional(
            "This variant has not been reported in patients; however functional studies "
            "confirmed significantly reduced enzymatic activity (PMID: 16754667)."
        )
