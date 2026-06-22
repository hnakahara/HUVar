"""PS3 SCV free-text mining patterns (clinvar_builder._FUNCTIONAL_POS / _NEG).

An SCV is counted toward PS3 when the combined comment/description text matches
a damaging-functional-assay indicator (_FUNCTIONAL_POS) AND does not match a
negation (_FUNCTIONAL_NEG). "loss of function" was removed as a positive trigger
(it over-fired on classification prose), and "assay" was added.
"""
from acmg_classifier.setup.clinvar_builder import _FUNCTIONAL_POS, _FUNCTIONAL_NEG


def _is_functional(text: str) -> bool:
    """Replicates the build-time counting rule (clinvar_builder.py:573)."""
    return bool(_FUNCTIONAL_POS.search(text) and not _FUNCTIONAL_NEG.search(text))


class TestPositive:
    def test_assay_now_matches(self):
        assert _is_functional("A minigene assay demonstrated aberrant splicing.")

    def test_functional_study_still_matches(self):
        assert _is_functional("Functional studies showed reduced enzyme activity.")

    def test_abolish_still_matches(self):
        assert _is_functional("The substitution abolishes catalytic activity in vitro.")


class TestLossOfFunctionRemoved:
    def test_bare_loss_of_function_no_longer_fires(self):
        # Classification prose mentioning "loss of function" without any other
        # functional-assay indicator must NOT count toward PS3.
        assert not _is_functional(
            "This variant is predicted to cause loss of function of the protein."
        )

    def test_loss_of_function_with_real_assay_still_fires(self):
        # A genuine functional study describing loss of function is still caught
        # via the "functional stud" / "assay" indicators.
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

    def test_experimental_studies_positive_still_fires(self):
        assert _is_functional("Experimental studies showed undetectable enzyme activity.")

    def test_plural_assays_with_real_result_still_fires(self):
        assert _is_functional("In vitro assays confirmed the mutation abolished enzyme activity.")

    def test_not_reported_then_positive_clause_still_fires(self):
        # Negation about individuals, separate clause with a real functional result.
        assert _is_functional(
            "This variant has not been reported in patients; however functional studies "
            "confirmed significantly reduced enzymatic activity."
        )
