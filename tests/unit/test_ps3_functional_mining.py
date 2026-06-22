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
