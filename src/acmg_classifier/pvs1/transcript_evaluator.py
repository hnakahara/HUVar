"""Evaluate transcript-level properties needed for PVS1 decision tree."""
from __future__ import annotations
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo


def has_alternative_transcript_rescue(annotation: AnnotationData) -> bool:
    """
    Return True if an alternative MANE/canonical transcript does NOT carry the LoF consequence.

    When a variant causes LoF on one transcript but is tolerated on another clinically
    relevant transcript, PVS1 strength is reduced (e.g., VeryStrong -> Strong).
    """
    primary = annotation.primary_consequence
    if primary is None:
        return False

    from acmg_classifier.models.enums import ConsequenceType
    lof_consequences = {
        ConsequenceType.FRAMESHIFT,
        ConsequenceType.STOP_GAINED,
        ConsequenceType.SPLICE_ACCEPTOR,
        ConsequenceType.SPLICE_DONOR,
        ConsequenceType.START_LOST,
        ConsequenceType.TRANSCRIPT_ABLATION,
    }

    # Count LoF transcripts vs MANE/canonical transcripts that escape LoF.
    # A "rescue" exists only when BOTH conditions hold: at least one
    # transcript carries the LoF, AND at least one clinically-relevant
    # (MANE/canonical) transcript carries a non-LoF consequence. We do not
    # rescue against arbitrary minor isoforms — those are unlikely to
    # produce enough protein to mitigate haploinsufficiency.
    lof_transcripts = 0
    non_lof_mane_canonical = 0

    for c in annotation.consequences:
        if c.consequence in lof_consequences:
            lof_transcripts += 1
        elif c.is_mane_select or c.is_canonical:
            non_lof_mane_canonical += 1

    return non_lof_mane_canonical > 0 and lof_transcripts > 0


_MIN_PLP_NULL = 3       # ClinVar P/LP null variants establishing LoF mechanism
_LOEUF_INTOLERANT = 0.35  # gnomAD LOEUF below which the gene is LoF-constrained


def gene_has_lof_mechanism(
    consequence: ConsequenceInfo,
    gnomad_loeuf: float | None,
    clinvar_plp_null: int = 0,
) -> bool:
    """Heuristic: is loss-of-function an established disease mechanism for the gene?

    Primary signal — the gene already has several P/LP null (nonsense/frameshift)
    variants in ClinVar (cf. Franklin's "pathogenic null variants reported").
    Many bona-fide LoF disease genes (tumour suppressors, recessive genes) are
    NOT population-constrained, so gnomAD LOEUF alone misclassifies them; LOEUF
    is therefore only a secondary signal.

    Established when: >=3 P/LP null ClinVar variants OR LOEUF < 0.35.
    When neither signal is present (no ClinVar nulls AND no LOEUF), LoF is treated
    as NOT established and PVS1 is not applied — consistent with ClinGen/Franklin.
    """
    if clinvar_plp_null >= _MIN_PLP_NULL:
        return True
    if gnomad_loeuf is not None and gnomad_loeuf < _LOEUF_INTOLERANT:
        return True
    return False
