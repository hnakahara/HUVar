"""Evaluate transcript- and gene-level properties needed for the PVS1 decision tree."""
from __future__ import annotations
from acmg_classifier.models.annotation import AnnotationData, ConsequenceInfo


def has_alternative_transcript_rescue(annotation: AnnotationData) -> bool:
    """
    Return True if an alternative MANE/canonical protein-coding transcript of the
    same gene does NOT carry the LoF consequence.

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

    # A "rescue" exists only when BOTH conditions hold: at least one transcript
    # carries the LoF, AND at least one clinically-relevant (MANE/canonical)
    # protein-coding transcript carries a non-LoF consequence. Minor isoforms
    # are ignored. Only transcripts of the SAME gene as the primary consequence
    # count: VEP also reports transcripts of neighbouring genes within 5 kb
    # (upstream/downstream_gene_variant, e.g. RB1-DT next to RB1), which cannot
    # rescue a loss of function of the gene under evaluation.
    gene = primary.gene_symbol
    lof_transcripts = 0
    non_lof_mane_canonical = 0

    for c in annotation.consequences:
        if c.gene_symbol != gene:
            continue
        if c.consequence in lof_consequences:
            lof_transcripts += 1
        elif (c.is_mane_select or c.is_canonical) and c.biotype == "protein_coding":
            non_lof_mane_canonical += 1

    return non_lof_mane_canonical > 0 and lof_transcripts > 0


# ClinVar route of the LoF-mechanism question (genes without a VCEP PVS1
# decision or a ClinGen HI score of 3): null variants (nonsense, frameshift,
# canonical +/-1,2 splice; ClinVar review status >=1 star) make up at least
# MIN_NULL_FRACTION of the gene's P/LP variants, with at least MIN_PLP P/LP
# variants in total. Selected against the VCEP PVS1 decisions and ClinGen HI
# scores (Genome Medicine revision; reviewer2/lof_rule_comparison.py): unlike a
# fixed count of null variants, the proportion excludes genes whose pathogenic
# variants are predominantly missense (gain-of-function / dominant-negative).
# Config.pvs1_min_plp / Config.pvs1_min_null_fraction override the defaults.
MIN_PLP = 5
MIN_NULL_FRACTION = 0.20


def gene_has_lof_mechanism(
    consequence: ConsequenceInfo | None,
    gnomad_loeuf: float | None = None,
    clinvar_plp_null: int = 0,
    hi_score: str | None = None,
    clinvar_plp_total: int = 0,
) -> bool:
    """Is loss-of-function an established disease mechanism for the gene?

    Established when the gene has a ClinGen Dosage Sensitivity haploinsufficiency
    score of 3, or when null variants make up >= 20% of the gene's ClinVar P/LP
    variants (>= 5 P/LP variants). (A VCEP that explicitly applies PVS1 is
    handled by the caller.)

    gnomAD LOEUF is NOT used: population constraint alone does not establish
    LoF as a disease mechanism (ACGS 2023), and recessive LoF genes are rarely
    constrained. The argument is retained for reporting.
    """
    return lof_mechanism_reason(gnomad_loeuf, clinvar_plp_null, hi_score, clinvar_plp_total)[0]


def lof_mechanism_reason(
    gnomad_loeuf: float | None,
    clinvar_plp_null: int,
    hi_score: str | None,
    clinvar_plp_total: int = 0,
    min_plp: int = MIN_PLP,
    min_null_fraction: float = MIN_NULL_FRACTION,
) -> tuple[bool, str]:
    """``(established, explanation)`` — see :func:`gene_has_lof_mechanism`."""
    loeuf_note = f"; LOEUF={gnomad_loeuf:.2f} (reference only)" if gnomad_loeuf is not None else ""
    if (hi_score or "").strip() == "3":
        return True, f"ClinGen haploinsufficiency score 3{loeuf_note}"
    total = max(clinvar_plp_total, clinvar_plp_null)
    frac = clinvar_plp_null / total if total else 0.0
    counts = (f"{clinvar_plp_null}/{total} ClinVar P/LP variants are null "
              f"({100 * frac:.0f}%)")
    if total >= min_plp and clinvar_plp_null > 0 and frac >= min_null_fraction:
        return True, (f"{counts}: >= {100 * min_null_fraction:.0f}% of >= {min_plp} "
                      f"P/LP variants{loeuf_note}")
    hi_note = f"ClinGen HI score {hi_score}" if hi_score else "no ClinGen HI score"
    why = (f"fewer than {min_plp} P/LP variants" if total < min_plp
           else f"< {100 * min_null_fraction:.0f}% null")
    return False, (
        f"LoF mechanism not established: {hi_note}, {counts} ({why}){loeuf_note}"
    )
