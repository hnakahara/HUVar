"""
PS4 -- variant prevalence in affected individuals significantly higher than controls.

Proband-count-only implementation (no case-control / OR data required).

Adjusted PS4 weightings (ClinGen SVI 2019, ACMG 2015):
- Strong       (PS4):            >=10 unrelated affected probands
- Moderate     (PS4_Moderate):   6-9 unrelated affected probands
- Supporting   (PS4_Supporting): 2-5 unrelated affected probands

Each ClinVar SCV submission that classified the variant as P/LP and reported
AffectedStatus="yes" is approximated as one unrelated proband. Affected
observations from Benign/VUS submitters are incidental findings and are NOT
counted (see clinvar_builder._parse_clinvarset). Both rarity (FAF95_popmax <
0.0001 or AC=0) AND >=2 affected observations are required to trigger any
strength.
"""
from __future__ import annotations
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_FAF95_RARE = 0.0001


def _ps4_strength(n_affected: int) -> CriterionStrength | None:
    if n_affected >= 10:
        return CriterionStrength.STRONG
    if n_affected >= 6:
        return CriterionStrength.MODERATE
    if n_affected >= 2:
        return CriterionStrength.SUPPORTING
    return None


class PS4Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # 1. Rarity gate (mirrors PM2): a common variant cannot satisfy PS4
        #    because case enrichment over controls is the entire point. Prefer
        #    FAF95_popmax (Karczewski 2020) over raw popmax/AF — FAF gives a
        #    conservative upper bound that is robust to small population
        #    sample sizes. An absent gnomAD record is treated as "rare".
        gd = annotation.gnomad
        if gd is None:
            rare = True
        else:
            # A failed gnomAD QC filter means AF is unreliable — refuse to
            # commit either way rather than risk a false PS4 trigger.
            if not gd.filter_pass:
                return CriteriaResult.not_met(ACMGCriterion.PS4, "gnomAD filter failed")
            faf = gd.faf95_popmax
            if faf is None:
                # Fall back to popmax_af, then AF — FAF is missing for very
                # rare variants where the upper bound is undefined.
                faf = gd.popmax_af or gd.af or 0.0
            rare = (faf == 0.0 or gd.ac == 0 or faf < _FAF95_RARE)

        if not rare:
            return CriteriaResult.not_met(
                ACMGCriterion.PS4, "Not rare enough in gnomAD for PS4"
            )

        # 2. Proxy proband count: number of ClinVar SCVs that BOTH classified
        #    the variant as P/LP AND reported AffectedStatus="yes". This is an
        #    approximation — one SCV ≈ one proband — chosen because true
        #    case/control statistics are unavailable from public data.
        from acmg_classifier.local_db.clinvar_sqlite import query_affected_cases
        n_affected = query_affected_cases(
            self._cfg.clinvar_sqlite,
            variant.chrom, variant.pos, variant.ref, variant.alt,
        )

        strength = _ps4_strength(n_affected)
        if strength is None:
            return CriteriaResult.not_met(
                ACMGCriterion.PS4,
                f"Insufficient affected probands in ClinVar ({n_affected} < 2)",
            )

        return CriteriaResult.met(
            ACMGCriterion.PS4,
            strength,
            evidence=f"Rare + {n_affected} unrelated affected SCV(s) ({strength.value})",
        )
