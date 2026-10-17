"""PM1 — located in mutational hotspot or critical functional domain."""
from __future__ import annotations

import re
from acmg_classifier.config import Config
from acmg_classifier.criteria.base import CriterionEvaluator
from acmg_classifier.criteria.pm1_hotspots import PM1_TOLERATED_SUBSTITUTIONS, PM1Hotspots
from acmg_classifier.models.annotation import AnnotationData
from acmg_classifier.models.criteria import CriteriaResult
from acmg_classifier.models.enums import ACMGCriterion, ConsequenceType, CriterionStrength
from acmg_classifier.models.variant import VariantRecord
from acmg_classifier.models.supplement import SupplementEntry

_RANK = {
    CriterionStrength.SUPPORTING: 1,
    CriterionStrength.MODERATE: 2,
    CriterionStrength.STRONG: 3,
}

# PM1 is a missense criterion by default. Only the few VCEPs that explicitly
# extend it to in-frame indels accept those: GAA (GN010, "missense substitution
# or in frame deletion"), IDUA (GN091, "missense substitutions or inframe
# deletions") and CYP1B1 (GN104, "missense variants and in-frame indels"). For
# every other gene an in-frame indel must NOT earn PM1 (e.g. RUNX1, whose VCEP
# scopes PM1 to missense affecting RHD residues).
_PM1_INDEL_GENES = frozenset({"GAA", "IDUA", "CYP1B1"})



# HNF1A (MDEP v3.1) PM1_Supporting also covers transcription-factor binding
# sites of the promoter (c. positions upstream of the start codon, inclusive).
_HNF1A_PROMOTER_SITES = (
    (-195, -187, "AP1 binding site"),
    (-227, -209, "overlapping HNF3 / NF-Y sites"),
    (-259, -238, "HNF1A binding site"),
    (-288, -276, "HNF4A binding site"),
)
_UPSTREAM_SNV = re.compile(r"c\.(-\d+)[ACGT]>[ACGT]$")


def _hnf1a_promoter_site(pc) -> str | None:
    if pc is None or pc.gene_symbol != "HNF1A" or not pc.hgvs_c:
        return None
    m = _UPSTREAM_SNV.search(pc.hgvs_c.split(":")[-1])
    if not m:
        return None
    pos = int(m.group(1))
    for lo, hi, name in _HNF1A_PROMOTER_SITES:
        if lo <= pos <= hi:
            return f"c.{pos} in the {name} (c.{lo} to c.{hi})"
    return None

class PM1Evaluator(CriterionEvaluator):
    def __init__(self, cfg: Config) -> None:
        self._cfg = cfg
        # Per-gene VCEP hotspot definitions (residue ranges / residues / strength)
        # mined from the cspec summaries. Authoritative where present.
        self._hotspots = PM1Hotspots(cfg.pm1_hotspots_tsv)

    def evaluate(
        self,
        variant: VariantRecord,
        annotation: AnnotationData,
        supplement: list[SupplementEntry] | None = None,
    ) -> CriteriaResult:
        # PM1 only applies to variants that change protein sequence at a
        # specific codon — missense and small in-frame indels. Truncating
        # variants (frameshift, stop_gain) are PVS1 territory, and silent
        # variants do not move an amino-acid position into a hotspot.
        pc = annotation.primary_consequence
        prom = _hnf1a_promoter_site(pc)
        if prom is not None:
            return CriteriaResult.met(
                ACMGCriterion.PM1, CriterionStrength.SUPPORTING,
                f"HNF1A VCEP: promoter {prom} -> PM1_Supporting",
            )
        if pc is None or pc.consequence not in (
            ConsequenceType.MISSENSE,
            ConsequenceType.INFRAME_INSERTION,
            ConsequenceType.INFRAME_DELETION,
        ):
            return CriteriaResult.not_met(ACMGCriterion.PM1, "Not a missense/in-frame variant")

        gene = pc.gene_symbol

        # PM1 is missense-only unless the gene's VCEP explicitly extends it to
        # in-frame indels — withhold PM1 for an in-frame indel in any other gene.
        if pc.consequence in (
            ConsequenceType.INFRAME_INSERTION, ConsequenceType.INFRAME_DELETION,
        ) and gene not in _PM1_INDEL_GENES:
            return CriteriaResult.not_met(
                ACMGCriterion.PM1, f"{gene}: PM1 applies to missense variants only"
            )

        # 1. VCEP declared PM1 not applicable for the gene (benign variation
        #    throughout / no defined hotspot — e.g. ABCA4, ATM, RASopathy genes).
        if self._hotspots.is_not_applicable(gene):
            return CriteriaResult.not_met(
                ACMGCriterion.PM1, f"{gene}: VCEP designates PM1 not applicable"
            )

        # 2. VCEP-curated hotspot regions are authoritative where defined: a hit
        #    awards PM1 at the VCEP strength; a miss withholds PM1 (do NOT fall
        #    back to any statistical hotspot estimate).
        if self._hotspots.has_gene(gene) or self._hotspots.has_cys_creating(gene):
            strength = self._hotspots.lookup(
                gene, pc.protein_position, _alt_aa(pc.amino_acids)
                if pc.consequence == ConsequenceType.MISSENSE else None,
            )
            evidence = f"{gene} residue {pc.protein_position} in VCEP PM1 hotspot"
            # Cys-creating missense: a substitution introducing a new cysteine in a
            # disulfide-bonded domain (FBN1 EGF/cbEGF/TB/hybrid) earns PM1_Moderate
            # regardless of the curated residue list.
            if (pc.consequence == ConsequenceType.MISSENSE
                    and _alt_aa(pc.amino_acids) == "C"
                    and self._hotspots.in_cys_creating_region(gene, pc.protein_position)):
                if strength is None or _RANK[strength] < _RANK[CriterionStrength.MODERATE]:
                    strength = CriterionStrength.MODERATE
                    evidence = (f"{gene} Cys-creating missense at residue "
                                f"{pc.protein_position} in a disulfide-bonded domain")
            if strength is None:
                return CriteriaResult.not_met(
                    ACMGCriterion.PM1, f"{gene}: not in a VCEP PM1 hotspot region"
                )
            # Substitutions the VCEP lists as tolerated at hotspot residues
            # (FBN1: N>S, G>A) do not receive PM1.
            sub = _ref_alt(pc.amino_acids) if pc.consequence == ConsequenceType.MISSENSE else None
            if sub and sub in PM1_TOLERATED_SUBSTITUTIONS.get(gene, frozenset()):
                return CriteriaResult.not_met(
                    ACMGCriterion.PM1,
                    f"{gene}: {sub[0]}>{sub[1]} at residue {pc.protein_position} "
                    "may be tolerated (VCEP caveat) — PM1 not applied",
                )
            return CriteriaResult.met(
                ACMGCriterion.PM1, strength=strength, evidence=evidence,
            )

        # 3. No curated VCEP data for the gene: PM1 only from curated evidence.
        #    The ClinVar hotspot heuristic (nearby P/LP missense clustering) was
        #    removed permanently: calibrated against the VCEP-defined regions it
        #    reached a residue-level F1 of only ~0.3 (Genome Medicine revision).
        return CriteriaResult.not_met(
            ACMGCriterion.PM1,
            f"{gene}: no VCEP PM1 definition (PM1 from curated evidence only)",
        )


def _alt_aa(amino_acids: str | None) -> str | None:
    """The 1-letter alternate residue from a VEP ``amino_acids`` "X/Y" missense
    field, or None when absent / not a single-residue substitution."""
    if not amino_acids or "/" not in amino_acids:
        return None
    alt = amino_acids.partition("/")[2].strip()
    return alt if len(alt) == 1 else None


def _ref_alt(amino_acids: str | None) -> tuple[str, str] | None:
    """``(ref, alt)`` 1-letter residues from a VEP ``amino_acids`` "X/Y" field."""
    if not amino_acids or "/" not in amino_acids:
        return None
    ref, _, alt = amino_acids.partition("/")
    ref, alt = ref.strip().upper(), alt.strip().upper()
    return (ref, alt) if len(ref) == 1 and len(alt) == 1 else None
