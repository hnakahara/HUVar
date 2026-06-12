from __future__ import annotations
from typing import Optional
from pydantic import BaseModel

from acmg_classifier.models.enums import ConsequenceType


class ConsequenceInfo(BaseModel):
    """Single VEP consequence for one transcript."""

    transcript_id: str
    gene_id: str
    gene_symbol: str
    consequence: ConsequenceType
    biotype: str
    is_canonical: bool = False
    is_mane_select: bool = False

    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    hgvs_g: Optional[str] = None
    exon: Optional[str] = None        # e.g. "5/27"
    intron: Optional[str] = None      # e.g. "4/26"
    domains: list[str] = []

    # Protein-level
    amino_acids: Optional[str] = None  # e.g. "R/H"
    codons: Optional[str] = None
    protein_position: Optional[int] = None
    codon_position: Optional[int] = None  # genomic codon start (for PS1/PM5)
    # SwissProt (preferred) / TrEMBL accession from VEP --uniprot; used as the
    # join key against the Brandes 2023 ESM1b LLR archive.
    uniprot_id: Optional[str] = None

    # Splice distance for BP7
    intron_distance_from_splice: Optional[int] = None


class GnomADData(BaseModel):
    """Frequency and constraint data from gnomAD."""

    af: Optional[float] = None
    an: Optional[int] = None
    ac: Optional[int] = None
    nhomalt: Optional[int] = None
    nhemi: Optional[int] = None
    popmax_af: Optional[float] = None
    popmax_pop: Optional[str] = None
    faf95_popmax: Optional[float] = None
    # Male (XY/hemizygous) allele frequency — used by BA1/BS1 for genes whose
    # VCEP defines the cutoff "in males" (X-linked: RPGR, RS1, ABCD1, SLC6A8,
    # OTC). None when the gnomAD DB predates the AF_XY column (graceful
    # fallback to the overall FAF in those evaluators).
    af_xy: Optional[float] = None
    # Female (XX) allele count and female homozygote count, used by BS2 for VCEPs
    # that count only females (e.g. TP53). Female carriers = ac_xx - nhomalt_xx.
    # None when the gnomAD DB predates these columns (graceful fallback: a
    # female-only BS2 gene then withholds BS2 rather than counting both sexes).
    ac_xx: Optional[int] = None
    nhomalt_xx: Optional[int] = None
    # GrpMax-population allele count / number (paired), used by the PM2
    # upper-95%-CI rule (Cardiomyopathy/HCM VCEP) to reconstruct the CI UPPER
    # bound of the highest-frequency subpopulation's AF — gnomAD only exposes the
    # FAF (CI lower bound). None when the DB predates these columns.
    ac_grpmax: Optional[int] = None
    an_grpmax: Optional[int] = None
    filter_pass: bool = True

    # Gene-level constraint (from gnomAD constraint table)
    pli: Optional[float] = None
    loeuf: Optional[float] = None
    # Missense z-score (mis.z_score in v4.1 / mis_z in v2.1.1) — used by PP2
    # alongside the ClinVar benign-missense rate to detect genes where missense
    # is a common disease mechanism.
    mis_z: Optional[float] = None


class AlphaMissenseData(BaseModel):
    """AlphaMissense pathogenicity score for a missense variant."""

    score: Optional[float] = None
    classification: Optional[str] = None  # "likely_pathogenic" / "ambiguous" / "likely_benign"


class RevelData(BaseModel):
    """REVEL ensemble pathogenicity score (Ioannidis et al. 2016) for a
    missense variant. Higher score ⇒ more pathogenic (0–1 scale)."""

    score: Optional[float] = None


class ESM1bData(BaseModel):
    """ESM1b log-likelihood ratio (LLR) score for a missense variant.

    LLR convention (Brandes et al., Nat Genet 2023): more negative LLR ⇒ more
    pathogenic; more positive LLR ⇒ more benign. Thresholds follow
    Bergquist et al. 2024 Table 2.
    """

    llr: Optional[float] = None


class SpliceScore(BaseModel):
    """Unified splice prediction result (SQUIRLS or SpliceAI)."""

    tool: str                              # "squirls" or "spliceai"
    max_delta: Optional[float] = None      # SpliceAI max delta score equivalent
    raw_score: Optional[float] = None      # tool-native score
    is_available: bool = False


class ClinVarRecord(BaseModel):
    """ClinVar classification for a specific variant or amino-acid position."""

    variation_id: Optional[str] = None
    clinical_significance: str
    review_status: str
    star_rating: int = 0
    gene_symbol: Optional[str] = None
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    amino_acid_change: Optional[str] = None


class RepeatMaskerRegion(BaseModel):
    """Whether the variant falls inside a repeat element."""

    in_repeat: bool = False
    repeat_class: Optional[str] = None
    repeat_name: Optional[str] = None


class AnnotationData(BaseModel):
    """Aggregated annotation for a single variant from all local data sources."""

    # VEP consequences (ordered: MANE Select first, then canonical, then rest)
    consequences: list[ConsequenceInfo] = []

    gnomad: Optional[GnomADData] = None
    alphamissense: Optional[AlphaMissenseData] = None
    esm1b: Optional[ESM1bData] = None
    revel: Optional[RevelData] = None
    # Primary splice predictor used for ACMG criteria (SpliceAI or SQUIRLS per config).
    splice: Optional[SpliceScore] = None
    # SQUIRLS score always stored separately for TSV reporting regardless of which
    # splice tool is configured as primary.
    squirls: Optional[SpliceScore] = None
    clinvar_vcf: list[ClinVarRecord] = []      # for PP5
    clinvar_sqlite: list[ClinVarRecord] = []   # for PS1 / PM5 (same-AA or same-codon)
    repeat: Optional[RepeatMaskerRegion] = None

    @property
    def primary_consequence(self) -> Optional[ConsequenceInfo]:
        """Pick the single transcript that ACMG criteria will be evaluated against.

        Priority order is RefSeq MANE Select > Ensembl MANE Select >
        RefSeq canonical > Ensembl canonical > first available. RefSeq (NM_)
        is preferred because clinical HGVS reporting and ClinVar use RefSeq
        accessions; ranking MANE Select first follows the 2023 ClinGen/HGVS
        recommendation to standardise on MANE for clinical variant reporting.
        """
        # MANE Select on a RefSeq transcript — the strongest clinical match.
        for c in self.consequences:
            if c.is_mane_select and c.transcript_id.startswith("NM_"):
                return c
        # Fall back to MANE Select on Ensembl (ENST) if RefSeq mapping is missing.
        for c in self.consequences:
            if c.is_mane_select:
                return c
        # No MANE record: use the canonical transcript, preferring RefSeq.
        for c in self.consequences:
            if c.is_canonical and c.transcript_id.startswith("NM_"):
                return c
        for c in self.consequences:
            if c.is_canonical:
                return c
        # Last resort: take whatever VEP returned first (preserves VEP's own ordering).
        return self.consequences[0] if self.consequences else None
