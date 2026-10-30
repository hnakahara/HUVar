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
    # MANE Plus Clinical: an additional transcript MANE designates for clinical
    # reporting of variants not captured by MANE Select (e.g. KRAS exon 4a,
    # MECP2 e2). Reported alongside the MANE Select evaluation; criteria are
    # still evaluated on ``primary_consequence``.
    is_mane_plus_clinical: bool = False

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


def _base_acc(acc: Optional[str]) -> str:
    """Version-stripped accession (NM_033360.4 -> NM_033360)."""
    return (acc or "").split(".")[0]


class ManePlusClinicalInfo(BaseModel):
    """Supplementary MANE Plus Clinical description of a variant.

    Reporting only: criteria are evaluated on the MANE Select (primary)
    transcript. ``consequence_differs`` flags variants whose consequence on the
    Plus Clinical transcript differs from the evaluated one (e.g. intronic on
    MANE Select but missense on the Plus Clinical exon)."""

    transcript_id: str
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    consequence: str = ""
    consequence_differs: bool = False
    # Classification with the criteria evaluated on this transcript (filled by
    # the pipeline; None when the Plus Clinical evaluation was not run).
    classification_2015: Optional[str] = None
    rules_2015: Optional[str] = None
    bayesian_score: Optional[int] = None
    classification_bayesian: Optional[str] = None
    # Criteria whose outcome differs from the MANE Select evaluation,
    # e.g. "PP3:NotMet->Moderate; PM1:Moderate->NotMet".
    criteria_changes: str = ""


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
    # Non-cancer subset allele frequency (gnomAD ``*_non_cancer`` INFO fields).
    # Used by PM2 for VCEPs that judge absence on the non-cancer subset (ENIGMA
    # BRCA1/2). None when the gnomAD DB predates this column → PM2 gracefully
    # falls back to the overall allele frequency.
    af_non_cancer: Optional[float] = None
    # Non-cancer-subset popmax FAF95 (95%-CI filtering AF), recomputed from the
    # gnomAD v3.1.2 per-group non-cancer AC/AN (gnomAD does not precompute it for
    # the subset). Used by BA1/BS1 for VCEPs that judge frequency on the
    # non-cancer subset (ENIGMA BRCA1/2). None when the DB predates this column →
    # BA1/BS1 gracefully fall back to the overall (full-release) FAF95.
    faf95_non_cancer: Optional[float] = None
    # True when the non-cancer companion DB was present and consulted for this
    # variant — so a None af_non_cancer means genuine ABSENCE in the non-cancer
    # subset (no record found), not "subset data unavailable". PM2 (ENIGMA
    # BRCA1/2) needs this to tell "present only in cancer cohorts → absent in
    # non-cancer → PM2 met" apart from an old DB with no non-cancer support.
    non_cancer_queried: bool = False
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
    # How the row was chosen among REVEL's per-transcript rows: "transcript"
    # (MANE ENST), "amino_acid" (MANE amino-acid change), "max" (no MANE match,
    # highest of several), "single" (one row) or "legacy" (5-column file).
    match: Optional[str] = None


class BayesDelData(BaseModel):
    """BayesDel meta-predictor score (Feng 2017) for a missense variant.

    Higher score ⇒ more pathogenic. Opt-in auxiliary predictor used only by the
    VCEPs that name it (e.g. ENIGMA BRCA1/2, TP53) and only under insilico_tool
    REVEL/ALPHAMISSENSE (never ESM1B — see Config.use_bayesdel)."""

    score: Optional[float] = None


class CADDData(BaseModel):
    """CADD score (Kircher 2014, PHRED-scaled) for a variant.

    Higher PHRED ⇒ more deleterious. Opt-in auxiliary predictor used only in the
    VCEP 2-of-3 combination rules and only under insilico_tool REVEL/ALPHAMISSENSE
    (never ESM1B — see Config.use_cadd)."""

    phred: Optional[float] = None


class ESM1bData(BaseModel):
    """ESM1b log-likelihood ratio (LLR) score for a missense variant.

    LLR convention (Brandes et al., Nat Genet 2023): more negative LLR ⇒ more
    pathogenic; more positive LLR ⇒ more benign. Thresholds follow
    Bergquist et al. 2024 Table 2.
    """

    llr: Optional[float] = None


class FuncVEPData(BaseModel):
    """FuncVEP-CTI score (Kayaalp et al., Nat Genet 2026) for a missense variant.

    CTI (clinical-threshold index) is on a 0-1 scale; higher ⇒ more pathogenic.
    PP3/BP4 strengths follow the paper's supplementary calibration table, which
    gives a score cutoff per ACMG point (+1..+4 / -1..-4). ``binary`` is the
    paper's single Damaging/Neutral cutoff (0.4196), reported for reference and
    not used to award evidence.

    Licence: PolyForm Strict 1.0.0 (non-commercial). The precomputed scores are
    read only when insilico_tool is FUNCVEP, so a commercial-use run never
    touches them.
    """

    score: Optional[float] = None
    # Gene symbol the chosen row belongs to (FuncVEP rows are per gene, so a
    # locus in overlapping genes carries several scores).
    gene: Optional[str] = None
    # How the row was chosen among the rows at this locus: "gene" (matched the
    # evaluated gene symbol), "single" (one row), "uniform" (several rows with
    # the same score) or "ambiguous" (conflicting scores, no gene match → no score).
    match: Optional[str] = None


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
    # Genomic coordinates of the comparator (when the query returns them); used
    # to match comparators against VCEP classifications in the eRepo.
    chrom: Optional[str] = None
    pos: Optional[int] = None
    ref: Optional[str] = None
    alt: Optional[str] = None


class RepeatMaskerRegion(BaseModel):
    """Whether the variant falls inside a repeat element."""

    in_repeat: bool = False
    repeat_class: Optional[str] = None
    repeat_name: Optional[str] = None


class TranscriptInSilico(BaseModel):
    """In-silico scores tied to one transcript/isoform (see AnnotationData.for_transcript)."""

    esm1b: Optional[ESM1bData] = None
    revel: Optional[RevelData] = None


class AnnotationData(BaseModel):
    """Aggregated annotation for a single variant from all local data sources."""

    # VEP consequences (ordered: MANE Select first, then canonical, then rest)
    consequences: list[ConsequenceInfo] = []

    gnomad: Optional[GnomADData] = None
    alphamissense: Optional[AlphaMissenseData] = None
    esm1b: Optional[ESM1bData] = None
    revel: Optional[RevelData] = None
    # FuncVEP-CTI. Populated only under insilico_tool FUNCVEP (PolyForm Strict,
    # non-commercial — see FuncVEPData), None otherwise.
    funcvep: Optional[FuncVEPData] = None
    # Opt-in auxiliary predictors, populated only under insilico_tool
    # REVEL/ALPHAMISSENSE when --with-bayesdel / --with-cadd are set (None
    # otherwise, including always under ESM1B — see AnnotationOrchestrator).
    bayesdel: Optional[BayesDelData] = None
    cadd: Optional[CADDData] = None
    # Primary splice predictor used for ACMG criteria (SpliceAI or SQUIRLS per config).
    splice: Optional[SpliceScore] = None
    # SQUIRLS score always stored separately for TSV reporting regardless of which
    # splice tool is configured as primary.
    squirls: Optional[SpliceScore] = None
    clinvar_vcf: list[ClinVarRecord] = []      # for PP5
    clinvar_sqlite: list[ClinVarRecord] = []   # for PS1 / PM5 (same-AA or same-codon)
    repeat: Optional[RepeatMaskerRegion] = None

    # Transcript the criteria are evaluated on. None = the default pick
    # (MANE Select first; see _default_primary). Set through for_transcript()
    # to evaluate on a MANE Plus Clinical transcript.
    evaluation_transcript: Optional[str] = None
    # Transcript-specific in-silico scores (ESM1b isoform / REVEL transcript row),
    # keyed by version-stripped accession. Filled by the annotation orchestrator
    # for the default transcript and each MANE Plus Clinical transcript so that
    # for_transcript() can swap in the scores of the evaluated isoform.
    transcript_insilico: dict[str, TranscriptInSilico] = {}

    def _default_primary(self) -> Optional[ConsequenceInfo]:
        """Default evaluation transcript.

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

    def _find_transcript(self, tx: str) -> Optional[ConsequenceInfo]:
        """Consequence on ``tx`` (exact id first, then version-stripped; RefSeq
        preferred among base matches)."""
        for c in self.consequences:
            if c.transcript_id == tx:
                return c
        base = _base_acc(tx)
        hits = [c for c in self.consequences if _base_acc(c.transcript_id) == base]
        hits.sort(key=lambda c: not c.transcript_id.startswith("NM_"))
        return hits[0] if hits else None

    @property
    def primary_consequence(self) -> Optional[ConsequenceInfo]:
        """The single transcript that ACMG criteria are evaluated against:
        :attr:`evaluation_transcript` when set and present, otherwise the
        default pick (MANE Select first; see :meth:`_default_primary`)."""
        if self.evaluation_transcript:
            c = self._find_transcript(self.evaluation_transcript)
            if c is not None:
                return c
        return self._default_primary()

    @property
    def rule_consequence(self) -> Optional[ConsequenceInfo]:
        """Consequence whose numbering VCEP position rules and ClinVar codon
        comparisons are read in.

        VCEP tables (PM1 hotspots, PM4 regions, BA1/BS1/BP4 variant lists,
        PVS1 exon/codon bands, PP3/BP4 domains) and ClinVar same-codon records
        are numbered on MANE Select. When the criteria are evaluated on another
        transcript (MANE Plus Clinical), those lookups still use the MANE Select
        consequence of the same gene — the same genomic residue — so a variant
        in a region shared by both transcripts gets the same rule; a variant in
        a Plus Clinical-only region has no MANE Select position and the rule is
        not met. Without an evaluation transcript this is the primary itself.
        """
        pc = self.primary_consequence
        if not self.evaluation_transcript:
            return pc
        d = self._default_primary()
        if d is not None and pc is not None and d.gene_symbol == pc.gene_symbol:
            return d
        return pc

    def for_transcript(self, transcript: Optional[str]) -> "AnnotationData":
        """Copy evaluated on ``transcript`` (None restores the default pick).

        In-silico scores are swapped for the isoform-specific ones recorded in
        :attr:`transcript_insilico` when available."""
        update: dict = {"evaluation_transcript": transcript or None}
        key = _base_acc(transcript) if transcript else _base_acc(
            getattr(self._default_primary(), "transcript_id", None))
        scores = self.transcript_insilico.get(key) if key else None
        if scores is not None:
            update["esm1b"] = scores.esm1b
            update["revel"] = scores.revel
        return self.model_copy(update=update)

    def rule_view(self) -> "AnnotationData":
        """Copy whose primary consequence is :attr:`rule_consequence` — for
        criteria that are entirely ClinVar-codon / VCEP-position based (PM1,
        PS1, PM5)."""
        if not self.evaluation_transcript:
            return self
        return self.for_transcript(None)

    @property
    def mane_plus_clinical_consequences(self) -> list[ConsequenceInfo]:
        """MANE Plus Clinical consequences of the default transcript's gene.

        VEP ``--merged`` reports each MANE transcript twice (RefSeq NM_ and
        Ensembl ENST). One entry per transcript is returned, RefSeq preferred:
        an ENST row is kept only when no NM_ Plus Clinical row exists for the
        gene. GNAS has two Plus Clinical transcripts, hence a list.
        """
        pc = self._default_primary()
        if pc is None:
            return []
        same_gene = [c for c in self.consequences
                     if c.is_mane_plus_clinical and c.gene_symbol == pc.gene_symbol
                     and c.transcript_id != pc.transcript_id]
        refseq = [c for c in same_gene if c.transcript_id.startswith("NM_")]
        return refseq or same_gene

    def mane_plus_clinical_summary(self) -> list[ManePlusClinicalInfo]:
        """Reporting view of :attr:`mane_plus_clinical_consequences`."""
        pc = self._default_primary()
        out = []
        for c in self.mane_plus_clinical_consequences:
            out.append(ManePlusClinicalInfo(
                transcript_id=c.transcript_id,
                hgvs_c=c.hgvs_c,
                hgvs_p=c.hgvs_p,
                consequence=c.consequence.value,
                consequence_differs=bool(pc and c.consequence != pc.consequence),
            ))
        return out
