from __future__ import annotations
from pathlib import Path
from typing import Optional
from pydantic import field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

from acmg_classifier.models.enums import Assembly, InSilicoTool, SpliceTool, SupplementMode


class Config(BaseSettings):
    """Runtime configuration — populated from CLI flags or environment variables.

    Backed by pydantic-settings so users can override any field via
    ACMG_* environment variables or a .env file. `extra="ignore"` lets
    callers pass extra kwargs (e.g. from a CLI option that maps to no
    Config field) without raising, which is essential for the
    forward-compatible CLI flags."""

    model_config = SettingsConfigDict(env_prefix="ACMG_", env_file=".env", extra="ignore")

    data_dir: Path = Path("./data")
    assembly: Assembly = Assembly.GRCH38
    workers: int = 4
    vep_batch_size: int = 500
    # Minimum ClinVar review stars for a PM5 same-codon comparator. Default 1
    # ("criteria provided"; excludes only 0-star no-assertion records). A
    # single-submitter (1-star) P/LP is often the only — and a legitimate —
    # PM5 anchor (e.g. comparators the eRepo truth set relies on), so requiring
    # 2 stars drops real PM5 calls. Precision is instead handled by the surgical
    # PM5 gates (benign-at-codon, Grantham, LP->Supporting, PM1 exclusion). Set
    # ACMG_PM5_MIN_STARS=2 for a stricter expert/multi-submitter-only policy.
    pm5_min_stars: int = 1
    # BS2 (observed in a healthy individual) gnomAD count thresholds, by the
    # gene's inheritance mode. Recessive counts homozygotes, X-linked counts
    # hemizygotes, dominant counts heterozygous carriers (AC - nhomalt). Defaults
    # follow the modal ClinGen VCEP BS2 counts (most fire on 1-2 homozygotes /
    # ~3 unaffected heterozygotes — "homozygous in a healthy adult"); the prior
    # flat 5 was too strict and missed benign evidence. Genes whose VCEP states
    # a higher bar (cancer panels: CDH1 >=10, TP53 >=8) carry a per-gene
    # `bs2_count` in disease_prevalence.tsv that overrides these. Override the
    # global defaults via ACMG_BS2_MIN_HOMALT / ACMG_BS2_MIN_HEMI / ACMG_BS2_MIN_HET.
    bs2_min_homalt: int = 2
    bs2_min_hemi: int = 2
    bs2_min_het: int = 3
    insilico_tool: InSilicoTool = InSilicoTool.ESM1B
    splice_tool: SpliceTool = SpliceTool.OPENSPLICEAI
    supplement_mode: SupplementMode = SupplementMode.MERGE
    spliceai_dir: Optional[Path] = None
    openspliceai_model_dir: Optional[Path] = None
    openspliceai_flanking_size: int = 2000
    # BP7 conservation gate: a synonymous/deep-intronic variant may only reach
    # BP7 when the nucleotide is NOT highly conserved (Walker 2023 / ACMG 2015).
    # phyloP100way is the conservation source (UCSC, commercial-use OK). A
    # position with phyloP >= this value is "highly conserved" and blocks BP7.
    # The gate is applied only when the phyloP file is present (graceful
    # degradation); otherwise BP7 falls back to the splice-only logic.
    bp7_phylop_max: float = 2.0

    @field_validator("data_dir")
    @classmethod
    def resolve_data_dir(cls, v: Path) -> Path:
        # Resolve to an absolute path on construction so that subprocess
        # invocations (VEP) and tabix lookups see a stable, fully-qualified
        # location regardless of the caller's working directory.
        return v.resolve()

    # ---- derived paths ----

    @property
    def assembly_dir(self) -> Path:
        return self.data_dir / self.assembly.value

    @property
    def vep_cache_dir(self) -> Path:
        return self.data_dir / "vep_cache"

    @property
    def genome_fasta(self) -> Path:
        names = {
            Assembly.GRCH38: "GRCh38.p14.fa",
            Assembly.GRCH37: "GRCh37.p13.fa",
        }
        return self.assembly_dir / "genome" / names[self.assembly]

    @property
    def gnomad_duckdb(self) -> Path:
        # GRCh38 uses the gnomAD v4.1 JOINT (combined exome+genome) build; v2.1.1
        # has no joint release, so exomes + genomes are both loaded and merged at
        # query time (per-field MAX — see gnomad_db._merge_rows).
        names = {
            Assembly.GRCH38: "gnomad_v4.1_joint.duckdb",
            Assembly.GRCH37: "gnomad_v2.1.1_exome_genome.duckdb",
        }
        return self.assembly_dir / "gnomad" / names[self.assembly]

    @property
    def gnomad_constraint_tsv(self) -> Path:
        names = {
            Assembly.GRCH38: "gnomad_v4.1_constraint.tsv",
            Assembly.GRCH37: "gnomad_v2.1.1_constraint.tsv",
        }
        return self.assembly_dir / "gnomad" / names[self.assembly]

    @property
    def clinvar_vcf(self) -> Path:
        return self.assembly_dir / "clinvar" / f"clinvar_{self.assembly.value}.vcf.gz"

    @property
    def clinvar_sqlite(self) -> Path:
        return self.assembly_dir / "clinvar" / f"clinvar_ps1_pm5_{self.assembly.value}.sqlite"

    @property
    def alphamissense_tsv(self) -> Path:
        names = {
            Assembly.GRCH38: "AlphaMissense_hg38.tsv.gz",
            Assembly.GRCH37: "AlphaMissense_hg19.tsv.gz",
        }
        return self.assembly_dir / "alphamissense" / names[self.assembly]

    @property
    def esm1b_sqlite(self) -> Path:
        """Brandes 2023 ESM1b LLR scores indexed by UniProt accession.

        Built from `ALL_hum_isoforms_ESM1b_LLR.zip` (see scripts/setup_data.py).
        Stored under `data_dir/esm1b/` rather than `assembly_dir/` because
        ESM1b scores are protein-coordinate and reused across assemblies.
        """
        return self.data_dir / "esm1b" / "esm1b_llr.sqlite"

    @property
    def squirls_db_dir(self) -> Path:
        # Retained for backward compatibility; SQUIRLS is no longer selectable
        # (its precomputed DB is no longer downloadable). See SpliceTool enum.
        names = {
            Assembly.GRCH38: "2203_hg38",
            Assembly.GRCH37: "2203_hg19",
        }
        return self.assembly_dir / "squirls" / names[self.assembly]

    # MMSplice GTF path — DISABLED along with the rest of the MMSplice
    # integration (dependency conflict). Retained, commented out, for later use.
    # @property
    # def mmsplice_gtf(self) -> Path:
    #     """Ensembl gene annotation GTF used by MMSplice (runtime splice scoring).
    #
    #     Pre-filtered to protein-coding genes by scripts/setup_data.py. The
    #     chromosome naming must match genome_fasta so MMSplice's dataloader can
    #     join annotation to reference sequence.
    #     """
    #     names = {
    #         Assembly.GRCH38: "Homo_sapiens.GRCh38.111.protein_coding.gtf",
    #         Assembly.GRCH37: "Homo_sapiens.GRCh37.87.protein_coding.gtf",
    #     }
    #     return self.assembly_dir / "mmsplice" / names[self.assembly]

    @property
    def _spliceai_base_dir(self) -> Path:
        return self.spliceai_dir.resolve() if self.spliceai_dir else self.assembly_dir / "spliceai"

    @property
    def spliceai_vcf(self) -> Optional[Path]:
        names = {
            Assembly.GRCH38: "spliceai_scores.raw.snv.hg38.vcf.gz",
            Assembly.GRCH37: "spliceai_scores.raw.snv.hg19.vcf.gz",
        }
        p = self._spliceai_base_dir / names[self.assembly]
        return p if p.exists() else None

    @property
    def spliceai_indel_vcf(self) -> Optional[Path]:
        names = {
            Assembly.GRCH38: "spliceai_scores.raw.indel.hg38.vcf.gz",
            Assembly.GRCH37: "spliceai_scores.raw.indel.hg19.vcf.gz",
        }
        p = self._spliceai_base_dir / names[self.assembly]
        return p if p.exists() else None

    @property
    def openspliceai_model_path(self) -> Path:
        """Path to the OpenSpliceAI model directory.

        If --openspliceai-model-dir is given, that path is used as-is.
        Otherwise defaults to data/<assembly>/openspliceai/<flanking_size>nt/,
        which matches the OSAI_MANE directory layout from the openspliceai repo.
        """
        if self.openspliceai_model_dir:
            return self.openspliceai_model_dir.resolve()
        return self.assembly_dir / "openspliceai" / f"{self.openspliceai_flanking_size}nt"

    @property
    def phylop_bigwig(self) -> Optional[Path]:
        """phyloP100way conservation bigWig for the BP7 conservation gate, or
        None when not downloaded (BP7 then uses its splice-only logic).

        UCSC tracks: hg38.phyloP100way.bw / hg19.phyloP100way.bw."""
        names = {
            Assembly.GRCH38: "hg38.phyloP100way.bw",
            Assembly.GRCH37: "hg19.phyloP100way.bw",
        }
        p = self.assembly_dir / "conservation" / names[self.assembly]
        return p if p.exists() else None

    @property
    def repeatmasker_bed(self) -> Path:
        names = {
            Assembly.GRCH38: "repeatmasker_dfam_hg38.bed.gz",
            Assembly.GRCH37: "repeatmasker_dfam_hg19.bed.gz",
        }
        return self.assembly_dir / "repeats" / names[self.assembly]

    @property
    def disease_prevalence_tsv(self) -> Path:
        return self.data_dir / "shared" / "disease_prevalence.tsv"

    @property
    def gene_inheritance_tsv(self) -> Path:
        """gene -> inheritance (AD/AR/XL...) map for inheritance-aware PM2 thresholds."""
        return self.data_dir / "shared" / "gene_inheritance.tsv"

    @property
    def pm1_hotspots_tsv(self) -> Path:
        """Per-gene PM1 hotspot regions/residues mined from the VCEP cspec
        summaries (see scripts/build_pm1_hotspots.py)."""
        return self.data_dir / "shared" / "pm1_hotspots.tsv"
