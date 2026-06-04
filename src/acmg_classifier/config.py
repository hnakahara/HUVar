from __future__ import annotations
from pathlib import Path
from typing import Optional
from pydantic import field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

from acmg_classifier.models.enums import Assembly, InSilicoTool, SpliceTool


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
    # Minimum ClinVar review stars for a PM5 same-codon comparator. Default 2
    # (multiple submitters / expert panel) excludes single-submitter (1-star)
    # P/LP assertions, which inflate PM5 over-assignment. This relies on the
    # ClinVar SQLite being built from the CURRENT RCV_release XML: the legacy
    # RCV_xml_old_format build was frozen at 2025-07, so VCEP expert-panel
    # re-reviews (e.g. PIK3CD Tyr524Asn -> 3-star) were missing and 2 stars
    # dropped real PM5 calls. With the up-to-date build those comparators are
    # 3-star and pass. Set ACMG_PM5_MIN_STARS=1 to also admit single-submitter
    # comparators (higher recall, lower precision).
    pm5_min_stars: int = 2
    insilico_tool: InSilicoTool = InSilicoTool.ALPHAMISSENSE
    splice_tool: SpliceTool = SpliceTool.NONE
    spliceai_dir: Optional[Path] = None

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
        names = {
            Assembly.GRCH38: "gnomad_v4.1_exomes.duckdb",
            Assembly.GRCH37: "gnomad_v2.1.1_exomes.duckdb",
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
