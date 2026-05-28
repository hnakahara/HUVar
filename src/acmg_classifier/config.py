from __future__ import annotations
from pathlib import Path
from typing import Optional
from pydantic import field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

from acmg_classifier.models.enums import Assembly, InSilicoTool, SpliceTool


class Config(BaseSettings):
    """Runtime configuration — populated from CLI flags or environment variables."""

    model_config = SettingsConfigDict(env_prefix="ACMG_", env_file=".env", extra="ignore")

    data_dir: Path = Path("./data")
    assembly: Assembly = Assembly.GRCH38
    workers: int = 4
    vep_batch_size: int = 500
    insilico_tool: InSilicoTool = InSilicoTool.ALPHAMISSENSE
    splice_tool: SpliceTool = SpliceTool.SQUIRLS
    spliceai_dir: Optional[Path] = None

    @field_validator("data_dir")
    @classmethod
    def resolve_data_dir(cls, v: Path) -> Path:
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
    def squirls_db_dir(self) -> Path:
        names = {
            Assembly.GRCH38: "squirls-2309-hg38",
            Assembly.GRCH37: "squirls-2309-hg19",
        }
        return self.assembly_dir / "squirls" / names[self.assembly]

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
