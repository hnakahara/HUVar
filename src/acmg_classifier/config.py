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
    # When True, genes whose VCEP defines the BA1/BS1 cutoff on the *point*
    # grpmax/popmax allele frequency (af_basis="popmax" in disease_prevalence.tsv
    # — e.g. RUNX1, GAA, MYOC) compare against the gnomAD popmax point estimate
    # rather than the more conservative FAF95. Set ACMG_POPMAX_AF_BASIS=false to
    # force every gene back to FAF95 (the universally conservative metric).
    popmax_af_basis: bool = True
    insilico_tool: InSilicoTool = InSilicoTool.ESM1B
    # Opt-in auxiliary missense predictors (BayesDel, CADD) consulted only by the
    # VCEPs that name them (e.g. ENIGMA BRCA1/2 / TP53 BayesDel; CADD in 2-of-3
    # combos). Both carry academic / non-commercial licence terms, so they are
    # gated to the licence-encumbered missense path: they run ONLY when
    # insilico_tool is REVEL or ALPHAMISSENSE and are skipped entirely under
    # ESM1B (the commercial-safe default). Enable per run with --with-bayesdel /
    # --with-cadd; the data must also have been staged with the matching
    # setup_data.py flag (--with-bayesdel / --with-cadd).
    use_bayesdel: bool = False
    use_cadd: bool = False
    splice_tool: SpliceTool = SpliceTool.OPENSPLICEAI
    supplement_mode: SupplementMode = SupplementMode.MERGE
    spliceai_dir: Optional[Path] = None
    openspliceai_model_dir: Optional[Path] = None
    # Optional override of the disease_prevalence.tsv path. Left None for CLI /
    # batch runs, which always use the conservative aggregated table under
    # data_dir. The HUVar app sets this to a per-CSpec overlay TSV (a copy of the
    # base table with one gene's row replaced by a specific disease's spec) to
    # evaluate a variant under that CSpec without touching the batch behaviour.
    disease_prevalence_tsv_override: Optional[Path] = None
    # OpenSpliceAI model context length (nt). Default 80 (the smallest OSAI_MANE
    # variant) for fast CPU inference — runtime scales with the flanking window,
    # so 80nt is ~25x lighter than 2000nt. Raise to 400/2000/10000 (all staged
    # by setup_data.py) for higher splice-prediction sensitivity on GPU or when
    # throughput is not a concern. Must match a downloaded model subdirectory.
    openspliceai_flanking_size: int = 80
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
    def gnomad_noncancer_duckdb(self) -> Optional[Path]:
        """Companion non-cancer-subset AF DB for PM2 (ENIGMA BRCA1/2), GRCh38 only.

        gnomAD v4.1 dropped the non-cancer subset, so the main v4.1 JOINT build
        carries af_non_cancer = NULL everywhere. This points to a small DB built
        from gnomAD v3.1.2 genomes (hg38 — same coordinates as v4.1) holding only
        (chrom, pos, ref, alt, af_non_cancer); gnomad_db consults it as a fallback
        when the main DB's af_non_cancer is NULL. v2.1.1 (GRCh37) already carries
        the non-cancer subset inline, so no companion is needed there. Returns
        None when the file is absent (PM2 then falls back to the overall AF)."""
        if self.assembly != Assembly.GRCH38:
            return None
        p = self.assembly_dir / "gnomad" / "gnomad_v3.1.2_non_cancer.duckdb"
        return p if p.exists() else None

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
    def revel_tsv(self) -> Path:
        """REVEL scores, converted to a 5-column tabix-indexed TSV
        (chrom, pos, ref, alt, REVEL) by scripts/setup_data.py.

        REVEL is distributed as a single file carrying both hg19 and GRCh38
        coordinates; setup_data builds one per-assembly TSV indexed on the
        matching position column, so the query side stays assembly-agnostic.
        """
        names = {
            Assembly.GRCH38: "revel_grch38.tsv.gz",
            Assembly.GRCH37: "revel_grch37.tsv.gz",
        }
        return self.assembly_dir / "revel" / names[self.assembly]

    @property
    def bayesdel_tsv(self) -> Path:
        """BayesDel precomputed scores, normalised by scripts/setup_data.py to a
        5-column tabix-indexed TSV (chrom, pos, ref, alt, BayesDel).

        Opt-in (setup_data.py --with-bayesdel) and consulted only under
        insilico_tool REVEL/ALPHAMISSENSE — never under ESM1B (licence gate).
        """
        names = {
            Assembly.GRCH38: "bayesdel_grch38.tsv.gz",
            Assembly.GRCH37: "bayesdel_grch37.tsv.gz",
        }
        return self.assembly_dir / "bayesdel" / names[self.assembly]

    @property
    def cadd_tsv(self) -> Path:
        """CADD precomputed scores (PHRED), normalised by scripts/setup_data.py
        to a 5-column tabix-indexed TSV (chrom, pos, ref, alt, CADD_PHRED).

        Opt-in (setup_data.py --with-cadd) and consulted only under insilico_tool
        REVEL/ALPHAMISSENSE — never under ESM1B (licence gate).
        """
        names = {
            Assembly.GRCH38: "cadd_grch38.tsv.gz",
            Assembly.GRCH37: "cadd_grch37.tsv.gz",
        }
        return self.assembly_dir / "cadd" / names[self.assembly]

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
    def openspliceai_annotation(self) -> Path:
        """Gene annotation table for `openspliceai variant -A`.

        The CLI's `-A grch38` / `grch37` keywords resolve to a relative path
        inside the openspliceai source tree and do NOT work from an installed
        package, so we pass an explicit file. setup_data.py downloads
        grch38.txt / grch37.txt (from the OpenSpliceAI repo) to this path."""
        name = {Assembly.GRCH38: "grch38.txt", Assembly.GRCH37: "grch37.txt"}[self.assembly]
        return self.assembly_dir / "openspliceai" / name

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
        # An explicit override (per-CSpec overlay set by the app) wins; otherwise
        # the conservative aggregated table under data_dir (CLI / batch default).
        if self.disease_prevalence_tsv_override is not None:
            return self.disease_prevalence_tsv_override
        return self.data_dir / "shared" / "disease_prevalence.tsv"

    @property
    def gene_inheritance_tsv(self) -> Path:
        """gene -> inheritance (AD/AR/XL...) map for inheritance-aware PM2 thresholds."""
        return self.data_dir / "shared" / "gene_inheritance.tsv"

    @property
    def tp53_codes_tsv(self) -> Path:
        """ClinGen TP53 VCEP precomputed per-missense PP3/BP4 codes (built from the
        VCEP supplementary table by scripts/build_tp53_codes.py). Carries the
        Align-GVGD-derived code this pipeline cannot compute itself; consulted only
        under the BayesDel licence gate (insilico_tool REVEL/AlphaMissense +
        --with-bayesdel). Absent file → TP53 auxiliary PP3/BP4 simply not applied."""
        return self.data_dir / "shared" / "tp53_pp3_bp4_codes.tsv"

    @property
    def pm1_hotspots_tsv(self) -> Path:
        """Per-gene PM1 hotspot regions/residues mined from the VCEP cspec
        summaries (see scripts/build_pm1_hotspots.py)."""
        return self.data_dir / "shared" / "pm1_hotspots.tsv"

    @property
    def gnomad_coverage_db(self) -> Path:
        """Per-locus gnomAD coverage DuckDB (mean read depth) for the PM2
        read-depth gate (ENIGMA BRCA1/2). Optional — absent file disables the
        gate (PM2 then does not enforce the depth requirement). Built by
        scripts/build_gnomad_coverage.py from the coverage summary downloaded via
        ``setup_data.py --with-gnomad-coverage``."""
        ver = "4.1" if self.assembly == Assembly.GRCH38 else "2.1.1"
        return self.assembly_dir / "gnomad" / f"gnomad_v{ver}_exomes_coverage.duckdb"

    @property
    def ps1_paralog_map_tsv(self) -> Path:
        """SCN paralogue amino-acid alignment for the PS1 analogous-residue route
        (SCN1A/2A/3A/8A), built by scripts/build_ps1_paralog_map.py. Absent file
        disables the SCN paralogue path."""
        return self.data_dir / "shared" / "ps1_paralog_map.tsv"

    @property
    def pm4_regions_tsv(self) -> Path:
        """Per-gene PM4 region / strength rules (Strong residues, allow/deny
        regions, region default, stop-loss strength) — see
        scripts/build_pm4_regions.py. Absent file → flat default PM4 behaviour."""
        return self.data_dir / "shared" / "pm4_regions.tsv"

    @property
    def vcep_pvs1_splice_exons_tsv(self) -> Path:
        """Optional reviewer-supplied per-(gene, skipped-exon) PVS1 splice-strength
        overrides. Absent file → flat per-gene splice defaults (unchanged
        behaviour). See scripts/build_vcep_pvs1_exons.py and pvs1/vcep_pvs1_exons.py."""
        return self.data_dir / "shared" / "vcep_pvs1_splice_exons.tsv"
