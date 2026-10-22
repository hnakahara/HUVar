# HUVar benchmark scripts

Scripts used for the analyses in the HUVar manuscript (Genome Medicine, revised version):
the ClinGen Evidence Repository (eRepo) benchmark, the agreement analyses with the eRepo
classifications, the rule calibrations (LoF mechanism for PVS1, PM1, PP2), the comparison
of disease-specific VCEP specifications, the resource measurements and the figures.

The input VCFs of the eRepo benchmark are distributed in [`inputs/`](../inputs):

| File | Variants | Description |
|---|---|---|
| `inputs/erepo_hg38_may_2026.vcf` | 11,545 records | eRepo variants (snapshot of 28 May 2026) converted to GRCh38 with TransVar; `ID`/`INFO NAME` = `GENE:HGVSc` |
| `inputs/erepo_hg19_may_2026.vcf` | 11,574 records | the same on GRCh37 |

Records that TransVar converted incorrectly were corrected from the genomic HGVS expression
of the ClinGen Allele Registry provided with each eRepo record (`prepare/vcffix_hg38.tsv`,
`prepare/vcffix_hg19.tsv`; `prepare/fix_vcf_from_erepo_hgvs.py`).

## Requirements

- HUVar installed from this repository (`pip install -e .`) with the data set up
  (`python scripts/setup_data.py`, see the main README). The curated-evidence supplement is
  expected in `$DATA/shared/` (copy `resources/GRCh38/erepo_manual_criteria_hg38.tsv` and
  `resources/GRCh37/erepo_manual_criteria_hg19.tsv`).
- For the VCEP-oriented and AlphaMissense configurations: the precomputed SpliceAI scores
  (`SPLICEAI_DIR`), REVEL, BayesDel, CADD and AlphaMissense, obtained under their own licences.
- Python packages for the analyses: pandas, numpy, scipy, statsmodels, matplotlib, openpyxl;
  for the figure PDFs also pypdf and reportlab (and the Arial fonts, `FONT_DIR`).
- For the comparison with other tools: the outputs of InterVar (v2.2.2) and BIAS-2015 (v2.1.1)
  run on the same VCFs, and the eRepo report (`erepo-original.csv`, downloaded from
  https://erepo.genome.network/evrepo/).

## Folder structure

| Folder | Contents |
|---|---|
| `prepare/` | correction of the benchmark VCFs; construction of the curated-evidence supplement from the eRepo (`make_manual_criteria_hg38.py` / `_hg19.py`, with `pvs1_rna_list.tsv` for PVS1 based on RNA studies) |
| `run/` | HUVar runs: eRepo benchmark (`run_benchmark.sh`), agreement analyses (`run_erepo_concordance.sh`), multi-specification genes (`run_cspec_test.sh`), PS1/PM5 comparator review status (`run_stars2_test.sh`), run time and memory (`measure_time_memory.sh`), CGP cohort (`run_cgp.sh`) |
| `analysis/` | performance metrics and the analyses reported in the manuscript (see below) |
| `figures/` | figure scripts (`make_figs.sh`, `build_suppfigs.py`) |

## Workflow

```bash
export DATA=/path/to/huvar/data SPLICEAI_DIR=/path/to/spliceai
cd benchmark/run

# 1. Classify the eRepo variants with the three predictor configurations,
#    with automated evidence only (--exclude-self-expert-panel) and with the curated supplement
bash run_benchmark.sh GRCh38 GRCh37            # -> erepo_results/

# 2. Category- and criterion-level performance (InterVar and BIAS-2015 outputs in the same folder)
cp ../../data/erepo_hg*.vcf erepo-original.csv <InterVar/BIAS outputs> erepo_results/
python ../analysis/run_all.py --input erepo_results --outdir results

# 3. Agreement of the five-tier classifications with the eRepo (overall, by class and by gene)
bash run_erepo_concordance.sh erepo_results results/concordance /path/to/HUVar hg38 hg19
python ../analysis/curated_evidence_discordance.py --dir results/concordance/hg38 \
    --disease-prevalence ../../resources/shared/disease_prevalence.tsv
```

### Analyses (`analysis/`)

| Script | Analysis |
|---|---|
| `run_all.py`, `verify_parsing.py`, `analyze_benchmark.py`, `analyze_huhvar_classification.py`, `make_*_table.py`, `make_criteria_metrics.py` | category- and criterion-level sensitivity, specificity, precision and F1; criterion assignment counts (Figures 2–4) |
| `erepo_concordance.py`, `erepo_concordance_diag.py`, `erepo_confusion.py`, `erepo_discordant_table.py`, `erepo_category_counts.py`, `make_fig_gene_concordance.py` | five-tier agreement with the eRepo, by class and by gene (Table 2, Supplementary Tables S11–S12, Supplementary Figure S7) |
| `erepo_codes_reclass.py` | classification of the eRepo records from their own evidence codes (2015 rules and point system) |
| `curated_evidence_discordance.py` | contribution of curated evidence to the discordance |
| `erepo_date_analysis.py`, `erepo_date_stratified.py` | effect of the date of the VCEP evaluation (dates from the ClinVar expert-panel records) |
| `lof_mechanism_genes.py`, `lof_rule_comparison.py`, `make_lof_table.py` | LoF-mechanism rules for PVS1 (Supplementary Tables S3–S4) |
| `pm1_pp2_calibration.py`, `pm1_region_search.py` | calibration of PM1 and PP2 against the VCEP decisions (Supplementary Table S8) |
| `bs2_erepo_eval.py`, `ep_import_validate.py`, `clinvar_familydata_survey.py` | BS2 rules, import of expert-panel criteria from ClinVar, ClinVar family data |
| `compare_cspec.py` | genes with several VCEP specifications (Supplementary Table S14) |
| `cgp_vcep_comparison.py` | CGP cohort, genes with and without VCEP specifications |

## Data not included

The comprehensive genomic profiling (CGP) cohort is derived from patients and the input VCF is
not distributed; the variants and their classifications are given in Supplementary Table S13
of the article. The NCI TP53 Database functional data, the eRepo report and the outputs of
the third-party tools must be obtained from their providers.
