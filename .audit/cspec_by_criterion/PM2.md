# PM2 — Applicable specializations across populated cspec specs
Total Applicable strength entries: 132

## GN001 v1.0.0 [None] genes=
VCEP: Standards and guidelines for the interpretation of sequence variants: a joint co
Strength: **Moderate**
Desc: Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes or Exome Aggregation Consortium. Caveat: Population data for indels may be poorly called by next generation sequencing.

## GN001 v1.0.0 [None] genes=
VCEP: Standards and guidelines for the interpretation of sequence variants: a joint co
Strength: **Supporting**
Desc: 

## GN002 v2.0.0 [Released] genes=MYH7
VCEP: ClinGen Cardiomyopathy Expert Panel Specifications to the ACMG/AMP Variant Inter
Strength: **Supporting**
Desc: The values used to calculate the PM2 thresholds were derived from studies in Northern European populations that have been relatively well-characterized with regards to disease prevalence and variant spectrum. These thresholds can be applied to any population where disease prevalence is considered comparable (1/500 or lower), where the most frequent pathogenic variant accounts for no more than 2% of cases (e.g., has an allele frequency of ≤0.02 in cases based on the upper bound of 95% CI), and where the penetrance of a pathogenic variant is expected to be at least 50% (Kelly _et al._ 2018[<sup>11</sup>](#pmid_29300372)). A threshold of **≤0.00004** in the subpopulation with the highest frequency when using the upper bound of the 95% CI activates this rule. 1. Alternatively, this is equivalent to the variant NOT being observed more than once (≤1 allele) in gnomAD v.2.1.1 in one of the non-founder populations (e.g., absence required from the Other and Ashkenazi Jewish subpopulations). 2. Applying a threshold of ≤0.00004 (upper bound of 95% CI of the allele frequency in gnomAD) is equivalent to the variant being seen in a single subpopulation and that subpopulation meets any of the following: * **Allele Count (AC) in Allele Number (AN)** * ≤1 in ≥120,000 * ≤2 in ≥160,000 * ≤3 in ≥195,000 * ≤4 in ≥230,000 gnomAD is the preferred database for this calculation, but currently only displays the filtering allele frequency (FAF), which is equivalent to a lower bound estimate of the 95% CI, when the upper bound is what is needed. * Confidence interval tools, such as [Confit-de-MAF](https://www.genecalculators.net/confit-de-maf.html), can be used to determine the upper bound of the 95% CI of the observed allele frequency. Due to current technical limitations of next generation sequencing technologies, minor allele frequencies for complex variants (e.g., large indels) may not be accurately represented in population databases. Caution should be used when a variant is only identified, or over-represented, in one of the smaller gnomAD populations, as the gnomAD allele frequencies may not accurately represent the true population frequency. Population databases may contain affected or pre-symptomatic individuals for diseases with reduced penetrance/variable onset.

## GN003 v3.2.0 [Released] genes=PTEN
VCEP: ClinGen PTEN Expert Panel Specifications to the ACMG/AMP Variant Interpretation 
Strength: **Supporting**
Desc: Absent in population * Databases present at \<0.00001 (0.001%) allele frequency in gnomAD or another large sequenced population. If multiple alleles are present within any subpopulation, allele frequency in that subpopulation must be \<0.00002 (0.002%).

## GN004 v1.0.0 [Released] genes=SHOC2,NRAS,RAF1,SOS1,SOS2,PTPN11,KRAS,MAP2K1,HRAS,RIT1,MAP2K2,BRAF
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Moderate**
Desc: The variant must be completely absent from all population databases.

## GN005 v2.0.0 [Released] genes=CDH23,COCH,GJB2,KCNQ4,MYO6,MYO7A,SLC26A4,TECTA,USH2A
VCEP: ClinGen Hearing Loss Expert Panel Specifications to the ACMG/AMP Variant Interpr
Strength: **Supporting**
Desc: Absent/Rare in population databases (absent or ≤0.00007 (0.007%) for autosomal recessive, ≤0.00002 (0.002%) for autosomal dominant). * Background: Rarity or absence in the general population is not robust evidence for pathogenicity, particularly for autosomal recessive disorders. However, the ACMG/AMP Guidelines were devised in such a way that absence or rarity were considered moderate evidence towards pathogenicity, and the framework requires multiple pieces of evidence to classify a variant as likely pathogenic or pathogenic.

## GN006 v2.0.0 [Approved For Release] genes=PAH
VCEP: ClinGen Phenylketonuria Expert Panel Specifications to the ACMG/AMP Variant Inte
Strength: **Supporting**
Desc: * Threshold \<0.0002 (0.02%) The 0.0002 cutoff is based on disease frequency of 1:12,000 and the most common PAH pathogenic variant, R408W, the ExAC frequency is 0.0006594 (ExAC MAF: 0.001109 74/66718 European Non-Finnish) and gnomAD overall: 0.0009056 (gnomAD MAF: 0.001728 219/126,700 European Non-Finnish).

## GN007 v3.1.0 [Released] genes=CDH1
VCEP: ClinGen CDH1 Expert Panel Specifications to the ACMG/AMP Variant Interpretation 
Strength: **Supporting**
Desc: ≤ One out of 100,000 alleles in gnomAD cohort; if present in ≥2 individuals within a subpopulation, must be present in ≤ One out of 50,000 alleles.

## GN008 v3.1.0 [Released] genes=RUNX1
VCEP: ClinGen Myeloid Malignancy Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: **PM2\_Supporting:** Minor allele frequency ≤ 0.00005 with at least 2000 alleles tested around and 20x coverage at the position. Caveat: \*We recommend evaluating PM2\_supporting using the GrpMax FAF when it is available in gnomAD v4.1.0. If a GrpMax FAF value is not available, we recommend requiring that all subpopulations meet the PM2\_supporting threshold.

## GN009 v2.4.0 [Released] genes=TP53
VCEP: ClinGen TP53 Expert Panel Specifications to the ACMG/AMP Variant Interpretation 
Strength: **Supporting**
Desc: This rule should be applied at supporting level. Variant should have an allele frequency of less than 0.00003 (0.003%) in gnomAD or another large sequenced population. If multiple alleles are present within any genetic ancestry group, allele frequency in that group must be \<0.00004 (0.004%). Genetic ancestry groups influenced by founder effects (such as Ashkenazi Jewish, Finnish, Amish, Middle Eastern, and “Remaining”) should be ignored. If the variant being assessed does not meet any population rule codes (PM2, BA1, BS1) **AND** has a total allele frequency >0.00003 with no single genetic ancestry group having multiple alleles with a frequency >0.00004, curators should recalculate the total allele frequency based on the number of alleles with variant allele fraction (VAF) >0.35 to assess whether PM2 may be met after excluding the low VAF alleles which are likely to represent clonal hematopoiesis of indeterminant potential (CHIP) contamination in the database. This can be done by visualizing the “allele balance” for heterozygotes under the genotype quality metrics for a given variant. By hovering over the histogram bars, the number of variant carriers for each bar between 0.35 and 0.65 can be totaled and this can be used to revise the allele count to determine the allele frequency that can be used to assess if PM2\_Supporting can be met. In general, the most recent version of gnomAD should be used when available; however, other population databases or earlier versions of gnomAD may be utilized if they are able to provide information the curator deems necessary for optimal variant classification (e.g, they would provide superior information for a particular variant type; have a larger sample size; or better representation for certain subpopulations, etc.)

## GN010 v2.0.0 [Released] genes=GAA
VCEP: ClinGen Lysosomal Storage Disorders Variant Curation Expert Panel Specifications
Strength: **Moderate**
Desc: Low frequency in population databases. * Minor allele frequency <0.1% (0.001) in all continental populations with >2000 alleles in gnomAD.

## GN011 v2.1.0 [Released] genes=ITGA2B,ITGB3
VCEP: ClinGen Platelet Disorders Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Moderate**
Desc: Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium

## GN011 v2.1.0 [Released] genes=ITGA2B,ITGB3
VCEP: ClinGen Platelet Disorders Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: Prevalence <1/10,000 (<0.0001) alleles in gnomAD.

## GN013 v1.2.0 [Pilot Rules In Prep] genes=LDLR
VCEP: ClinGen Familial Hypercholesterolemia Expert Panel Specifications to the ACMG/AM
Strength: **Moderate**
Desc: Variant has a PopMax MAF ≤0.0002 (0.02%) in gnomAD. Consider exceptions for known founder variants.

## GN014 v1.0.0 [Released] genes=SLC19A3
VCEP: ClinGen Mitochondrial Disease Nuclear and Mitochondrial Expert Panel Specificati
Strength: **Moderate**
Desc: <0.00005 (<0.0050%)

## GN014 v1.0.0 [Released] genes=PDHA1
VCEP: ClinGen Mitochondrial Disease Nuclear and Mitochondrial Expert Panel Specificati
Strength: **Moderate**
Desc: 0.0000092 (<0.00092%)

## GN014 v1.0.0 [Released] genes=POLG
VCEP: ClinGen Mitochondrial Disease Nuclear and Mitochondrial Expert Panel Specificati
Strength: **Moderate**
Desc: <0.0005 (<0.05% )

## GN014 v1.0.0 [Released] genes=ETHE1
VCEP: ClinGen Mitochondrial Disease Nuclear and Mitochondrial Expert Panel Specificati
Strength: **Moderate**
Desc: <0.00002 (<0.0020%)

## GN015 v1.0.0 [Released] genes=
VCEP: ClinGen Mitochondrial Disease Nuclear and Mitochondrial Expert Panel Specificati
Strength: **Supporting**
Desc: Frequency <0.00002 (0.002%, 1/50,000) from controls

## GN016 v2.0.0 [Released] genes=CDKL5
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN016 v2.0.0 [Released] genes=FOXG1
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN016 v2.0.0 [Released] genes=MECP2
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN016 v2.0.0 [Released] genes=SLC9A6
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN016 v2.0.0 [Released] genes=TCF4
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN016 v2.0.0 [Released] genes=UBE3A
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN017 v3.1.0 [Released] genes=HNF1A
VCEP: ClinGen Monogenic Diabetes Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: gnomAD Grpmax FAF ≤ 1:333,000 (≤ 0.000003 or 0.0003%)

## GN018 v1.1.0 [Released] genes=AKT3,MTOR,PIK3CA,PIK3R2
VCEP: ClinGen Brain Malformations Expert Panel Specifications to the ACMG/AMP Variant 
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample ( ≥1).

## GN019 v2.1.0 [Released] genes=MYOC
VCEP: ClinGen Glaucoma Expert Panel Specifications to the ACMG/AMP Variant Interpretat
Strength: **Supporting**
Desc: Allele frequency ≤ 0.0001 in population databases.

## GN020 v1.5.0 [Released] genes=ATM
VCEP: ClinGen Hereditary Breast, Ovarian and Pancreatic Cancer Expert Panel Specificat
Strength: **Supporting**
Desc: Frequency **≤.001%** in gnomAD v4 dataset If n=1 in a single sub population, that is sufficiently rare and PM2\_supporting would apply.

## GN021 v2.1.0 [Released] genes=ACADVL
VCEP: ClinGen ACADVL Expert Panel Specifications to the ACMG/AMP Variant Interpretatio
Strength: **Supporting**
Desc: Variants with a highest population minor allele frequency (MAF) \<0.001 (0.1%) in any continental population with >2000 alleles in gnomAD will meet PM2\_supporting. * Calculated using the Prevalence of 1:100,000, Allelic Contribution of 0.2, Genetic Contribution of 1, and Penetrance of 0.75 to allow for mild VLCADD that may develop in adulthood. This was multiplied by 1.5 to account for mildly pathogenic variants being present in carriers within the population databases. * It is acceptable for an ACADVL variant to be present in controls because VLCAD deficiency is a recessive condition. It is also possible for homozygous ACADVL variants to be present in population databases due to later onset of the condition. If homozygous variants are present, the number should be noted and discussed with an expert.

## GN022 v1.0.0 [Pilot Rules Submitted] genes=FBN1
VCEP: ClinGen FBN1 Expert Panel Specifications to the ACMG/AMP Variant Interpretation 
Strength: **Supporting**
Desc: * Threshold: <5.0E-6 (<0.0005%). * Use the highest ethnic population allele frequency. * Caveat: PVS1 + PM2_Supportive may reach Likely Pathogenic. * Caveat: Do not use Finnish, Ashkenazi Jewish, or “Other” populations in gnomAD. * Minimum amount of studied alleles should be 2000.

## GN023 v1.0.0 [Released] genes=MYO15A,OTOF
VCEP: ClinGen Hearing Loss Expert Panel Specifications to the ACMG/AMP Variant Interpr
Strength: **Supporting**
Desc: Absent/Rare in population databases (absent or ≤0.00007 (0.007%) for autosomal recessive).

## GN024 v1.4.0 [Released] genes=DICER1
VCEP: ClinGen DICER1 and miRNA-Processing Gene Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Allele frequency \<0.000005 across gnomAD with no more than one allele in any subpopulation and at least 20x coverage.

## GN025 v2.0.0 [Released] genes=GATM
VCEP: ClinGen Cerebral Creatine Deficiency Syndromes Expert Panel Specifications to th
Strength: **Supporting**
Desc: **Allele frequency \<0.000055 (\<0.0055%) in all populations in gnomAD.** CCDS VCEP notes: It is acceptable for a GATM variant to be present in controls, if heterozygous, because AGAT-D is a recessive disorder. Homozygotes should not be seen in a population database, such as gnomAD, because the penetrance of this condition in individuals with biallelic pathogenic variants is expected to be 100%. GATM specifications: * All subpopulations in gnomAD must have a maximum allele frequency less than 0.000055 (based on the prevalence of the most common suspected pathogenic variants, c.484+1G>T and p.Arg169Ter) (see Appendix 4). Use the current version recommended by SVI; version number will be stated in classification summary. * Note – PM2 will NOT be used at moderate strength; PM2 will only be applied as a Supporting criterion. * If homozygotes are observed, the variant will meet BS2 (assuming 100% penetrance for an individual with 2 pathogenic variants in trans).

## GN026 v2.0.0 [Released] genes=GAMT
VCEP: ClinGen Cerebral Creatine Deficiency Syndromes Expert Panel Specifications to th
Strength: **Supporting**
Desc: **Allele frequency \<0.0004 (\<0.04%) in all populations in gnomAD.** It is acceptable for a GAMT variant to be present in controls, if heterozygous, because GAMT-D is a recessive disorder. Homozygotes should not be seen in a population database, such as gnomAD, because the penetrance of this condition in individuals with biallelic pathogenic variants is expected to be 100% and the condition is expected to present with severe symptoms early in life. GAMT specifications: * All subpopulations in gnomAD v4.0 must have a maximum allele frequency less than 0.0004 (the highest population minor allele frequency of the most common pathogenic GAMT variant, c.327G>A, in gnomAD). Any variant with a frequency below this cutoff will meet PM2\_Supporting. * If homozygotes are observed, or variant is confirmed in trans with a known pathogenic variant, the variant will meet BS2 (assuming 100% penetrance for an individual with 2 pathogenic variants in trans).

## GN027 v2.0.0 [Released] genes=SLC6A8
VCEP: ClinGen Cerebral Creatine Deficiency Syndromes Expert Panel Specifications to th
Strength: **Supporting**
Desc: Applicable when Grpmax Filtering Allele Frequency is ≤0.00002 (0.002%) AND 0 homo- or hemizygotes are present in the most current version of gnomAD available at the time of curation

## GN032 v6.0.0 [Released] genes=TCF4
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN033 v6.0.0 [Released] genes=SLC9A6
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN034 v6.0.0 [Released] genes=CDKL5
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN035 v6.0.0 [Released] genes=FOXG1
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN036 v6.0.0 [Released] genes=MECP2
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN037 v7.0.0 [Released] genes=UBE3A
VCEP: ClinGen Rett and Angelman-like Disorders Expert Panel Specifications to the ACMG
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Use if absent, zero observations in control databases.

## GN038 v2.3.0 [Released] genes=SHOC2
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN039 v2.3.0 [Pilot Rules In Prep] genes=NRAS
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN040 v2.3.0 [Released] genes=RAF1
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN041 v2.3.0 [Pilot Rules In Prep] genes=SOS1
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN042 v2.3.0 [Pilot Rules In Prep] genes=SOS2
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN043 v2.3.0 [Released] genes=PTPN11
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN044 v2.3.0 [Released] genes=KRAS
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN045 v2.3.0 [Released] genes=MAP2K1
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN046 v2.3.0 [Released] genes=HRAS
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN047 v2.3.0 [Released] genes=RIT1
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN048 v2.3.0 [Released] genes=MAP2K2
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN049 v2.3.0 [Released] genes=BRAF
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN067 v2.0.0 [Released] genes=SCN1A
VCEP: ClinGen Epilepsy Sodium Channel Expert Panel Specifications to the ACMG/AMP Vari
Strength: **Supporting**
Desc: One or fewer alleles, if a minimum of 10,000 alleles assessed in population databases, such as the Genome Aggregation Database (gnomAD). Caveat: Population data for indels may be poorly called by next generation sequencing.

## GN068 v2.0.0 [Released] genes=SCN2A
VCEP: ClinGen Epilepsy Sodium Channel Expert Panel Specifications to the ACMG/AMP Vari
Strength: **Supporting**
Desc: One or fewer alleles, if a minimum of 10,000 alleles assessed in population databases, such as the Genome Aggregation Database (gnomAD). Caveat: Population data for indels may be poorly called by next generation sequencing.

## GN069 v2.1.0 [Released] genes=SCN3A
VCEP: ClinGen Epilepsy Sodium Channel Expert Panel Specifications to the ACMG/AMP Vari
Strength: **Supporting**
Desc: One or fewer alleles, if a minimum of 10,000 alleles assessed in population databases, such as the Genome Aggregation Database (gnomAD). Caveat: Population data for indels may be poorly called by next generation sequencing.

## GN070 v2.0.0 [Released] genes=SCN8A
VCEP: ClinGen Epilepsy Sodium Channel Expert Panel Specifications to the ACMG/AMP Vari
Strength: **Supporting**
Desc: One or fewer alleles, if a minimum of 10,000 alleles assessed in population databases, such as the Genome Aggregation Database (gnomAD). Caveat: Population data for indels may be poorly called by next generation sequencing.

## GN071 v2.0.0 [Released] genes=F8
VCEP: ClinGen Coagulation Factor Deficiency Expert Panel Specifications to the ACMG/AM
Strength: **Supporting**
Desc: Variant must be absent in males in population databases, such as gnomAD.

## GN076 v2.0.0 [Released] genes=SCN1B
VCEP: ClinGen Epilepsy Sodium Channel Expert Panel Specifications to the ACMG/AMP Vari
Strength: **Supporting**
Desc: One or fewer alleles, if a minimum of 10,000 alleles assessed in population databases, such as the Genome Aggregation Database (gnomAD). Caveat: Population data for indels may be poorly called by next generation sequencing.

## GN077 v1.2.0 [Released] genes=PALB2
VCEP: ClinGen Hereditary Breast, Ovarian and Pancreatic Cancer Expert Panel Specificat
Strength: **Supporting**
Desc: Frequency ≤ 1/300,000 (**0.000333%**) in gnomAD v4 dataset

## GN078 v1.1.0 [Released] genes=VHL
VCEP: ClinGen VHL Expert Panel Specifications to the ACMG/AMP Variant Interpretation G
Strength: **Supporting**
Desc: PM2\_Supporting can be applied for variants either absent from gnomAD or with \<= 0.00000156 (0.000156%) GroupMax Filtering Allele Frequency in gnomAD (based on gnomAD v4 release). If no GroupMax Filtering Allele Frequency is calculated (ex. due to a single variant present), PM2\_Supporting may also be applied.

## GN079 v1.1.0 [Released] genes=GP1BA
VCEP: ClinGen Platelet Disorders Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: gnomAD MAF of less than or equal to 0.0001114.

## GN080 v2.1.0 [Released] genes=F9
VCEP: ClinGen Coagulation Factor Deficiency Expert Panel Specifications to the ACMG/AM
Strength: **Supporting**
Desc: Variant must be absent in males in population databases, such as gnomAD.

## GN081 v1.0.0 [Released] genes=VWF
VCEP: ClinGen von Willebrand Disease  Expert Panel Specifications to the ACMG/AMP Vari
Strength: **Supporting**
Desc: Use code for variants with a popmax MAF of \<0.0001 in gnomAD.

## GN082 v1.1.0 [Released] genes=GP1BB
VCEP: ClinGen Platelet Disorders Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: gnomAD MAF of less than or equal to 0.00006517.

## GN083 v1.1.0 [Released] genes=GP9
VCEP: ClinGen Platelet Disorders Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: gnomAD MAF of less than or equal to 0.0000329.

## GN084 v1.1.0 [Pilot Rules In Prep] genes=SERPINC1
VCEP: ClinGen Thrombosis Expert Panel Specifications to the ACMG/AMP Variant Interpret
Strength: **Supporting**
Desc: Use code for variants with a popmax MAF of \<0.00002 in gnomAD.

## GN085 v4.0.0 [Released] genes=HNF4A
VCEP: ClinGen Monogenic Diabetes Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: gnomAD Grpmax FAF ≤ 1:333,000 (≤ 0.000003 or 0.0003%)

## GN086 v3.1.0 [Released] genes=GCK
VCEP: ClinGen Monogenic Diabetes Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: gnomAD Grpmax FAF ≤ 1:333,000 (≤ 0.000003 or 0.0003%)

## GN087 v1.4.0 [Released] genes=MRAS
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN088 v1.3.0 [Released] genes=RMRP
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: Absent in population databases (or at extremely low frequency if recessive). * Downgraded to PM2\_Supporting. * gnomAD popmax filtering allele frequency \<0.0000447 The applicability of PM2 to suspected founder variants with allele frequencies exceeding the PM2 threshold will be evaluated on a case-by-case basis by the VCEP.

## GN089 v2.1.0 [Released] genes=APC
VCEP: ClinGen InSiGHT Hereditary Colorectal Cancer/Polyposis Expert Panel Specificatio
Strength: **Supporting**
Desc: Rare in controls, defined by an allele frequency ≤ 0.0003% (0.000003) if the allele count is > 1 OR by an allele frequency \< 0.001% (0.00001) if the allele count is ≤ 1.

## GN090 v1.0.0 [Released] genes=VWF
VCEP: ClinGen von Willebrand Disease  Expert Panel Specifications to the ACMG/AMP Vari
Strength: **Supporting**
Desc: Use code for variants with a popmax MAF of \<0.005 in gnomAD.

## GN091 v1.2.0 [Released] genes=IDUA
VCEP: ClinGen Lysosomal Diseases Expert Panel Specifications to the ACMG/AMP Variant I
Strength: **Supporting**
Desc: This criterion will be applied at the supporting level based on [guidance](https://www.clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf) from the ClinGen Sequence Variant Interpretation Working Group. Minor allele frequency \<0.025% (0.00025) in any continental population with >2000 alleles in the most recent version of gnomAD (version # will be stated in the written summary). Variants may be observed in the homozygous state because MPS1 can present in adulthood, and some variants may be hypomorphic. However, the presence and number of homozygotes should be noted.

## GN092 v1.2.0 [Released] genes=BRCA1
VCEP: ClinGen ENIGMA BRCA1 and BRCA2 Expert Panel Specifications to the ACMG/AMP Varia
Strength: **Supporting**
Desc: Absent from controls in an outbred population, from gnomAD v2.1 (non-cancer, exome only subset) and gnomAD v3.1 (non-cancer). Region around the variant must have an average read depth ≥25. See Appendix G for details.

## GN094 v1.3.0 [Released] genes=LZTR1
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD). For variants in _LZTR1_, PM2\_P ≤0.0025% may be applied to support AR disease.

## GN095 v1.0.0 [Released] genes=MYBPC3
VCEP: ClinGen Cardiomyopathy Expert Panel Specifications to the ACMG/AMP Variant Inter
Strength: **Supporting**
Desc: The values used to calculate the PM2 thresholds were derived from studies in Northern European populations that have been relatively well-characterized with regards to disease prevalence and variant spectrum. These thresholds can be applied to any population where disease prevalence is considered comparable (1/500 or lower), where the most frequent pathogenic variant accounts for no more than 2% of cases (e.g., has an allele frequency of ≤0.02 in cases based on the upper bound of 95% CI), and where the penetrance of a pathogenic variant is expected to be at least 50% (Kelly _et al._ 2018[<sup>15</sup>](#pmid_29300372)). A threshold of **≤0.00004** in the subpopulation with the highest frequency when using the upper bound of the 95% CI activates this rule. 1. Alternatively, this is equivalent to the variant NOT being observed more than once (≤1 allele) in gnomAD v.2.1.1 in one of the non-founder populations (e.g., absence required from the Other and Ashkenazi Jewish subpopulations). 2. Applying a threshold of ≤0.00004 (upper bound of 95% CI of the allele frequency in gnomAD) is equivalent to the variant being seen in a single subpopulation and that subpopulation meets any of the following: * **Allele Count (AC) in Allele Number (AN)** * ≤1 in ≥120,000 * ≤2 in ≥160,000 * ≤3 in ≥195,000 * ≤4 in ≥230,000 gnomAD is the preferred database for this calculation, but currently only displays the filtering allele frequency (FAF), which is equivalent to a lower bound estimate of the 95% CI, when the upper bound is what is needed. * Confidence interval tools, such as [Confit-de-MAF](https://www.genecalculators.net/confit-de-maf.html), can be used to determine the upper bound of the 95% CI of the observed allele frequency. Due to current technical limitations of next generation sequencing technologies, minor allele frequencies for complex variants (e.g., large indels) may not be accurately represented in population databases. Caution should be used when a variant is only identified, or over-represented, in one of the smaller gnomAD populations, as the gnomAD allele frequencies may not accurately represent the true population frequency. Population databases may contain affected or pre-symptomatic individuals for diseases with reduced penetrance/variable onset.

## GN097 v1.2.0 [Released] genes=BRCA2
VCEP: ClinGen ENIGMA BRCA1 and BRCA2 Expert Panel Specifications to the ACMG/AMP Varia
Strength: **Supporting**
Desc: Absent from controls in an outbred population, from gnomAD v2.1 (non-cancer, exome only subset) and gnomAD v3.1 (non-cancer). Region around the variant must have an average read depth ≥25. See Appendix G for details.

## GN098 v1.0.0 [Released] genes=TNNI3
VCEP: ClinGen Cardiomyopathy Expert Panel Specifications to the ACMG/AMP Variant Inter
Strength: **Supporting**
Desc: The values used to calculate the PM2 thresholds were derived from studies in Northern European populations that have been relatively well-characterized with regards to disease prevalence and variant spectrum. These thresholds can be applied to any population where disease prevalence is considered comparable (1/500 or lower), where the most frequent pathogenic variant accounts for no more than 2% of cases (e.g., has an allele frequency of ≤0.02 in cases based on the upper bound of 95% CI), and where the penetrance of a pathogenic variant is expected to be at least 50% (Kelly _et al._ 2018[<sup>11</sup>](#pmid_29300372)). A threshold of **≤0.00004** in the subpopulation with the highest frequency when using the upper bound of the 95% CI activates this rule. 1. Alternatively, this is equivalent to the variant NOT being observed more than once (≤1 allele) in gnomAD v.2.1.1 in one of the non-founder populations (e.g., absence required from the Other and Ashkenazi Jewish subpopulations). 2. Applying a threshold of ≤0.00004 (upper bound of 95% CI of the allele frequency in gnomAD) is equivalent to the variant being seen in a single subpopulation and that subpopulation meets any of the following: * **Allele Count (AC) in Allele Number (AN)** * ≤1 in ≥120,000 * ≤2 in ≥160,000 * ≤3 in ≥195,000 * ≤4 in ≥230,000 gnomAD is the preferred database for this calculation, but currently only displays the filtering allele frequency (FAF), which is equivalent to a lower bound estimate of the 95% CI, when the upper bound is what is needed. * Confidence interval tools, such as [Confit-de-MAF](https://www.genecalculators.net/confit-de-maf.html), can be used to determine the upper bound of the 95% CI of the observed allele frequency. Due to current technical limitations of next generation sequencing technologies, minor allele frequencies for complex variants (e.g., large indels) may not be accurately represented in population databases. Caution should be used when a variant is only identified, or over-represented, in one of the smaller gnomAD populations, as the gnomAD allele frequencies may not accurately represent the true population frequency. Population databases may contain affected or pre-symptomatic individuals for diseases with reduced penetrance/variable onset.

## GN099 v1.0.0 [Released] genes=TNNT2
VCEP: ClinGen Cardiomyopathy Expert Panel Specifications to the ACMG/AMP Variant Inter
Strength: **Supporting**
Desc: The values used to calculate the PM2 thresholds were derived from studies in Northern European populations that have been relatively well-characterized with regards to disease prevalence and variant spectrum. These thresholds can be applied to any population where disease prevalence is considered comparable (1/500 or lower), where the most frequent pathogenic variant accounts for no more than 2% of cases (e.g., has an allele frequency of ≤0.02 in cases based on the upper bound of 95% CI), and where the penetrance of a pathogenic variant is expected to be at least 50% (Kelly _et al._ 2018[<sup>11</sup>](#pmid_29300372)). A threshold of **≤0.00004** in the subpopulation with the highest frequency when using the upper bound of the 95% CI activates this rule. 1. Alternatively, this is equivalent to the variant NOT being observed more than once (≤1 allele) in gnomAD v.2.1.1 in one of the non-founder populations (e.g., absence required from the Other and Ashkenazi Jewish subpopulations). 2. Applying a threshold of ≤0.00004 (upper bound of 95% CI of the allele frequency in gnomAD) is equivalent to the variant being seen in a single subpopulation and that subpopulation meets any of the following: * **Allele Count (AC) in Allele Number (AN)** * ≤1 in ≥120,000 * ≤2 in ≥160,000 * ≤3 in ≥195,000 * ≤4 in ≥230,000 gnomAD is the preferred database for this calculation, but currently only displays the filtering allele frequency (FAF), which is equivalent to a lower bound estimate of the 95% CI, when the upper bound is what is needed. * Confidence interval tools, such as [Confit-de-MAF](https://www.genecalculators.net/confit-de-maf.html), can be used to determine the upper bound of the 95% CI of the observed allele frequency. Due to current technical limitations of next generation sequencing technologies, minor allele frequencies for complex variants (e.g., large indels) may not be accurately represented in population databases. Caution should be used when a variant is only identified, or over-represented, in one of the smaller gnomAD populations, as the gnomAD allele frequencies may not accurately represent the true population frequency. Population databases may contain affected or pre-symptomatic individuals for diseases with reduced penetrance/variable onset.

## GN100 v1.0.0 [Released] genes=TPM1
VCEP: ClinGen Cardiomyopathy Expert Panel Specifications to the ACMG/AMP Variant Inter
Strength: **Supporting**
Desc: The values used to calculate the PM2 thresholds were derived from studies in Northern European populations that have been relatively well-characterized with regards to disease prevalence and variant spectrum. These thresholds can be applied to any population where disease prevalence is considered comparable (1/500 or lower), where the most frequent pathogenic variant accounts for no more than 2% of cases (e.g., has an allele frequency of ≤0.02 in cases based on the upper bound of 95% CI), and where the penetrance of a pathogenic variant is expected to be at least 50% (Kelly _et al._ 2018[<sup>10</sup>](#pmid_29300372)). A threshold of **≤0.00004** in the subpopulation with the highest frequency when using the upper bound of the 95% CI activates this rule. 1. Alternatively, this is equivalent to the variant NOT being observed more than once (≤1 allele) in gnomAD v.2.1.1 in one of the non-founder populations (e.g., absence required from the Other and Ashkenazi Jewish subpopulations). 2. Applying a threshold of ≤0.00004 (upper bound of 95% CI of the allele frequency in gnomAD) is equivalent to the variant being seen in a single subpopulation and that subpopulation meets any of the following: * **Allele Count (AC) in Allele Number (AN)** * ≤1 in ≥120,000 * ≤2 in ≥160,000 * ≤3 in ≥195,000 * ≤4 in ≥230,000 gnomAD is the preferred database for this calculation, but currently only displays the filtering allele frequency (FAF), which is equivalent to a lower bound estimate of the 95% CI, when the upper bound is what is needed. * Confidence interval tools, such as [Confit-de-MAF](https://www.genecalculators.net/confit-de-maf.html), can be used to determine the upper bound of the 95% CI of the observed allele frequency. Due to current technical limitations of next generation sequencing technologies, minor allele frequencies for complex variants (e.g., large indels) may not be accurately represented in population databases. Caution should be used when a variant is only identified, or over-represented, in one of the smaller gnomAD populations, as the gnomAD allele frequencies may not accurately represent the true population frequency. Population databases may contain affected or pre-symptomatic individuals for diseases with reduced penetrance/variable onset.

## GN101 v1.0.0 [Released] genes=ACTC1
VCEP: ClinGen Cardiomyopathy Expert Panel Specifications to the ACMG/AMP Variant Inter
Strength: **Supporting**
Desc: The values used to calculate the PM2 thresholds were derived from studies in Northern European populations that have been relatively well-characterized with regards to disease prevalence and variant spectrum. These thresholds can be applied to any population where disease prevalence is considered comparable (1/500 or lower), where the most frequent pathogenic variant accounts for no more than 2% of cases (e.g., has an allele frequency of ≤0.02 in cases based on the upper bound of 95% CI), and where the penetrance of a pathogenic variant is expected to be at least 50% (Kelly _et al._ 2018[<sup>10</sup>](#pmid_29300372)). A threshold of **≤0.00004** in the subpopulation with the highest frequency when using the upper bound of the 95% CI activates this rule. 1. Alternatively, this is equivalent to the variant NOT being observed more than once (≤1 allele) in gnomAD v.2.1.1 in one of the non-founder populations (e.g., absence required from the Other and Ashkenazi Jewish subpopulations). 2. Applying a threshold of ≤0.00004 (upper bound of 95% CI of the allele frequency in gnomAD) is equivalent to the variant being seen in a single subpopulation and that subpopulation meets any of the following: * **Allele Count (AC) in Allele Number (AN)** * ≤1 in ≥120,000 * ≤2 in ≥160,000 * ≤3 in ≥195,000 * ≤4 in ≥230,000 gnomAD is the preferred database for this calculation, but currently only displays the filtering allele frequency (FAF), which is equivalent to a lower bound estimate of the 95% CI, when the upper bound is what is needed. * Confidence interval tools, such as [Confit-de-MAF](https://www.genecalculators.net/confit-de-maf.html), can be used to determine the upper bound of the 95% CI of the observed allele frequency. Due to current technical limitations of next generation sequencing technologies, minor allele frequencies for complex variants (e.g., large indels) may not be accurately represented in population databases. Caution should be used when a variant is only identified, or over-represented, in one of the smaller gnomAD populations, as the gnomAD allele frequencies may not accurately represent the true population frequency. Population databases may contain affected or pre-symptomatic individuals for diseases with reduced penetrance/variable onset.

## GN102 v1.0.0 [Released] genes=MYL2
VCEP: ClinGen Cardiomyopathy Expert Panel Specifications to the ACMG/AMP Variant Inter
Strength: **Supporting**
Desc: The values used to calculate the PM2 thresholds were derived from studies in Northern European populations that have been relatively well-characterized with regards to disease prevalence and variant spectrum. These thresholds can be applied to any population where disease prevalence is considered comparable (1/500 or lower), where the most frequent pathogenic variant accounts for no more than 2% of cases (e.g., has an allele frequency of ≤0.02 in cases based on the upper bound of 95% CI), and where the penetrance of a pathogenic variant is expected to be at least 50% (Kelly _et al._ 2018[<sup>10</sup>](#pmid_29300372)). A threshold of **≤0.00004** in the subpopulation with the highest frequency when using the upper bound of the 95% CI activates this rule. 1. Alternatively, this is equivalent to the variant NOT being observed more than once (≤1 allele) in gnomAD v.2.1.1 in one of the non-founder populations (e.g., absence required from the Other and Ashkenazi Jewish subpopulations). 2. Applying a threshold of ≤0.00004 (upper bound of 95% CI of the allele frequency in gnomAD) is equivalent to the variant being seen in a single subpopulation and that subpopulation meets any of the following: * **Allele Count (AC) in Allele Number (AN)** * ≤1 in ≥120,000 * ≤2 in ≥160,000 * ≤3 in ≥195,000 * ≤4 in ≥230,000 gnomAD is the preferred database for this calculation, but currently only displays the filtering allele frequency (FAF), which is equivalent to a lower bound estimate of the 95% CI, when the upper bound is what is needed. * Confidence interval tools, such as [Confit-de-MAF](https://www.genecalculators.net/confit-de-maf.html), can be used to determine the upper bound of the 95% CI of the observed allele frequency. Due to current technical limitations of next generation sequencing technologies, minor allele frequencies for complex variants (e.g., large indels) may not be accurately represented in population databases. Caution should be used when a variant is only identified, or over-represented, in one of the smaller gnomAD populations, as the gnomAD allele frequencies may not accurately represent the true population frequency. Population databases may contain affected or pre-symptomatic individuals for diseases with reduced penetrance/variable onset.

## GN103 v1.0.0 [Released] genes=MYL3
VCEP: ClinGen Cardiomyopathy Expert Panel Specifications to the ACMG/AMP Variant Inter
Strength: **Supporting**
Desc: The values used to calculate the PM2 thresholds were derived from studies in Northern European populations that have been relatively well-characterized with regards to disease prevalence and variant spectrum. These thresholds can be applied to any population where disease prevalence is considered comparable (1/500 or lower), where the most frequent pathogenic variant accounts for no more than 2% of cases (e.g., has an allele frequency of ≤0.02 in cases based on the upper bound of 95% CI), and where the penetrance of a pathogenic variant is expected to be at least 50% (Kelly _et al._ 2018[<sup>10</sup>](#pmid_29300372)). A threshold of **≤0.00004** in the subpopulation with the highest frequency when using the upper bound of the 95% CI activates this rule. 1. Alternatively, this is equivalent to the variant NOT being observed more than once (≤1 allele) in gnomAD v.2.1.1 in one of the non-founder populations (e.g., absence required from the Other and Ashkenazi Jewish subpopulations). 2. Applying a threshold of ≤0.00004 (upper bound of 95% CI of the allele frequency in gnomAD) is equivalent to the variant being seen in a single subpopulation and that subpopulation meets any of the following: * **Allele Count (AC) in Allele Number (AN)** * ≤1 in ≥120,000 * ≤2 in ≥160,000 * ≤3 in ≥195,000 * ≤4 in ≥230,000 gnomAD is the preferred database for this calculation, but currently only displays the filtering allele frequency (FAF), which is equivalent to a lower bound estimate of the 95% CI, when the upper bound is what is needed. * Confidence interval tools, such as [Confit-de-MAF](https://www.genecalculators.net/confit-de-maf.html), can be used to determine the upper bound of the 95% CI of the observed allele frequency. Due to current technical limitations of next generation sequencing technologies, minor allele frequencies for complex variants (e.g., large indels) may not be accurately represented in population databases. Caution should be used when a variant is only identified, or over-represented, in one of the smaller gnomAD populations, as the gnomAD allele frequencies may not accurately represent the true population frequency. Population databases may contain affected or pre-symptomatic individuals for diseases with reduced penetrance/variable onset.

## GN104 v1.0.0 [Released] genes=CYP1B1
VCEP: ClinGen Glaucoma Expert Panel Specifications to the ACMG/AMP Variant Interpretat
Strength: **Supporting**
Desc: Allele frequency ≤ 0.0005 in population databases.

## GN105 v1.0.0 [Released] genes=ABCD1
VCEP: ClinGen Peroxisomal Disorders Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: * PM2\_Supporting can be applied if the variant is absent in hemizygotes AND has a maximum allele frequency of \<0.00017% (0.0000017) in heterozygotes in the latest version of gnomAD. Use the highest population MAF from a non-bottleneck population from the latest version of gnomAD.

## GN106 v1.0.0 [Released] genes=RPGR
VCEP: ClinGen X-linked Inherited Retinal Disease Expert Panel Specifications to the AC
Strength: **Supporting**
Desc: Allele frequency in males ≤ 0.00005 (≤5x10<sup>-5</sup>) in population databases. * Highest allele frequency in a subpopulation should be used to assess this.

## GN112 v1.0.0 [Released] genes=KCNQ1
VCEP: ClinGen Potassium Channel Arrhythmia Expert Panel Specifications to the ACMG/AMP
Strength: **Supporting**
Desc: Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes or Exome Aggregation Consortium. * Supporting level only * Maximum allele frequency in gnomAD (in one of the 5 continental populations; African/African-American, East Asian, European non-Finnish, Latino/Admixed-American, or South Asian) \<0.00001 (0.001%)

## GN113 v2.3.0 [Released] genes=FOXN1
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: gnomAD Grpmax filtering allele frequency ≤0.00002412

## GN114 v2.1.0 [Pilot Rules Submitted] genes=ADA
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: gnomAD popmax filtering allele frequency \<0.0001742 * An additional requirement is that **no homozygotes** have been observed in gnomAD.

## GN115 v2.0.0 [Pilot Rules Submitted] genes=MLH1
VCEP: ClinGen InSiGHT Hereditary Colorectal Cancer/Polyposis Expert Panel Specificatio
Strength: **Supporting**
Desc: Absent/extremely rare allele frequency \<0.00002 (\<1 in 50,000 alleles ) in gnomAD v4 dataset

## GN116 v2.2.0 [Released] genes=DCLRE1C
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: gnomAD popmax filtering allele frequency \<0.00003266 * An additional requirement is that **no homozygotes** have been observed in gnomAD.

## GN119 v2.2.0 [Released] genes=IL7R
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: gnomAD popmax filtering allele frequency \<0.00004129. * An additional requirement is that **no homozygotes** have been observed in gnomAD.

## GN120 v1.0.0 [Released] genes=RPE65
VCEP: ClinGen Leber Congenital Amaurosis/early onset Retinal Dystrophy Expert Panel Sp
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Used if the gnomAD PopMax Filtering Allele Frequency (FAF) is ≤ 2.0 x 10<sup>-4</sup>.

## GN121 v2.3.0 [Released] genes=JAK3
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: gnomAD popmax filtering allele frequency \<0.000115 * An additional requirement is that **no homozygotes** have been observed in gnomAD.

## GN122 v1.0.0 [Pilot Rules In Prep] genes=CTLA4
VCEP: ClinGen Antibody Deficiencies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: * Supporting level only * Met for total allele frequency lower than 1.43 x 10<sup>-7</sup> (0.000000143) across all populations in gnomAD v4.1.0. * Threshold is based on the experts’ estimate of CTLA-4 insufficiency prevalence of 1/200,000 – 1/1,000,000 people and 45-70% penetrance. The lower end of the prevalence estimate (1 in 1,000,000) and the higher end of the penetrance estimate (70%) were used for this calculation. Allelic heterogeneity of 1 and genetic heterogeneity of 1 were also assumed for the calculation.

## GN123 v2.2.0 [Released] genes=RAG1
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: gnomAD popmax filtering allele frequency \<0.000102 * An additional requirement is that **no homozygotes** have been observed in gnomAD.

## GN124 v2.2.0 [Released] genes=RAG2
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: gnomAD popmax filtering allele frequency \<0.0000588 * An additional requirement is that **no homozygotes** have been observed in gnomAD.

## GN125 v2.0.0 [Released] genes=BMPR2
VCEP: ClinGen Pulmonary Hypertension Expert Panel Specifications to the ACMG/AMP Varia
Strength: **Supporting**
Desc: Present at \<0.01% among gnomAD controls, using the subpopulation with the highest frequency and at least 1,000 allele counts. Caveat: Population data for indels may be poorly called by next generation sequencing.

## GN126 v1.0.0 [Released] genes=RS1
VCEP: ClinGen X-linked Inherited Retinal Disease Expert Panel Specifications to the AC
Strength: **Supporting**
Desc: At low frequency in males in population databases. Use \<2.0x10<sup>-6 </sup> for cut off. This is defined relative to the BA1 cutoff.

## GN127 v1.3.0 [Pilot Rules In Prep] genes=RRAS2
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN128 v1.3.0 [Released] genes=PPP1CB
VCEP: ClinGen RASopathy Expert Panel Specifications to the ACMG/AMP Variant Interpreta
Strength: **Supporting**
Desc: The variant must be absent from controls (gnomAD).

## GN129 v2.2.0 [Released] genes=IL2RG
VCEP: ClinGen Severe Combined Immunodeficiency Disease  Expert Panel Specifications to
Strength: **Supporting**
Desc: Strength modification based on an abnormal result in at least one approved _in vitro_ assay.

## GN135 v1.1.0 [Released] genes=ACVRL1
VCEP: ClinGen Hereditary Hemorrhagic Telangiectasia Expert Panel Specifications to the
Strength: **Supporting**
Desc: \<6 total alleles in gnomAD or \<0.00004 (0.004%) in gnomAD subpopulations.

## GN136 v1.1.0 [Released] genes=ENG
VCEP: ClinGen Hereditary Hemorrhagic Telangiectasia Expert Panel Specifications to the
Strength: **Supporting**
Desc: \<6 total alleles in gnomAD or \<0.00004 (0.004%) in gnomAD subpopulations.

## GN137 v2.0.0 [Pilot Rules Submitted] genes=MSH2
VCEP: ClinGen InSiGHT Hereditary Colorectal Cancer/Polyposis Expert Panel Specificatio
Strength: **Supporting**
Desc: Absent/extremely rare allele frequency \<0.00002 (\<1 in 50,000 alleles ) in gnomAD v4 dataset

## GN138 v2.0.0 [Pilot Rules Submitted] genes=MSH6
VCEP: ClinGen InSiGHT Hereditary Colorectal Cancer/Polyposis Expert Panel Specificatio
Strength: **Supporting**
Desc: Absent/extremely rare allele frequency \<0.00002 (\<1 in 50,000 alleles ) in gnomAD v4 dataset

## GN139 v2.0.0 [Pilot Rules Submitted] genes=PMS2
VCEP: ClinGen InSiGHT Hereditary Colorectal Cancer/Polyposis Expert Panel Specificatio
Strength: **Supporting**
Desc: Absent/extremely rare allele frequency \<0.00002 (\<1 in 50,000 alleles ) in gnomAD v4 dataset

## GN141 v1.0.0 [Released] genes=PIK3CD
VCEP: ClinGen Antibody Deficiencies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: * Downgraded to PM2\_Supporting * Applicable to variants with a total allele frequency \<0.00000132 across all populations in gnomAD v4.1.0. * Maximum credible population allele frequency threshold determined using Whiffin/Ware calculator ([https://www.cardiodb.org/allelefrequencyapp/](https://www.cardiodb.org/allelefrequencyapp/)) and the following estimated parameters (with the prevalence estimated for autosomal dominant PIK3CD-related immune disease): * Prevalence: 1 in 4000 (a conservative estimate for primary immunodeficiency diseases from PMID: 17577648, PMID: 23201919) * Allelic heterogeneity: 1 * Genetic heterogeneity: 1 * Penetrance: 0.95 (a conservative estimate based on multiple reports of incomplete but nearly complete penetrance in PMID: 27555459, PMID: 36749229, PMID: 37390899)

## GN146 v1.0.0 [Released] genes=NEB
VCEP: ClinGen Congenital Myopathies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: PM2\_Supporting may be applied if the minor allele frequency in population databases of at least 2000 alleles is ≤ 0.0000559. 1 allele is allowed.

## GN147 v2.0.0 [Released] genes=ACTA1
VCEP: ClinGen Congenital Myopathies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: PM2\_Supporting may be applied if the minor allele frequency in population databases of at least 2000 alleles is absent (1 allele allowed) for autosomal dominant

## GN148 v1.0.0 [Released] genes=DNM2
VCEP: ClinGen Congenital Myopathies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: PM2\_Supporting may be applied if the minor allele frequency in population databases of at least 2000 alleles is absent (1 allele allowed)

## GN149 v1.0.0 [Released] genes=MTM1
VCEP: ClinGen Congenital Myopathies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: PM2\_Supporting may be applied if the minor allele frequency in population databases of at least 2000 alleles is absent (1 observation allowed in females only)

## GN150 v2.0.0 [Released] genes=RYR1
VCEP: ClinGen Congenital Myopathies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: PM2\_Supporting may be applied if the minor allele frequency in population databases of at least 2000 alleles is absent (1 allele allowed) for autosomal dominant

## GN156 v1.0.0 [Released] genes=OTC
VCEP: ClinGen Urea Cycle Disorders Expert Panel Specifications to the ACMG/AMP Variant
Strength: **Supporting**
Desc: Applicable for variants in OTC with Grpmax Filtering Allele Frequency \<0.000015 (0.0015%) AND ≤1 homo- or hemizygote in the most current version of gnomAD available at the time of curation. Rationale: The most common pathogenic variant in population databases is p.Arg40Cys, which is associated with late onset OTC Deficiency (PMID: 23209112, 7860066, 11260212, others) and present in 17 heterozygotes and 6 hemizygotes in gnomAD(v4.0.0) (Mino Allele frequency=0.001586% in European populations). Other commonly reported pathogenic variants (p.Arg277Trp, p.Arg141Gln, p.Arg141Ter) are rare or absent in population databases, therefore a threshold of 0.0015% is set for PM2\_Supporting.

## GN158 v1.0.0 [Released] genes=GALT
VCEP: ClinGen Galactosemia Expert Panel Specifications to the ACMG/AMP Variant Interpr
Strength: **Supporting**
Desc: * Can only be used at the supporting level * Use the PopMax filtering allele frequency (FAF) from gnomAD. * Use if variant is present at ≤ 0.0005 or 0.05% * One order of magnitude below BS1 * A curated list of _GALT_ variants known to be pathogenic despite not fulfilling the PM2 criterion has been developed (see PM2 exception list)

## GN160 v1.0.0 [Released] genes=PIK3R1
VCEP: ClinGen Antibody Deficiencies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: * Used at PM2\_Supporting strength. * Applicable to variants with a total allele frequency \<0.00000132 across all populations in gnomAD v4.1.0. * Population allele frequency threshold has been determined using the Whiffin/Ware calculator ([https://www.cardiodb.org/allelefrequencyapp/](https://www.cardiodb.org/allelefrequencyapp/)) and the following estimated parameters (with the prevalence estimated for primary immunodeficiency and the inheritance tailored to the autosomal dominant mode of inheritance): * Prevalence: 1 in 4000 * Allelic heterogeneity: 1 * Genetic heterogeneity: 1 * Penetrance: 0.95

## GN164 v1.0.0 [Pilot Rules In Prep] genes=ABCA4
VCEP: ClinGen ABCA4 Expert Panel Specifications to the ACMG/AMP Variant Interpretation
Strength: **Supporting**
Desc: Total MAF \<0.0001 in gnomAD.

## GN167 v1.0.0 [Released] genes=GUCY2D
VCEP: ClinGen Leber Congenital Amaurosis/early onset Retinal Dystrophy Expert Panel Sp
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Used if the gnomAD total allele frequency is ≤ 4.0 x 10<sup>-4</sup>.

## GN169 v1.0.0 [Released] genes=ACTA1
VCEP: ClinGen Congenital Myopathies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: PM2\_Supporting may be applied if the minor allele frequency in population databases of at least 2000 alleles is ≤ 0.000005 for autosomal recessive

## GN170 v1.0.0 [Released] genes=HBB
VCEP: ClinGen Hemoglobinopathy Expert Panel Specifications to the ACMG/AMP Variant Int
Strength: **Supporting**
Desc: Allele frequency \<0.0001 (0.01%) in gnomAD.

## GN173 v1.0.0 [Released] genes=HBA2
VCEP: ClinGen Hemoglobinopathy Expert Panel Specifications to the ACMG/AMP Variant Int
Strength: **Supporting**
Desc: Allele frequency \<0.0001 (0.01%) in gnomAD.

## GN179 v2.0.0 [Released] genes=RYR1
VCEP: ClinGen Congenital Myopathies Expert Panel Specifications to the ACMG/AMP Varian
Strength: **Supporting**
Desc: PM2\_Supporting may be applied if the minor allele frequency in population databases of at least 2000 alleles is ≤ 0.00000697 for autosomal recessive

## GN180 v2.0.0 [Released] genes=DYSF
VCEP: ClinGen Limb Girdle Muscular Dystrophy Expert Panel Specifications to the ACMG/A
Strength: **Supporting**
Desc: Apply if the Grpmax variant allele frequency / upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency in gnomAD is \<0.0001. Do not use data for which the variant does not pass quality control filters. * If only 1 or 2 variant alleles are present in the Grpmax population, use the Grpmax variant allele frequency * If at least 3 variant alleles are present in the Grpmax population, use the upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency Grpmax refers to the gnomAD subpopulation with the highest variant allele frequency. Use large, non-bottlenecked genetic ancestry groups for the Grpmax; avoid using the Amish, Ashkenazi Jewish, European Finnish, and Remaining Individuals groups as well as the genomes-only data for the Middle Eastern group. The upper bound of the 95% CI must be calculated using variant allele numbers and counts from gnomAD. Confidence interval tools, such as Confit-de-MAF ([https://www.genecalculators.net/confit-de-maf.html](https://www.genecalculators.net/confit-de-maf.html)), can be used. Use the gnomAD version with the largest allele number. For larger deletions or duplications that may not be well represented in gnomAD (e.g., single- or multi-exon events), also confirm the variant is not common in gnomAD SVs, gnomAD CNVs, or the Database of Genomic Variants (DGV) ([https://dgv.tcag.ca/dgv/app/home](https://dgv.tcag.ca/dgv/app/home)).

## GN184 v2.0.0 [Released] genes=SGCB
VCEP: ClinGen Limb Girdle Muscular Dystrophy Expert Panel Specifications to the ACMG/A
Strength: **Supporting**
Desc: Apply if the Grpmax variant allele frequency / upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency in gnomAD is \<0.00009. Do not use data for which the variant does not pass quality control filters. * If only 1 or 2 variant alleles are present in the Grpmax population, use the Grpmax variant allele frequency * If at least 3 variant alleles are present in the Grpmax population, use the upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency Grpmax refers to the gnomAD subpopulation with the highest variant allele frequency. Use large, non-bottlenecked genetic ancestry groups for the Grpmax; avoid using the Amish, Ashkenazi Jewish, European Finnish, and Remaining Individuals groups as well as the genomes-only data for the Middle Eastern group. The upper bound of the 95% CI must be calculated using variant allele numbers and counts from gnomAD. Confidence interval tools, such as Confit-de-MAF ([https://www.genecalculators.net/confit-de-maf.html](https://www.genecalculators.net/confit-de-maf.html)), can be used. Use the gnomAD version with the largest allele number. For larger deletions or duplications that may not be well represented in gnomAD (e.g., single- or multi-exon events), also confirm the variant is not common in gnomAD SVs, gnomAD CNVs, or the Database of Genomic Variants (DGV) ([https://dgv.tcag.ca/dgv/app/home](https://dgv.tcag.ca/dgv/app/home)).

## GN185 v2.0.0 [Released] genes=SGCG
VCEP: ClinGen Limb Girdle Muscular Dystrophy Expert Panel Specifications to the ACMG/A
Strength: **Supporting**
Desc: Apply if the Grpmax variant allele frequency / upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency in gnomAD is \<0.00009. Do not use data for which the variant does not pass quality control filters. * If only 1 or 2 variant alleles are present in the Grpmax population, use the Grpmax variant allele frequency * If at least 3 variant alleles are present in the Grpmax population, use the upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency Grpmax refers to the gnomAD subpopulation with the highest variant allele frequency. Use large, non-bottlenecked genetic ancestry groups for the Grpmax; avoid using the Amish, Ashkenazi Jewish, European Finnish, and Remaining Individuals groups as well as the genomes-only data for the Middle Eastern group. The upper bound of the 95% CI must be calculated using variant allele numbers and counts from gnomAD. Confidence interval tools, such as Confit-de-MAF ([https://www.genecalculators.net/confit-de-maf.html](https://www.genecalculators.net/confit-de-maf.html)), can be used. Use the gnomAD version with the largest allele number. For larger deletions or duplications that may not be well represented in gnomAD (e.g., single- or multi-exon events), also confirm the variant is not common in gnomAD SVs, gnomAD CNVs, or the Database of Genomic Variants (DGV) ([https://dgv.tcag.ca/dgv/app/home](https://dgv.tcag.ca/dgv/app/home)).

## GN186 v2.0.0 [Released] genes=SGCD
VCEP: ClinGen Limb Girdle Muscular Dystrophy Expert Panel Specifications to the ACMG/A
Strength: **Supporting**
Desc: Apply if the Grpmax variant allele frequency / upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency in gnomAD is \<0.00009. Do not use data for which the variant does not pass quality control filters. * If only 1 or 2 variant alleles are present in the Grpmax population, use the Grpmax variant allele frequency * If at least 3 variant alleles are present in the Grpmax population, use the upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency Grpmax refers to the gnomAD subpopulation with the highest variant allele frequency. Use large, non-bottlenecked genetic ancestry groups for the Grpmax; avoid using the Amish, Ashkenazi Jewish, European Finnish, and Remaining Individuals groups as well as the genomes-only data for the Middle Eastern group. The upper bound of the 95% CI must be calculated using variant allele numbers and counts from gnomAD. Confidence interval tools, such as Confit-de-MAF ([https://www.genecalculators.net/confit-de-maf.html](https://www.genecalculators.net/confit-de-maf.html)), can be used. Use the gnomAD version with the largest allele number. For larger deletions or duplications that may not be well represented in gnomAD (e.g., single- or multi-exon events), also confirm the variant is not common in gnomAD SVs, gnomAD CNVs, or the Database of Genomic Variants (DGV) ([https://dgv.tcag.ca/dgv/app/home](https://dgv.tcag.ca/dgv/app/home)).

## GN187 v2.0.0 [Released] genes=CAPN3
VCEP: ClinGen Limb Girdle Muscular Dystrophy Expert Panel Specifications to the ACMG/A
Strength: **Supporting**
Desc: Apply if the Grpmax variant allele frequency / upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency in gnomAD is \<0.0001. Do not use data for which the variant does not pass quality control filters. * If only 1 or 2 variant alleles are present in the Grpmax population, use the Grpmax variant allele frequency * If at least 3 variant alleles are present in the Grpmax population, use the upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency Grpmax refers to the gnomAD subpopulation with the highest variant allele frequency. Use large, non-bottlenecked genetic ancestry groups for the Grpmax; avoid using the Amish, Ashkenazi Jewish, European Finnish, and Remaining Individuals groups as well as the genomes-only data for the Middle Eastern group. The upper bound of the 95% CI must be calculated using variant allele numbers and counts from gnomAD. Confidence interval tools, such as Confit-de-MAF ([https://www.genecalculators.net/confit-de-maf.html](https://www.genecalculators.net/confit-de-maf.html)), can be used. Use the gnomAD version with the largest allele number. For larger deletions or duplications that may not be well represented in gnomAD (e.g., single- or multi-exon events), also confirm the variant is not common in gnomAD SVs, gnomAD CNVs, or the Database of Genomic Variants (DGV) ([https://dgv.tcag.ca/dgv/app/home](https://dgv.tcag.ca/dgv/app/home)).

## GN188 v2.0.0 [Released] genes=ANO5
VCEP: ClinGen Limb Girdle Muscular Dystrophy Expert Panel Specifications to the ACMG/A
Strength: **Supporting**
Desc: Apply if the Grpmax variant allele frequency / upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency in gnomAD is \<0.0001. Do not use data for which the variant does not pass quality control filters. * If only 1 or 2 variant alleles are present in the Grpmax population, use the Grpmax variant allele frequency * If at least 3 variant alleles are present in the Grpmax population, use the upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency Grpmax refers to the gnomAD subpopulation with the highest variant allele frequency. Use large, non-bottlenecked genetic ancestry groups for the Grpmax; avoid using the Amish, Ashkenazi Jewish, European Finnish, and Remaining Individuals groups as well as the genomes-only data for the Middle Eastern group. The upper bound of the 95% CI must be calculated using variant allele numbers and counts from gnomAD. Confidence interval tools, such as Confit-de-MAF ([https://www.genecalculators.net/confit-de-maf.html](https://www.genecalculators.net/confit-de-maf.html)), can be used. Use the gnomAD version with the largest allele number. For larger deletions or duplications that may not be well represented in gnomAD (e.g., single- or multi-exon events), also confirm the variant is not common in gnomAD SVs, gnomAD CNVs, or the Database of Genomic Variants (DGV) ([https://dgv.tcag.ca/dgv/app/home](https://dgv.tcag.ca/dgv/app/home)).

## GN189 v2.0.0 [Released] genes=SGCA
VCEP: ClinGen Limb Girdle Muscular Dystrophy Expert Panel Specifications to the ACMG/A
Strength: **Supporting**
Desc: Apply if the Grpmax variant allele frequency / upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency in gnomAD is \<0.00009. Do not use data for which the variant does not pass quality control filters. * If only 1 or 2 variant alleles are present in the Grpmax population, use the Grpmax variant allele frequency * If at least 3 variant alleles are present in the Grpmax population, use the upper bound of the 95% confidence interval (95% CI) of the Grpmax variant allele frequency Grpmax refers to the gnomAD subpopulation with the highest variant allele frequency. Use large, non-bottlenecked genetic ancestry groups for the Grpmax; avoid using the Amish, Ashkenazi Jewish, European Finnish, and Remaining Individuals groups as well as the genomes-only data for the Middle Eastern group. The upper bound of the 95% CI must be calculated using variant allele numbers and counts from gnomAD. Confidence interval tools, such as Confit-de-MAF ([https://www.genecalculators.net/confit-de-maf.html](https://www.genecalculators.net/confit-de-maf.html)), can be used. Use the gnomAD version with the largest allele number. For larger deletions or duplications that may not be well represented in gnomAD (e.g., single- or multi-exon events), also confirm the variant is not common in gnomAD SVs, gnomAD CNVs, or the Database of Genomic Variants (DGV) ([https://dgv.tcag.ca/dgv/app/home](https://dgv.tcag.ca/dgv/app/home)).

## GN208 v1.0.0 [Released] genes=AIPL1
VCEP: ClinGen Leber Congenital Amaurosis/early onset Retinal Dystrophy Expert Panel Sp
Strength: **Supporting**
Desc: Absent/rare from controls in an ethnically-matched cohort population sample. * Used if the gnomAD total allele frequency is ≤ 4.0 x 10<sup>-4</sup>.

