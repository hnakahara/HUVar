#!/usr/bin/env python3
"""
Write the criteria that require manual curation (10 criteria plus conditional BS2) from the original (eRepo) and the VCF
in long format (variant_id / criterion / strength / evidence).

Matching: GENE:HGVSc of the original == ID column of the VCF -> chr:pos:ref:alt of the VCF is used as variant_id.
The strength is the one applied in the eRepo (the default strength of the code when there is no suffix).
"""
import re
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--up",
    default="inputs",
    help="Input directory"
)
parser.add_argument(
    "--out",
    default="outputs",
    help="Output directory"
)

parser.add_argument(
    "--pvs1-rna",
    default="pvs1_rna_list.tsv",
    help="PVS1 based on RNA observation (gene / hgvs_c / strength)"
)
parser.add_argument(
    "--disease-prevalence",
    default="disease_prevalence.tsv",
    help="HUVar resources/shared/disease_prevalence.tsv (BS2 applicability per gene)"
)

args = parser.parse_args()

UP = args.up
OUT = args.out

# criteria supplied by curation (difficult to automate)
# PP1 was added in the 2026-10 revision (automated PP1 from ClinVar text was removed; PP1 is curated only)
TARGET = ['PS2', 'PM3', 'PM6', 'PP4', 'BS3', 'BS4', 'BP2', 'BP5', 'PS4', 'PP1', 'BS2']
# BS2 is added only for genes whose VCEP does not allow BS2 based on gnomAD (bs2 == not_applicable in
# disease_prevalence.tsv). For these genes BS2 rests on internal cohorts or clinical information and
# cannot be automated (2026-10 revision). BS2 for genes that can be evaluated from gnomAD is
# evaluated automatically by HUVar and is therefore not included in the supplement.
_dp = pd.read_csv(args.disease_prevalence, sep='\t', dtype=str)
BS2_MANUAL_GENES = set(_dp.loc[_dp['bs2'].fillna('') == 'not_applicable', 'gene_symbol'])

SUFFIX_MAP = {'Supporting':'Supporting', 'Moderate':'Moderate', 'Strong':'Strong',
              'Very':'VeryStrong', 'Very Strong':'VeryStrong', 'VeryStrong':'VeryStrong',
              # HUVar has no Standalone strength, so "_Stand" codes other than BA1
              # (e.g. BP5_Stand in the eRepo) are treated as VeryStrong
              'Stand':'VeryStrong', 'Standalone':'VeryStrong'}

def default_weight(code):
    if code.startswith('PVS'): return 'VeryStrong'
    if code.startswith('PS'):  return 'Strong'
    if code.startswith('PM'):  return 'Moderate'
    if code.startswith('PP'):  return 'Supporting'
    if code.startswith('BA'):  return 'Standalone'
    if code.startswith('BS'):  return 'Strong'
    if code.startswith('BP'):  return 'Supporting'
    return 'Supporting'

def parse_codes(s):
    """'PM3_Strong,PP4_Moderate,...' → {'PM3':'Strong','PP4':'Moderate',...}"""
    out = {}
    if pd.isna(s): return out
    for tok in str(s).split(','):
        tok = tok.strip()
        if not tok: continue
        if '_' in tok:
            code, suf = tok.split('_', 1)
            code = code.strip(); suf = suf.strip()
            w = SUFFIX_MAP.get(suf, default_weight(code))
        else:
            code = tok; w = default_weight(code)
        out[code] = w
    return out

def norm(s):
    return '' if pd.isna(s) else re.sub(r'\s+', '', str(s))

def chrom_key(vid):
    c = vid.split(':')[0].replace('chr', '')
    special = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
    if c in special: return special[c]
    return int(c) if c.isdigit() else 99

# ---- read and match ----
orig = pd.read_csv(f"{UP}/erepo-original.csv")
vcf  = pd.read_csv(f"{UP}/erepo_hg19.vcf", sep='\t', comment='#', header=None,
                   names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'], dtype=str)

orig['key'] = orig['Gene'].map(norm).str.upper() + '|' + orig['HGVSc'].map(norm)
vcf['vid']  = vcf['CHROM'] + ':' + vcf['POS'] + ':' + vcf['REF'] + ':' + vcf['ALT']
vcf['key']  = (vcf['ID'].str.split(':', n=1).str[0].map(norm).str.upper() + '|'
               + vcf['ID'].str.split(':', n=1).str[1].map(norm))
vcf = vcf.drop_duplicates(subset='key', keep='first')

m = orig.merge(vcf[['key', 'vid']], on='key', how='left')
unmapped = m['vid'].isna().sum()

# ---- expand the target criteria to long format ----
rows = []
for _, r in m.iterrows():
    if pd.isna(r['vid']):
        continue                       # skip variants without coordinates
    codes = parse_codes(r['Applied Evidence Codes (Met)'])
    for code in TARGET:
        if code == 'BS2' and str(r['Gene']).strip().upper() not in BS2_MANUAL_GENES:
            continue
        if code in codes:
            rows.append({'variant_id': r['vid'], 'criterion': code,
                         'strength': codes[code], 'evidence': 'eRepo'})

# PVS1 based on RNA observation (2026-10 revision): RNA / splicing assays cannot be
# evaluated automatically, so eRepo PVS1 applied on RNA evidence (explicit
# "PVS1(RNA)" / "PVS1_RNA", or an RNA experiment cited as the basis of PVS1) is
# supplied as manual evidence at the eRepo strength, irrespective of HUVar's call.
import os
n_rna = 0
if args.pvs1_rna and os.path.exists(args.pvs1_rna):
    rl = pd.read_csv(args.pvs1_rna, sep='\t', dtype=str)
    rmap = dict(zip(rl['gene'].map(norm).str.upper() + '|' + rl['hgvs_c'].map(norm), rl['strength']))
    for _, r in m.iterrows():
        if pd.isna(r['vid']) or r['key'] not in rmap:
            continue
        rows.append({'variant_id': r['vid'], 'criterion': 'PVS1',
                     'strength': rmap[r['key']], 'evidence': 'eRepo (PVS1 RNA)'})
        n_rna += 1
else:
    print(f"[warn] PVS1 RNA list not found: {args.pvs1_rna}")

out = pd.DataFrame(rows, columns=['variant_id', 'criterion', 'strength', 'evidence'])

# remove duplicate rows and sort by coordinate and criterion
out = out.drop_duplicates()
crit_order = {c: i for i, c in enumerate(TARGET + ['PVS1'])}
out['_c'] = out['variant_id'].map(chrom_key)
out['_p'] = out['variant_id'].str.split(':').str[1].astype(int)
out['_o'] = out['criterion'].map(crit_order)
out = out.sort_values(['_c', '_p', '_o']).drop(columns=['_c', '_p', '_o'])

# check variant x criterion pairs with conflicting strengths (from duplicate records)
conflict = out.groupby(['variant_id', 'criterion'])['strength'].nunique()
conflict = conflict[conflict > 1]

out.to_csv(f"{OUT}/erepo_manual_criteria_hg19.tsv", sep='\t', index=False)

print(f"Matched variants: {m['vid'].notna().sum()} / excluded without coordinates: {unmapped}")
print(f"Rows written: {len(out)}  (variant x criterion)")
print(f"Variants with curated evidence: {out['variant_id'].nunique()}")
print("\nRows per criterion:")
print(out['criterion'].value_counts().reindex(TARGET + ['PVS1']))
print(f"PVS1 (RNA) rows: {n_rna}")
print("\nRows per strength:")
print(out['strength'].value_counts())
print(f"\nvariant x criterion pairs with conflicting strengths: {len(conflict)}")
print("\n--- preview ---")
print(out.head(8).to_string(index=False))
