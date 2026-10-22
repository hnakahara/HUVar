import re, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
U='/root/.claude/uploads/cf4f42ee-70a5-5d68-882c-91f8640dd0a6/'
m=pd.read_csv(U+'b831c9d9-lof_rules_metrics.tsv',sep='\t')
summ=open(U+'921c600b-lof_rules_summary.txt').read()
inh={}
blk=summ.split('by inheritance (%):')[1].split('\n\n')[0].splitlines()
for ln in blk[2:]:
    mm=re.match(r'(.+?)\s+([\d.]+)%\s+([\d.]+)%\s+([\d.]+)%\s+([\d.]+)%\s*$',ln)
    if mm: inh[mm.group(1).strip()]=[float(x)/100 for x in mm.groups()[1:]]
ninh=[int(x) for x in re.search(r'\(n genes\)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',summ).groups()]
def parse(rule):
    loeuf='';frac='';mn=''
    if 'LOEUF_v4' in rule: loeuf='v4.1 < 0.45'
    if 'LOEUF_v2' in rule: loeuf='v2.1.1 < 0.35'
    f=re.search(r'frac>=(\d+)%',rule); frac=int(f.group(1))/100 if f else ''
    p=re.search(r'P/LP>=(\d+)',rule); mn=int(p.group(1)) if p else ''
    if rule.startswith('count>=3'): return 'ClinVar ≥3 P/LP null variants (original rule)','', '', ''
    desc=[]
    if loeuf: desc.append('LOEUF '+loeuf)
    if frac!='': desc.append(f'null ≥{int(frac*100)}% of P/LP (≥{mn} P/LP)')
    return ' OR '.join(desc), loeuf, frac, mn
DSM,LW=2519,538
wb=Workbook(); ws=wb.active; ws.title='LoF rules'
F=Font(name='Arial',size=10); B=Font(name='Arial',size=10,bold=True)
title='Supplementary Table S[X]. Rules for establishing loss of function (LoF) as a disease mechanism (PVS1) for genes without a VCEP PVS1 decision, compared against ClinGen haploinsufficiency (HI) scores and VCEP PVS1 decisions.'
ws['A1']=title; ws['A1'].font=Font(name='Arial',size=11,bold=True)
hdr1=['Rule','Rule components','','','Reference: ClinGen HI score 3 (n = 420) vs HI 0/1/40 (n = 500)','','','','Reference: VCEP PVS1 applicable (n = 100) vs not applicable (n = 31)','','','','Recessive genes (HI 30, n = 597)','ClinGen GDV Definitive/Strong/Moderate genes (n = 2,519)','','Recognised among GDV D/S/M genes, by inheritance','','','','ClinGen GDV Limited or weaker (n = 538)','']
hdr2=['','LOEUF cut-off','Null fraction ≥','Min. P/LP variants','Sensitivity','Specificity','PPV','Balanced accuracy','Sensitivity','Specificity','PPV','Balanced accuracy','Recognised (fraction)','Genes recognised','Fraction',f'AD (n = {ninh[0]})',f'AR (n = {ninh[1]})',f'AD/AR (n = {ninh[2]})',f'XL (n = {ninh[3]})','Genes recognised','Fraction']
ws.append([]); ws.append(hdr1); ws.append(hdr2)
for c in range(1,len(hdr2)+1):
    for r in (3,4):
        cell=ws.cell(r,c); cell.font=B; cell.alignment=Alignment(wrap_text=True,vertical='center',horizontal='center')
        cell.fill=PatternFill('solid',fgColor='D9E1F2')
for a,b in [(2,4),(5,8),(9,12),(14,15),(16,19),(20,21)]:
    ws.merge_cells(start_row=3,start_column=a,end_row=3,end_column=b)
ws.merge_cells('A3:A4')
thin=Side(style='thin',color='999999')
adopted='frac>=20% (P/LP>=5)'
r0=5
for i,row in m.iterrows():
    desc,lo,fr,mn=parse(row.rule)
    vals=[desc,lo or '—',fr if fr!='' else '—',mn if mn!='' else ('—' if not row.rule.startswith('count') else '—'),
          row.hi_sens,row.hi_spec,row.hi_ppv,row.hi_ba,row.vcep_sens,row.vcep_spec,row.vcep_ppv,row.vcep_ba,
          row.hi30_recognised,int(row.gdv_dsm_genes),None,*inh.get(row.rule,['']*4),int(row.gdv_limited_or_weaker_genes),None]
    ws.append(vals); r=ws.max_row
    ws.cell(r,15).value=f'=N{r}/{DSM}'; ws.cell(r,21).value=f'=T{r}/{LW}'
    for c in range(1,len(vals)+1):
        cell=ws.cell(r,c); cell.font=F; cell.border=Border(bottom=thin)
        if c in (3,): cell.number_format='0%'
        elif c in (5,6,7,8,9,10,11,12,13): cell.number_format='0.000'
        elif c in (15,16,17,18,19,21): cell.number_format='0.0%'
        if c>1: cell.alignment=Alignment(horizontal='center')
    if row.rule==adopted:
        for c in range(1,len(vals)+1):
            ws.cell(r,c).font=B; ws.cell(r,c).fill=PatternFill('solid',fgColor='FFF2CC')
last=ws.max_row
notes=['Notes',
 'Rule: genes are recognised when the stated condition is met; "OR" combines LOEUF with the ClinVar null-fraction rule. The original rule (≥3 P/LP null variants) is shown first. The rule adopted in the revised version (null variants ≥20% of ≥5 P/LP variants) is highlighted.',
 'Null variants: nonsense, frameshift and canonical ±1/2 splice-site variants. P/LP: Pathogenic, Likely pathogenic or Pathogenic/Likely pathogenic in ClinVar (release 2026-05-30), review status ≥1 star, all variant types. Null fraction = null P/LP variants / all P/LP variants of the gene; at least one null variant required.',
 'LOEUF: gnomAD v4.1 constraint (cut-off 0.45, gnomAD recommendation for v4.1.1) and gnomAD v2.1.1 constraint (cut-off 0.35), representative transcript per gene (MANE Select, else canonical). The gnomAD v4.1 constraint table has no values for chromosome X genes, hence 0% of X-linked genes with LOEUF v4.1 alone.',
 'References: ClinGen Dosage Sensitivity haploinsufficiency score 3 (sufficient evidence) as positive versus scores 0, 1 and 40 (no evidence, little evidence, dosage sensitivity unlikely) as negative; scores 2 and 30 excluded. Because PVS1 also applies to recessive LoF genes, which may lack evidence for haploinsufficiency, specificity against the HI reference is a lower bound. VCEP: PVS1 designated applicable versus not applicable in the ClinGen CSpec specifications (131 genes).',
 'Recessive genes: genes with HI score 30 (gene associated with an autosomal recessive phenotype). GDV: ClinGen Gene-Disease Validity, best classification per gene; inheritance from the HUVar gene inheritance table (NIH Clinical Genomic Database and curated additions).',
 'In HUVar, genes with a VCEP PVS1 decision or a ClinGen HI score of 3 are evaluated by those routes; the rule applies to the remaining genes. Per-gene data: sheet "Per gene".']
ws.cell(last+2,1,notes[0]).font=B
for k,t in enumerate(notes[1:]):
    c=ws.cell(last+3+k,1,t); c.font=F
ws.column_dimensions['A'].width=46
for c in range(2,len(hdr2)+1): ws.column_dimensions[get_column_letter(c)].width=12
ws.row_dimensions[3].height=48; ws.row_dimensions[4].height=42
ws.freeze_panes='B5'
# VCEP null fraction sheet
w2=wb.create_sheet('Null fraction, VCEP genes')
w2.append(['ClinVar null fraction of genes with a VCEP PVS1 decision']); w2['A1'].font=Font(name='Arial',size=11,bold=True)
w2.append([])
w2.append(['VCEP PVS1','Genes with ClinVar P/LP','Median','IQR lower','IQR upper','< 5%','< 20%'])
for ln in summ.splitlines():
    mm=re.match(r'\s+PVS1 (\S+)\s+n=\s*(\d+)\s+median ([\d.]+)\s+IQR ([\d.]+)-([\d.]+)\s+<5%: (\d+)\s+<20%: (\d+)',ln)
    if mm: g=mm.groups(); w2.append([g[0].replace('_',' '),int(g[1]),float(g[2]),float(g[3]),float(g[4]),int(g[5]),int(g[6])])
for row in w2.iter_rows(min_row=3):
    for c in row: c.font=B if c.row==3 else F
for r in (4,5):
    for c in (3,4,5): w2.cell(r,c).number_format='0.000'
w2.column_dimensions['A'].width=20
for c in 'BCDEFG': w2.column_dimensions[c].width=14
w2.append([]); w2.append(['Median gnomAD v4.1 LOEUF of GDV Definitive/Strong/Moderate genes'])
w2.cell(w2.max_row,1).font=B
w2.append(['Inheritance','Genes with LOEUF','Median LOEUF'])
for c in w2[w2.max_row]: c.font=B
for ln in summ.splitlines():
    mm=re.match(r'median LOEUF v4, (\S+)\s*: ([\d.]+) \(n=(\d+)\)',ln)
    if mm: w2.append([mm.group(1),int(mm.group(3)),float(mm.group(2))])
for row in w2.iter_rows(min_row=10):
    for c in row: c.font=F
# per gene
g=pd.read_csv(U+'a7e815f5-lof_rules_per_gene.tsv',sep='\t',dtype={'hi_score':str})
w3=wb.create_sheet('Per gene')
cols=['gene','vcep_pvs1','hi_score','gdv_best','inheritance','clinvar_plp_null','clinvar_plp_all','null_fraction','loeuf_v4','loeuf_v2']
names=['Gene','VCEP PVS1','ClinGen HI score','ClinGen GDV (best)','Inheritance','ClinVar P/LP null','ClinVar P/LP all','Null fraction','LOEUF gnomAD v4.1','LOEUF gnomAD v2.1.1','LoF established (HUVar revised)','Route']
w3.append(names)
for c in w3[1]: c.font=B; c.fill=PatternFill('solid',fgColor='D9E1F2')
for _,r in g.iterrows():
    hi=str(r.hi_score).strip() if pd.notna(r.hi_score) else ''
    n,t=int(r.clinvar_plp_null),int(r.clinvar_plp_all)
    cv= t>=5 and n>0 and n/t>=0.2
    vc=r.vcep_pvs1 if pd.notna(r.vcep_pvs1) else ''
    route='VCEP' if vc=='applicable' else 'VCEP: not applicable' if vc=='not_applicable' else 'HI score 3' if hi=='3' else 'ClinVar null fraction' if cv else ''
    est='Yes' if route in ('VCEP','HI score 3','ClinVar null fraction') else 'No'
    w3.append([r.gene,vc.replace('_',' '),hi,r.gdv_best if pd.notna(r.gdv_best) else '',r.inheritance if pd.notna(r.inheritance) else '',n,t,
               (n/t if t else None),(r.loeuf_v4 if pd.notna(r.loeuf_v4) else None),(r.loeuf_v2 if pd.notna(r.loeuf_v2) else None),est,route])
for row in w3.iter_rows(min_row=2):
    for c in row: c.font=F
    row[7].number_format='0.000'; row[8].number_format='0.000'; row[9].number_format='0.000'
w3.freeze_panes='B2'; w3.auto_filter.ref=f'A1:L{w3.max_row}'
for i,wd in enumerate([12,16,10,22,10,10,10,10,10,10,14,20],1): w3.column_dimensions[get_column_letter(i)].width=wd
wb.save('out/Supplementary_Table_LoF_rules.xlsx')
print('ok', last, w3.max_row)
