"""Reviewer 3 Comment 15: eRepo concordance (exact / confidence / significant / opposite), overall and per gene."""
import sys, re, pandas as pd
from pathlib import Path
d=Path(sys.argv[1]); out=Path(sys.argv[2]); out.mkdir(exist_ok=True,parents=True)
asm=sys.argv[3] if len(sys.argv)>3 else 'hg38'
e=pd.read_csv(d/'erepo-original.csv')
e['key']=e['Gene'].astype(str)+':'+e['HGVSc'].astype(str)
amap={'Pathogenic':'P','Likely Pathogenic':'LP','Uncertain Significance':'VUS','Likely Benign':'LB','Benign':'B'}
e['truth']=e['Assertion'].map(amap)
g=e.groupby('key')['truth'].nunique()
e=e[e.key.isin(g[g==1].index)].drop_duplicates('key')
vcf=[l.rstrip('\n').split('\t') for l in open(d/f'erepo_{asm}.vcf') if not l.startswith('#')]
v=pd.DataFrame([(r[0],r[1],r[3],r[4],r[2]) for r in vcf],columns=['chrom','pos','ref','alt','key'])
v['variant_id']=v.chrom+':'+v.pos+':'+v.ref+':'+v.alt
truth=v.merge(e[['key','Gene','truth','Applied Evidence Codes (Met)']],on='key').drop_duplicates('variant_id')
cmap={'Pathogenic':'P','LikelyPathogenic':'LP','VUS':'VUS','UncertainSignificance':'VUS','LikelyBenign':'LB','Benign':'B'}
rank={'P':2,'LP':1,'VUS':0,'LB':-1,'B':-2}
def cat(t,p):
    if t==p: return 'exact'
    a,b=rank[t],rank[p]
    if a*b>0: return 'confidence'
    if a==0 or b==0: return 'significant'
    return 'opposite'
summ=[];genes=[]
for cfg in ['default','revel','default_manual','revel_manual','spliceai','spliceai_manual']:
    f=d/f'erepo_{asm}_huhvar_{cfg}_results.tsv'
    if not f.exists(): continue
    r=pd.read_csv(f,sep='\t',usecols=['variant_id','classification_2015','classification_bayesian'],dtype=str)
    m=truth.merge(r,on='variant_id')
    for fw in ['classification_bayesian','classification_2015']:
        m['pred']=m[fw].map(cmap)
        mm=m.dropna(subset=['pred']).copy()
        mm['cat']=[cat(t,p) for t,p in zip(mm.truth,mm.pred)]
        c=mm['cat'].value_counts()
        row=dict(config=cfg,framework=fw.split('_')[1],n=len(mm))
        for k in ['exact','confidence','significant','opposite']:
            row[k]=int(c.get(k,0)); row[k+'_pct']=round(100*c.get(k,0)/len(mm),1)
        summ.append(row)
        if fw=='classification_bayesian':
            mm.to_csv(out/f'concordance_{asm}_{cfg}_per_variant.tsv',sep='\t',index=False)
            pg=mm.groupby('Gene')['cat'].value_counts().unstack(fill_value=0)
            for k in ['exact','confidence','significant','opposite']:
                if k not in pg: pg[k]=0
            pg['n']=pg[['exact','confidence','significant','opposite']].sum(1)
            pg['exact_pct']=(100*pg.exact/pg.n).round(1)
            pg['exact_or_conf_pct']=(100*(pg.exact+pg.confidence)/pg.n).round(1)
            pg['config']=cfg
            genes.append(pg.reset_index())
S=pd.DataFrame(summ); S.to_csv(out/f'concordance_{asm}_summary.tsv',sep='\t',index=False)
G=pd.concat(genes); G.to_csv(out/f'concordance_{asm}_per_gene.tsv',sep='\t',index=False)
print(S.to_string(index=False))
