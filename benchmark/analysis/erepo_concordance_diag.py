import sys,re,pandas as pd,collections
from pathlib import Path
d=Path(sys.argv[1]); o=Path(sys.argv[2]); cfg=sys.argv[3]; asm=sys.argv[4] if len(sys.argv)>4 else 'hg38'
pv=pd.read_csv(o/f'concordance_{asm}_{cfg}_per_variant.tsv',sep='\t',dtype=str)
codes=['PVS1','PS1','PS2','PS3','PS4','PM1','PM2','PM3','PM4','PM5','PM6','PP1','PP2','PP3','PP4','BA1','BS1','BS2','BS3','BS4','BP1','BP2','BP3','BP4','BP5','BP7']
r=pd.read_csv(d/f'erepo_{asm}_huhvar_{cfg}_results.tsv',sep='\t',dtype=str,usecols=['variant_id']+codes+[c+'_strength' for c in codes])
m=pv.merge(r,on='variant_id')
def parse(s):
    out={}
    for t in str(s).split(','):
        t=t.strip()
        if not t or t=='nan': continue
        b=re.match(r'([A-Z]+\d)',t)
        if not b: continue
        b=b.group(1)
        sfx=t[len(b):].lstrip('_').lower()
        st={'':None,'very':'VeryStrong','verystrong':'VeryStrong','very strong':'VeryStrong','strong':'Strong','moderate':'Moderate','supporting':'Supporting','stand alone':'StandAlone','standalone':'StandAlone'}.get(sfx,sfx)
        out[b]=st
    return out
m['er']=m['Applied Evidence Codes (Met)'].map(parse)
# Bayesian points for eRepo codes
dflt={'PVS':8,'PS':4,'PM':2,'PP':1,'BA':None,'BS':-4,'BP':-1}
pts={'VeryStrong':8,'Strong':4,'Moderate':2,'Supporting':1}
def score(dd):
    s=0
    for c,st in dd.items():
        if c=='BA1': s-=8; continue  # BA1 = -8 points (scored, no override)
        sign=-1 if c.startswith('B') else 1
        base=dflt[re.match(r'[A-Z]+',c).group(0)]
        p=pts.get(st, abs(base) if base else 0)
        s+=sign*p
    return s
def cls(s):
    if s<=-99: return 'B'
    return 'P' if s>=10 else 'LP' if s>=6 else 'VUS' if s>=0 else 'LB' if s>=-6 else 'B'
m['erepo_codes_cls']=[cls(score(x)) for x in m.er]
rank={'P':2,'LP':1,'VUS':0,'LB':-1,'B':-2}
def cat(t,p):
    if t==p: return 'exact'
    a,b=rank[t],rank[p]
    return 'confidence' if a*b>0 else ('significant' if a==0 or b==0 else 'opposite')
m['codes_cat']=[cat(t,p) for t,p in zip(m.truth,m.erepo_codes_cls)]
print('eRepo codes -> HUVar Bayesian point rules vs eRepo assertion:'); print(m.codes_cat.value_counts(normalize=True).round(3).to_dict())
# criterion-level diffs among discordant (significant/opposite)
dis=m[m['cat'].isin(['significant','opposite'])]
print('discordant n=',len(dis))
fn=collections.Counter(); fp=collections.Counter(); sd=collections.Counter()
for _,x in dis.iterrows():
    for c in codes:
        h=x[c]=='1'; e=c in x.er
        if e and not h: fn[c]+=1
        if h and not e: fp[c]+=1
        if h and e and x.er[c] and x[c+'_strength']!=x.er[c]: sd[c]+=1
t=pd.DataFrame({'missing_in_HUVar':fn,'extra_in_HUVar':fp,'strength_diff':sd}).fillna(0).astype(int)
t['total']=t.sum(1); print(t.sort_values('total',ascending=False).head(15).to_string())
# direction
print(pd.crosstab(dis.truth,dis.pred))
# variants where codes reproduce truth but HUVar discordant
print('of discordant, eRepo codes reproduce assertion exactly/conf:', (dis.codes_cat.isin(['exact','confidence'])).mean().round(3))
t.to_csv(o/f'diag_{cfg}_criteria.tsv',sep='\t')
m[['variant_id','Gene','truth','pred','cat','erepo_codes_cls','codes_cat','Applied Evidence Codes (Met)']].to_csv(o/f'diag_{cfg}_variants.tsv',sep='\t',index=False)
