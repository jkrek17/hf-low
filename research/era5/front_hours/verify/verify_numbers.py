import json, numpy as np, pandas as pd
from scipy import stats
R='/home/claude/hf-low/research/era5/'
fm=pd.read_csv(R+'front_hours/results/front_monthly.csv')
cat=pd.read_csv(R+'hf_history/results/era5_hf_catalog.csv')
tr=pd.read_csv(R+'hf_history/results/era5_hf_catalog_tracks.csv')
ev=cat[(cat.role=='event')&(cat.basin=='pac')&(cat.season.between(2004,2025))]
tr=tr[tr.track.isin(ev.track)].copy()
tmin=tr.loc[tr.groupby('track').msl.idxmin()].set_index('track').time  # first min
tr['tmin']=tr.track.map(tmin)
tr['in']=tr.basin.notna()&(tr.basin.astype(str).str.strip()!='')
tr['season']=tr.track.map(ev.set_index('track').season)
b=tr[tr['in']&(tr.g800>=71.7)&(tr.lat>=40)&(tr.lat<45)&((tr.lon%360)>=160)&((tr.lon%360)<170)].copy()
b['dep']=b.time<=b.tmin
S=list(range(2004,2026))
H=pd.Series({s:6*(b.season==s).sum() for s in S})
Hd=pd.Series({s:6*((b.season==s)&b.dep).sum() for s in S})
Hm=pd.Series({s:6*((b.season==s)&~b.dep).sum() for s in S})
C=pd.Series({s:b[b.season==s].track.nunique() for s in S})
# archive
d=json.load(open('/home/claude/hf-low/docs/data/hf-lows.json'));lf=d['lowFields']
ia,ib=lf.index('basin'),lf.index('season');ifx=lf.index('fixes')
A={s:0 for s in S}
for l in d['lows']:
    if l[ia]!='pac' or l[ib] not in A: continue
    for f in l[ifx]:
        if f[3]=='HF' and 40<=f[1]<45 and 160<=(f[2]%360)<170: A[l[ib]]+=1
A=pd.Series(A)
def fr(col,months_by_offset):
    out={}
    for s in S:
        v=[fm[(fm.year==s+o)&(fm.month==m)][col].iloc[0] for o,m in months_by_offset]
        out[s]=np.mean(v)
    return pd.Series(out)
son=[(0,9),(0,10),(0,11)];djf=[(0,12),(1,1),(1,2)]
F=fr('front',son);FK=fr('front_ke',son);FO=fr('front_oe',son);FD=fr('front',djf)
out=[]
def P(*a):
    s=' '.join(str(x) for x in a);print(s);out.append(s)
y19=list(range(2007,2026))
P('1. hours',y19[0],'-',y19[-1],'mean %.2f sd %.2f min %d max %d meanstorms %.2f'%(H[y19].mean(),H[y19].std(ddof=1),H[y19].min(),H[y19].max(),C[y19].mean()))
P('  hours by season',H[y19].to_dict())
def resid(v,X):
    X=np.column_stack([np.ones(len(v)),X]);return v-X@np.linalg.lstsq(X,v,rcond=None)[0]
def sp(x,y,seasons,cov=None,perm=20000):
    x=np.asarray(x[seasons],float);y=np.asarray(y[seasons],float)
    rx,ry=stats.rankdata(x),stats.rankdata(y)
    if cov is not None:
        rx=resid(rx,cov);ry=resid(ry,cov)
    rho=np.corrcoef(rx,ry)[0,1]
    rng=np.random.default_rng(1);n=0
    for _ in range(perm):
        if abs(np.corrcoef(rx,rng.permutation(ry))[0,1])>=abs(rho)-1e-12:n+=1
    p=(n+1)/(perm+1)
    ps=stats.spearmanr(x,y).pvalue if cov is None else float('nan')
    return rho,p,ps
def rep(name,*a,**k):
    r,p,ps=sp(*a,**k);P('%-40s rho=%+.3f perm_p=%.3f scipy_p=%.3f'%(name,r,p,ps))
rep('P1 FRONT SON vs H 07-25',F,H,y19)
s22=list(range(2004,2026));dm=np.array([[1.0 if s<2007 else 0.0] for s in s22])
rep('S1 22 seasons dummy<2007',F,H,s22,cov=dm)
rep('S2 FRONT_KE',FK,H,y19);rep('S3 FRONT_OE',FO,H,y19);rep('S4 FRONT DJF',FD,H,y19)
rep('S5 deepening',F,Hd,y19);rep('S6 mature/decay',F,Hm,y19);rep('S7 storm count',F,C,y19)
rep('S8 archive, dummy>=2013',F,A,y19,cov=np.array([[1.0 if s>=2013 else 0.0] for s in y19]))
rep('S9 year removed',F,H,y19,cov=np.array([[float(s)] for s in y19]))
P('  archive box HF fixes by season',A[y19].to_dict())
P('  (S8 no-dummy rho for reference: %+.3f)'%stats.spearmanr(F[y19],A[y19])[0])
for n in (19,22):
    P('3. MDE n=%d: %.3f'%(n,np.tanh((stats.norm.ppf(.975)+stats.norm.ppf(.8))/np.sqrt(n-3))))
def ac(x):
    x=np.asarray(x,float);x=x-x.mean();return (x[1:]*x[:-1]).sum()/(x*x).sum()
P('4. lag1 FRONT_SON %.3f  hours %.3f  corr(SON,DJF) %.3f'%(ac(F[y19]),ac(H[y19]),np.corrcoef(F[y19],FD[y19])[0,1]))
open('verify_output.txt','w').write('\n'.join(out)+'\n')
