import os
import numpy as np, pandas as pd
d=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','results','scat_footprints.csv.gz'))
f=d[d.cover>=0.40]
a=f.sort_values('cover',ascending=False).drop_duplicates('track').reset_index(drop=True)
print(len(a),a.basin.value_counts().to_dict(),a.season.nunique(),a.sensor.value_counts().to_dict())
h=a[a.idx>=71.7];print('hi',len(h),h.basin.value_counts().to_dict())
def rk(x):return pd.Series(x).rank().values
def sp(x,y):return np.corrcoef(rk(x),rk(y))[0,1]
def auc(s,y):
    r=rk(s);p=y.astype(bool);n1=p.sum();n0=len(y)-n1
    return (r[p].sum()-n1*(n1+1)/2)/(n1*n0)
def p3(a):
    out=[]
    lo=a[a.idx<71.7]
    for b,g in a[a.idx>=71.7].groupby('basin'):
        l=lo[lo.basin==b]; lm=l.minp.values; ls=l.scat_max.values
        for mp in g.minp.values:
            m=np.abs(lm-mp)<=5
            pass
    return out
def p3(a):
    diffs=[]
    for b in ['atl','pac']:
        l=a[(a.basin==b)&(a.idx<71.7)];lm=l.minp.values;ls=l.scat_max.values
        for _,r in a[(a.basin==b)&(a.idx>=71.7)].iterrows():
            m=np.abs(lm-r.minp)<=5
            if m.any():diffs.append(r.scat_max-ls[m].mean())
    return np.mean(diffs),len(diffs)
def stats(a):
    return (sp(a[a.basin=='atl'].era5_max_fp,a[a.basin=='atl'].scat_max),sp(a[a.basin=='pac'].era5_max_fp,a[a.basin=='pac'].scat_max),auc(a.era5_max_fp.values,(a.scat_max>=45).values),p3(a)[0])
print(stats(a),p3(a))
rng=np.random.default_rng(5);seas=a.season.unique();grp={s:a[a.season==s] for s in seas};B=[]
for _ in range(1000):
    c=pd.concat([grp[s] for s in rng.choice(seas,len(seas))],ignore_index=True);B.append(stats(c))
print(np.percentile(np.array(B),[2.5,97.5],axis=0))
bins=pd.cut(a.idx,[-1e9,54,64,71.7,1e9],right=False)
for k,g in a.groupby(bins):print(k,len(g),g.scat_max.median(),(g.scat_max>=45).mean(),(g.scat_max>=50).mean())
e=np.sort(a.era5_max_fp.values);s=np.sort(a.scat_max.values)
for x in(54,64,71.7):print(x,np.quantile(s,(e<x).mean()))
g=f[f.era5_max_fp>=54]
for k,t in g.groupby('sensor'):print(k,len(t),np.exp(np.log(t.scat_max/t.era5_max_fp).mean()))
