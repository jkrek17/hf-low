import numpy as np, pandas as pd, glob, datetime as dt
R='/home/claude/hf-low/research/era5/'
out=[]
def rec(c,r,q,m): out.append((c,r,q,m))
tr=pd.read_csv(R+'hf_history/results/all_tracks.csv.gz')
t=tr[tr.season>=2004]
for b,hq,sq in (('atl',1104,4032),('pac',896,4104)):
    x=t[t.basin==b].gust800_kt
    h=(x>=71.7).sum(); s=((x>=54)&(x<71.7)).sum()
    rec(f'1. {b} HF n',h,hq,h==hq); rec(f'1. {b} SF n',s,sq,s==sq)
S=pd.read_csv(R+'hf_vs_storm/results/storms.csv')
T=pd.read_csv(R+'hf_vs_storm/results/times_storm.csv')
st=pd.concat([pd.read_csv(f) for f in sorted(glob.glob('/tmp/era5work/hf_vs_storm/stats/*.csv'))],ignore_index=True)
rec('2. anchor rows',len(st),2400,len(st)==2400)
vc=st.anchor.value_counts().to_dict(); rec('2. rows per anchor',vc,'800 each',all(v==800 for v in vc.values()) and len(vc)==3)
nt=T.time.nunique(); rec('2. distinct times (times_storm.csv)',nt,1894,nt==1894)
nst=st.time.nunique(); rec('2. distinct times (stats)',nst,1894,nst==1894)
m=st[st.g800_cat>0]; dd=(m.g800-m.g800_cat).abs()
rec('2. max |g800-g800_cat| (n=%d)'%len(m),round(dd.max(),3),'<=0.5',dd.max()<=0.5)
rec('2. n g800_cat>0 with |d|>0.5',int((dd>0.5).sum()),0,(dd>0.5).sum()==0)
# C1 scalars
a=st.merge(S[['track','mon']],on='track',how='left')
sc=pd.read_csv(R+'hf_vs_storm/results/storm/scalars.csv'); sc=sc[sc.cmp=='C1']
for b in ('atl','pac'):
  for v in ['pc','msl_grad','gmax_r','g48_rmax','d2m_500','gust_factor']:
    H=a[(a.basin==b)&(a.anchor=='HF_onset')].dropna(subset=[v])
    Sf=a[(a.basin==b)&(a.anchor=='SF_peak')].dropna(subset=[v])
    hc=H.groupby('mon').size(); hm=H.groupby('mon')[v].mean(); sm=Sf.groupby('mon')[v].mean()
    common=hc.index.intersection(sm.index); w=hc[common]/hc[common].sum()
    sfw=(w*sm[common]).sum(); hfm=H[v].mean()
    q=sc[(sc.basin==b)&(sc.scalar==v)].iloc[0]
    # also weighted HF over common strata for reference
    for nm,r,qq in (('hf',hfm,q.hf),('sf',sfw,q.sf),('diff',hfm-sfw,q['diff'])):
        tol=0.05*abs(qq) if True else 0
        rec(f'3. {b} {v} {nm} (nHF={len(H)},nSF={len(Sf)})',round(r,4),round(qq,4),abs(r-qq)<=tol)
# large
cl=pd.read_csv(R+'hf_vs_storm/results/large/boxtests.csv'); cl=cl[cl.cmp=='C1']
seasons=list(range(2004,2026)); F={}
for s in seasons:
    d=np.load(f'/tmp/era5work/large/season_{s}.npz'); F[s]=d
dates=[];fld={'z500':[],'u250':[]}
for s in seasons:
    d=F[s]; dates+=list(d['dates'])
    for k in fld: fld[k].append(d[k])
for k in fld: fld[k]=np.concatenate(fld[k])
dates=np.array(dates); di={d:i for i,d in enumerate(dates)}
mmdd=np.array([d[5:] for d in dates])
# day list along Aug20-May31
def key(md):
    mo=int(md[:2]); return (mo-8)%12*100+int(md[3:])  # Aug=0..
ul=sorted(set(mmdd),key=key)
ui={m:i for i,m in enumerate(ul)}
lat=-87.1875+5.625*(20+np.arange(12)); lon=5.625*np.arange(64)
anom={}
for k,A in fld.items():
    clim=np.full((len(ul),12,64),np.nan)
    for m,i in ui.items(): clim[i]=np.nanmean(A[mmdd==m],axis=0)
    sm=np.full_like(clim,np.nan)
    for i in range(len(ul)):
        lo=max(0,i-15); hi=min(len(ul),i+16); sm[i]=np.nanmean(clim[lo:hi],axis=0)
    idx=np.array([ui[m] for m in mmdd]); anom[k]=A-sm[idx]
def boxmask(la0,lo0):
    dl=np.abs((lon-lo0+180)%360-180)
    return (np.abs(lat-la0)[:,None]<=15)&(dl[None,:]<=30)
BOX={'atl':(53.5,318.9),'pac':(42.25,171.5)}
def tm(x):
    x=str(int(float(x))); return dt.datetime(int(x[:4]),int(x[4:6]),int(x[6:8]),int(x[8:10]))
def lagdate(t0,k):
    t=t0-dt.timedelta(days=k)
    if t.hour<12: t=t-dt.timedelta(days=1)
    return t.strftime('%Y-%m-%d')
def boxmean(k,i,mask):
    f=anom[k][i]; w=np.cos(np.deg2rad(lat))[:,None]*np.ones((1,64))
    w=np.where(mask&np.isfinite(f),w,0); 
    if w.sum()==0: return np.nan
    return np.nansum(np.where(w>0,f,0)*w)/w.sum()
for b in ('atl','pac'):
    mask=boxmask(*BOX[b]); rows=[]
    for _,r in S[S.basin==b].iterrows():
        if r.grp=='HF': tt=r.onset_time
        elif r.grp=='SF': tt=r.peak_time
        else: continue
        if pd.isna(tt): continue
        t0=tm(tt)
        if t0.month in (6,7,8): continue
        vals={}
        ok=True
        for k in ('z500','u250'):
            bm=[]
            for kk in range(0,11):
                d=lagdate(t0,kk)
                if d not in di: bm.append(np.nan)
                else: bm.append(boxmean(k,di[d],mask))
            bm=np.array(bm)
            vals[k+'_k0']=bm[0]; vals[k+'_L']=bm[4:11].mean() if np.all(np.isfinite(bm[4:11])) else np.nan
        # exclude if any lag missing in 0,4..10
        rows.append(dict(grp=r.grp,mon=r.mon,**vals))
    D=pd.DataFrame(rows).dropna()
    for k in ('z500','u250'):
      for cfg,col in (('L','_L'),('k0','_k0')):
        c=k+col
        H=D[D.grp=='HF']; Sf=D[D.grp=='SF']
        hc=H.groupby('mon').size(); hm=H.groupby('mon')[c].mean(); sm=Sf.groupby('mon')[c].mean()
        cm=hc.index.intersection(sm.index); w=hc[cm]/hc[cm].sum()
        diff=(w*(hm[cm]-sm[cm])).sum(); hf=(w*hm[cm]).sum(); sf=(w*sm[cm]).sum()
        q=cl[(cl.basin==b)&(cl['var']==k)&(cl.cfg==cfg)].iloc[0]
        tol=0.5 if k=='z500' else 0.05
        rec(f'4. {b} {k} {cfg} n={len(H)}/{len(Sf)} (quoted n {q.n_hf}/{q.n_sf})',f'diff {diff:.3f} hf {hf:.3f} sf {sf:.3f}',f'diff {q["diff"]:.3f} hf {q.hf:.3f} sf {q.sf:.3f}',
            abs(diff-q['diff'])<=tol and abs(hf-q.hf)<=tol and abs(sf-q.sf)<=tol)
md='| claim | recomputed | quoted | match? |\n|---|---|---|---|\n'+'\n'.join(f'| {c} | {r} | {q} | {"YES" if m else "NO"} |' for c,r,q,m in out)
open('/tmp/verify/table.md','w').write(md); print(md)
