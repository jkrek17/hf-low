import numpy as np, pandas as pd, json
R='/home/claude/hf-low/research/era5/'
LAGS=[-48,-36,-24,-12,0,12,24]
NB=2000
rng=np.random.default_rng(20261012)
SEASONS=np.arange(2004,2026)
# pre-draw season-count matrix once (shared across tests; 2000 x 22)
CNT=np.stack([np.bincount(rng.integers(0,22,22),minlength=22) for _ in range(NB)]).astype(float)

def strat_diff(num_h,n_h,num_s,n_s):
    # arrays [..., strata]; returns diff along last axis, nan if no common strata
    ok=(n_h>0)&(n_s>0)
    mh=np.where(ok,num_h/np.where(n_h>0,n_h,1),0); ms=np.where(ok,num_s/np.where(n_s>0,n_s,1),0)
    w=np.where(ok,n_h,0)
    sw=w.sum(-1)
    with np.errstate(invalid='ignore',divide='ignore'):
        return (w*(mh-ms)).sum(-1)/sw

def run_test(val,grp,strat,seas):
    """val array, grp 1=HF 0=SF, strat int 0..11, seas int. returns dict"""
    m=~np.isnan(val); val,grp,strat,seas=val[m],grp[m],strat[m],seas[m]
    si=seas-2004
    def agg(g):
        k=g
        S=np.zeros((22,12));N=np.zeros((22,12))
        sel=grp==k
        np.add.at(S,(si[sel],strat[sel]),val[sel]); np.add.at(N,(si[sel],strat[sel]),1)
        return S,N
    Sh,Nh=agg(1);Ss,Ns=agg(0)
    d0=strat_diff(Sh.sum(0),Nh.sum(0),Ss.sum(0),Ns.sum(0))
    reps=strat_diff(CNT@Sh,CNT@Nh,CNT@Ss,CNT@Ns)
    reps=reps[~np.isnan(reps)]
    n=len(reps)
    p=min(1.0,2*min(((reps>0).sum()+1)/(n+1),((reps<0).sum()+1)/(n+1)))
    se=reps.std(ddof=1)
    h=val[grp==1];s=val[grp==0]
    ok=(Nh.sum(0)>0)&(Ns.sum(0)>0); w=np.where(ok,Nh.sum(0),0)
    mh=np.where(ok,Sh.sum(0)/np.where(Nh.sum(0)>0,Nh.sum(0),1),0); ms=np.where(ok,Ss.sum(0)/np.where(Ns.sum(0)>0,Ns.sum(0),1),0)
    hfw=(w*mh).sum()/w.sum(); sfw=(w*ms).sum()/w.sum()
    sd=np.sqrt((h.var(ddof=1)+s.var(ddof=1))/2)
    return dict(hf=hfw,sf=sfw,hf_raw=h.mean(),sf_raw=s.mean(),diff=d0,se=se,p=p,mde80=2.8*se,sd=sd,d_std=d0/sd,mde_std=2.8*se/sd,n_hf=len(h),n_sf=len(s))

def bh(p):
    p=np.asarray(p);n=len(p);o=np.argsort(p);q=np.empty(n)
    r=p[o]*n/(np.arange(n)+1)
    r=np.minimum.accumulate(r[::-1])[::-1]
    q[o]=np.minimum(r,1);return q

def first_lag(qdict):
    # qdict lag->q for lags -48..0
    L=[-48,-36,-24,-12,0];res=99
    for i in range(len(L)-1,-1,-1):
        if qdict[L[i]]<0.05: res=L[i]
        else: break
    return res

out=[]
def P(*a):
    s=' '.join(str(x) for x in a);print(s);out.append(s)

def compare(ref,mine,key_cols,P_):
    """ref/mine dataframes with same keys; returns mismatches"""
    mism=[];matched=0
    mg=ref.merge(mine,on=key_cols,suffixes=('_r','_m'))
    for _,r in mg.iterrows():
        bad=[]
        if r.n_hf_r!=r.n_hf_m or r.n_sf_r!=r.n_sf_m: bad.append(f"n ref({r.n_hf_r},{r.n_sf_r}) mine({r.n_hf_m},{r.n_sf_m})")
        if abs(r.hf_r-r.hf_m)>0.01: bad.append(f"hf ref {r.hf_r:.4f} mine {r.hf_m:.4f}")
        if abs(r.sf_r-r.sf_m)>0.01: bad.append(f"sf ref {r.sf_r:.4f} mine {r.sf_m:.4f}")
        if abs(r.diff_r-r.diff_m)>0.02*r.sd_r: bad.append(f"diff ref {r.diff_r:.4f} mine {r.diff_m:.4f} (sd {r.sd_r:.3f})")
        pr,pm=r.p_r,r.p_m
        okp=abs(pr-pm)<=0.02 or ((pr<0.05)==(pm<0.05) and pr<0.01 and pm<0.01)
        if not okp: bad.append(f"p ref {pr:.4f} mine {pm:.4f}")
        if (r.q_r<0.05)!=(r.q_m<0.05): bad.append(f"q ref {r.q_r:.4f} mine {r.q_m:.4f}")
        if bad: mism.append((tuple(r[k] for k in key_cols),bad))
        else: matched+=1
    return matched,len(mg),mism

# ---------------- Stage A
st=pd.read_csv(R+'hf_vs_storm/results/storms.csv')
fx=pd.read_csv(R+'intensity/results/fixes_2004.csv.gz')
fx['dt']=pd.to_datetime(fx.time.astype(str),format='%Y%m%d%H')
fx=fx[fx.track.isin(st.track)]
P('storms in storms.csv:',len(st),' with fixes:',st.track.isin(fx.track).sum(),' missing from fixes:',st.loc[~st.track.isin(fx.track),'track'].tolist())
cnt=fx.groupby('track').size()
short=cnt[cnt<2].index
P('tracks with <2 00/12 fixes dropped:',len(short))
fx=fx[~fx.track.isin(short)].sort_values(['track','dt'])
idx=fx.groupby('track').g800.idxmax()  # idxmax returns first max in sorted order -> earliest
# ensure earliest tie: idxmax gives first occurrence
anc=fx.loc[idx.values,['track','dt']].rename(columns={'dt':'adt'})
fx=fx.merge(anc,on='track')
fx['lag']=((fx.dt-fx.adt)/pd.Timedelta(hours=1)).astype(int)
fx=fx[fx.lag.isin(LAGS)]
d=fx.drop(columns=['basin','season']).merge(st[['track','grp','basin','season','mon']],on='track')
d['hf']=(d.grp=='HF').astype(int)
rows=[]
for b in ['atl','pac']:
    for v in ['msl','dp12']:
        for lag in LAGS:
            x=d[(d.basin==b)&(d.lag==lag)]
            r=run_test(x[v].values.astype(float),x.hf.values,(x.mon.values-1),x.season.values)
            r.update(sel='all',basin=b,var=v,lag=lag);rows.append(r)
A=pd.DataFrame(rows);A['q']=bh(A.p.values)
A.to_csv('/tmp/verify_lc/stage_a_mine.csv',index=False)
ref=pd.read_csv(R+'hf_lifecycle/results/stage_a_tests.csv');ref=ref[ref.sel=='all']
m,t,mm=compare(ref,A,['basin','var','lag'],P)
P(f'STAGE A: {m}/{t} rows match (28 expected); mismatches:');[P('  ',k,b) for k,b in mm]
for _,r in A.iterrows():
    P('A',r.basin,r['var'],r.lag,f'hf={r.hf:.3f} sf={r.sf:.3f} diff={r["diff"]:.3f} p={r.p:.4f} q={r.q:.4f} mde_sd={r.mde_std:.3f} n={r.n_hf}/{r.n_sf}')
# n check also vs raw ref columns: mde_std
mg=ref.merge(A,on=['basin','var','lag'],suffixes=('_r','_m'))
P('A mde_std max abs diff:',(mg.mde_std_r-mg.mde_std_m).abs().max(),' se max abs diff',(mg.se_r-mg.se_m).abs().max(),' sd max diff',(mg.sd_r-mg.sd_m).abs().max())
# first lag A
jA=json.load(open(R+'hf_lifecycle/results/stage_a_first_lag.json'))
P('STAGE A first persistent lag:')
for b in ['atl','pac']:
    for v in ['msl','dp12']:
        q={r.lag:r.q for _,r in A[(A.basin==b)&(A['var']==v)].iterrows()}
        fl=first_lag(q);rf=jA[f'{b}_{v}']['first_persistent_lag']
        P('  ',b,v,'mine',fl,'ref',rf,'OK' if fl==rf else 'MISMATCH')

# ---------------- Stage B
V=pd.read_csv(R+'hf_lifecycle/results/lc_values.csv.gz')
S=pd.read_csv(R+'hf_lifecycle/results/lc_sample.csv')
P('lc_values rows',len(V),'tracks',V.track.nunique(),'sample tracks',S.track.nunique(),'values tracks in sample',V.track.isin(S.track).all())
# coverage
cov=[]
for (b,g,l),x in V.groupby(['basin','grp','lag']):
    cov.append(dict(basin=b,grp='HF' if g==0 else 'SF',lag=l,n_storms=S[(S.basin==b)&(S.grp==g)].track.nunique(),n_track_lag=len(x),n_snapped=int(x.snapped.sum())))
C=pd.DataFrame(cov);rc=pd.read_csv(R+'hf_lifecycle/results/stage_b_coverage.csv')
mc=rc.merge(C,on=['basin','grp','lag'],suffixes=('_r','_m'))
bad=mc[(mc.n_storms_r!=mc.n_storms_m)|(mc.n_track_lag_r!=mc.n_track_lag_m)|(mc.n_snapped_r!=mc.n_snapped_m)]
P(f'COVERAGE: {len(mc)-len(bad)}/{len(rc)} rows match (ref rows {len(rc)}, mine {len(C)})')
for _,r in bad.iterrows(): P('  MISMATCH',r.to_dict())
Vs=V[V.snapped==True].copy()
Vs['hf']=(Vs.grp==0).astype(int)
rows=[]
for b in ['atl','pac']:
    for v in ['pc','msl_grad','g48_rmax','gmax_r']:
        for lag in LAGS:
            x=Vs[(Vs.basin==b)&(Vs.lag==lag)]
            r=run_test(x[v].values.astype(float),x.hf.values,(x.mon.values-1),x.season.values)
            r.update(basin=b,var=v,lag=lag);rows.append(r)
B=pd.DataFrame(rows);B['q']=bh(B.p.values);B.to_csv('/tmp/verify_lc/stage_b_mine.csv',index=False)
refb=pd.read_csv(R+'hf_lifecycle/results/stage_b_tests.csv')
m,t,mm=compare(refb,B,['basin','var','lag'],P)
P(f'STAGE B: {m}/{t} rows match (56 expected; ref rows {len(refb)}); mismatches:');[P('  ',k,b) for k,b in mm]
for _,r in B.iterrows():
    P('B',r.basin,r['var'],r.lag,f'hf={r.hf:.3f} sf={r.sf:.3f} diff={r["diff"]:.3f} p={r.p:.4f} q={r.q:.4f} mde_sd={r.mde_std:.3f} n={r.n_hf}/{r.n_sf}')
jB=json.load(open(R+'hf_lifecycle/results/stage_b_first_lag.json'))
P('STAGE B first persistent lag:')
for b in ['atl','pac']:
    for v in ['pc','msl_grad','g48_rmax','gmax_r']:
        q={r.lag:r.q for _,r in B[(B.basin==b)&(B['var']==v)].iterrows()}
        fl=first_lag(q);rf=jB[f'{b}_{v}']['first_persistent_lag']
        P('  ',b,v,'mine',fl,'ref',rf,'OK' if fl==rf else 'MISMATCH')

# ---------------- g800 reproduction
allfx=pd.read_csv(R+'intensity/results/fixes_2004.csv.gz',usecols=['track','time','g800'])
V2=V.merge(allfx.rename(columns={'g800':'tab'}),on=['track','time'],how='left')
V2=V2[(V2.time%100).isin([0,12])&(V2.tab>0)&V2.g800.notna()]
dd=(V2.g800-V2.tab).abs()
P(f'G800: rows with 00/12 UTC time, table g800>0, recomputed present: {len(V2)}; within 0.5 kt: {(dd<=0.5).sum()} ({(dd<=0.5).mean()*100:.2f}%); max abs diff {dd.max():.2f}')
for nm,sub in [('lag0',V2.lag==0),('lag!=0',V2.lag!=0),('snapped',V2.snapped==True)]:
    P(f'   {nm}: n={sub.sum()} share={(dd[sub]<=0.5).mean()*100:.2f}%')
V3=V[V.g800.notna()];P('rows with recomputed g800:',len(V3),' of which 00/12 time:',(V3.time%100).isin([0,12]).sum())
open('/tmp/verify_lc/out.txt','w').write('\n'.join(out)+'\n')
