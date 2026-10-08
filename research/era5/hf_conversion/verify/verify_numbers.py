import numpy as np, pandas as pd
R='/home/claude/hf-low/research/era5/'
H=R+'hf_conversion/results/'
fx=pd.read_csv(R+'intensity/results/fixes_2004.csv.gz')
env=pd.read_csv(R+'intensity/results/env_2004.csv.gz')
pairs=pd.read_csv(H+'pairs.csv')
print("fixes",len(fx),"dtypes hf24",fx.hf24.dtype, fx.hf24.unique()[:5])
el=fx[(fx.g800<71.7)&(fx.ndr24>=1.0)&fx.hf24.notna()].copy()
el['hf24']=el.hf24.astype(str).str.lower().isin(['true','1','1.0'])
print("eligible fixes",len(el),"storms",el.track.nunique())
cv=el[el.hf24]; print("converting fixes",len(cv),"storms",cv.track.nunique())
for b,g in el.groupby('basin'): print(b,"elig fixes",len(g),"storms",g.track.nunique(),"conv fixes",g.hf24.sum(),"conv storms",g[g.hf24].track.nunique())
# pairs checks
print("\nPAIRS rows",len(pairs))
m=pairs.merge(el[['track','time','hf24','ndr24','g800','basin','lat','msl','age','speed','dp12']],on=['track','time'],how='left',suffixes=('','_f'))
print("members not eligible:",m.hf24.isna().sum())
print("label mismatch:",((m.hf24==True)&(m.grp==1)).sum()+((m.hf24==False)&(m.grp==0)).sum())
print("ndr24 mismatch vs fixes:",(abs(m.ndr24-m.ndr24_f)>1e-6).sum(), "basin mism",(m.basin!=m.basin_f).sum(),"lat mism",(abs(m.lat-m.lat_f)>1e-6).sum(), "msl mism",(abs(m.msl-m.msl_f)>1e-6).sum())
print("dup (track,grp):",pairs.duplicated(['track','grp']).sum(), "dup (track,time):",pairs.duplicated(['track','time']).sum())
print("grp counts",pairs.grp.value_counts().to_dict(), "pair sizes", pairs.groupby('pair').size().value_counts().to_dict())
w=pairs.pivot(index='pair',columns='grp',values=['basin','mon','ndr24','track','time'])
print("basin differs",(w['basin'][0]!=w['basin'][1]).sum(),"mon differs",(w['mon'][0]!=w['mon'][1]).sum(),
 "max|dndr|",(w['ndr24'][0]-w['ndr24'][1]).abs().max(),"n>=0.25",((w['ndr24'][0]-w['ndr24'][1]).abs()>=0.25).sum())
print("mon vs time month mismatch",(pairs.mon!=(pairs.time//10000%100)).sum(), "season check", (pairs.season!=np.where(pairs.mon>=6,pairs.time//1000000,pairs.time//1000000-1)).sum())
print("pairs per basin",pairs[pairs.grp==0].basin.value_counts().to_dict())
print("noncv reused (time,track):",pairs[pairs.grp==1].duplicated(['track','time']).sum())
print("tracks in both groups",len(set(pairs[pairs.grp==0].track)&set(pairs[pairs.grp==1].track)))
print("unique storms per grp",pairs.groupby('grp').track.nunique().to_dict())

# env join
E=env.set_index(['track','time'])
fx['t']=pd.to_datetime(fx.time.astype(str),format='%Y%m%d%H')
lagmap={}
key=lambda tr,t:(tr,t)
fxi=fx.set_index(['track','time'])
def shift(df,h):
    t=pd.to_datetime(df.time.astype(str),format='%Y%m%d%H')-pd.Timedelta(hours=h)
    return t.dt.strftime('%Y%m%d%H').astype(np.int64)
def lagged_env(df,h,var):
    tt=shift(df,h)
    idx=pd.MultiIndex.from_arrays([df.track.values,tt.values])
    return env.set_index(['track','time'])[var].reindex(idx).values
# note: a lag fix 'present' = row exists in env; value NaN if missing
VARS=['jet250','eady','sstgrad','sst_t500','tcwv','flux','vadv500','div300']
rng=np.random.default_rng(1)
SEASONS=np.arange(2004,2026)
def strat(val,grp,strata,seas,nboot=2000,seed=7):
    d=pd.DataFrame(dict(v=val,g=grp,s=strata,y=seas)).dropna(subset=['v'])
    if len(d)==0: return None
    cv=d[d.g==1]; nc=d[d.g==0]
    skeys=sorted(d.s.unique()); si={k:i for i,k in enumerate(skeys)}
    ys=SEASONS; yi={y:i for i,y in enumerate(ys)}
    A=np.zeros((4,len(ys),len(skeys)))  # nconv,sumconv,nnon,sumnon  (grp 0 = converter!)
    for _,r in d.iterrows():
        pass
    return d,skeys,si,yi
def stratdiff_tables(d,skeys,si,yi):
    A=np.zeros((4,len(yi),len(skeys)))
    for conv,k in ((0,0),(1,2)):
        g=d[d.g==conv]
        np.add.at(A[k],(g.y.map(yi).values,g.s.map(si).values),1)
        np.add.at(A[k+1],(g.y.map(yi).values,g.s.map(si).values),g.v.values)
    return A
def stat_from(T):  # T shape (...,4,S) totals
    nc,sc,nn,sn=T[...,0,:],T[...,1,:],T[...,2,:],T[...,3,:]
    ok=(nc>0)&(nn>0)
    with np.errstate(all='ignore'):
        dif=np.where(ok,sc/nc-sn/nn,0.0)
    w=np.where(ok,nc,0.0)
    ws=w.sum(-1)
    with np.errstate(all='ignore'):
        return np.where(ws>0,(w*dif).sum(-1)/ws,np.nan)
import os
NB=int(os.environ.get('NB','2000'))
DD=0
def test(val,grp,strata,seas,nboot=None,seed=11):
    nboot=nboot or NB
    # grp: 0=converter 1=non
    val=np.asarray(val,float)
    d=pd.DataFrame(dict(v=val,g=grp,s=strata,y=seas)).dropna(subset=['v'])
    out=dict(n_conv=int((d.g==0).sum()),n_non=int((d.g==1).sum()))
    if out['n_conv']<2 or out['n_non']<2: return out|dict(diff=np.nan,p=np.nan)
    skeys=sorted(d.s.unique()); si={k:i for i,k in enumerate(skeys)}; yi={y:i for i,y in enumerate(SEASONS)}
    A=stratdiff_tables(d,skeys,si,yi)  # (4,Y,S)
    A=np.transpose(A,(1,0,2))  # (Y,4,S)
    obs=stat_from(A.sum(0))
    r=np.random.default_rng(seed)
    cnt=r.multinomial(len(SEASONS),np.ones(len(SEASONS))/len(SEASONS),size=nboot)  # (B,Y)
    T=np.einsum('by,yks->bks',cnt,A)
    bs=stat_from(T); bs=bs[~np.isnan(bs)]
    n=len(bs)
    p=min(1,2*min(((bs>0).sum()+1)/(n+1),((bs<0).sum()+1)/(n+1)))
    se=bs.std(ddof=1)
    c=d[d.g==0].v; nn=d[d.g==1].v
    TT=A.sum(0); okS=(TT[0]>0)&(TT[2]>0); cres=TT[1][okS].sum()/TT[0][okS].sum()
    sd=np.sqrt((c.var(ddof=DD)+nn.var(ddof=DD))/2)
    return out|dict(conv=cres,non=cres-obs,conv_raw=c.mean(),non_raw=nn.mean(),diff=obs,p=p,se=se,mde=2.8*se,sd=sd,d_std=obs/sd,mde_std=2.8*se/sd,lo=np.percentile(bs,2.5),hi=np.percentile(bs,97.5),nb=n)
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)
    r=p[o]*n/(np.arange(n)+1); r=np.minimum.accumulate(r[::-1])[::-1]; q[o]=np.minimum(r,1); return q
def cmp(mine,ref,keys,name,cols=('n_conv','n_non')):
    print(f"\n=== {name}: mine rows {len(mine)}, ref rows {len(ref)}")
    mm=mine.merge(ref,on=keys,suffixes=('_m','_r'),how='outer',indicator=True)
    print("unmatched keys:",(mm._merge!='both').sum())
    bad=0;ok=0
    for _,r in mm[mm._merge=='both'].iterrows():
        msgs=[]
        for c in cols:
            if r[c+'_m']!=r[c+'_r']: msgs.append(f"{c} {r[c+'_m']} vs {r[c+'_r']}")
        if 'conv_r' in r:
            if not(abs(r.conv_m-r.conv_r)<0.01 or (np.isnan(r.conv_m) and np.isnan(r.conv_r))): msgs.append(f"conv {r.conv_m:.4f} vs {r.conv_r:.4f}")
            if not(abs(r.non_m-r.non_r)<0.01 or (np.isnan(r.non_m) and np.isnan(r.non_r))): msgs.append(f"non {r.non_m:.4f} vs {r.non_r:.4f}")
            sd=r.sd_r
            if not abs(r.diff_m-r.diff_r)<=0.02*sd: msgs.append(f"diff {r.diff_m:.5g} vs {r.diff_r:.5g} (sd {sd:.4g}, dev {abs(r.diff_m-r.diff_r)/sd:.3f} SD)")
            if abs(r.sd_m-r.sd_r)>1e-3*max(1,abs(r.sd_r)): msgs.append(f"sd {r.sd_m:.5g} vs {r.sd_r:.5g}")
            pm,pr=r.p_m,r.p_r
            if not(abs(pm-pr)<=0.02 or ((pm>=.05)==(pr>=.05) and pm<.01 and pr<.01) or ((pm<.05)==(pr<.05) and pm<.01 and pr<.01)): msgs.append(f"p {pm:.4f} vs {pr:.4f}")
            if abs(r.mde_std_m-r.mde_std_r)>0.02: msgs.append(f"mde_std {r.mde_std_m:.4f} vs {r.mde_std_r:.4f}")
            if (r.q_m<0.05)!=(r.q_r<0.05): msgs.append(f"q-pass {r.q_m:.4f} vs {r.q_r:.4f}")
        if msgs: bad+=1; print("MISMATCH",[r[k] for k in keys],msgs)
        else: ok+=1
    print(f"matched {ok}, mismatched {bad}")
    if 'q_m' in mm:
        x=mm[mm._merge=='both']
        print(" max|q diff|=%.4f max|p diff|=%.4f max|diff dev| in SD=%.4f"%((x.q_m-x.q_r).abs().max(),(x.p_m-x.p_r).abs().max(),((x.diff_m-x.diff_r).abs()/x.sd_r).max()))
        print(" max|se rel diff|=%.3f"%(((x.se_m-x.se_r).abs()/x.se_r).max()))

# ---------- Stage A
rows=[];chg=[];lat=[];
P=pairs.copy()
P['latband']=np.clip((P.lat-20)//5,0,10)
lagH={'t0':0,'t0-12h':12,'t0-24h':24}
vals={}
for v in VARS:
    for ln,h in lagH.items():
        vals[(v,ln)]=lagged_env(P,h,v) if h else E[v].reindex(pd.MultiIndex.from_arrays([P.track.values,P.time.values])).values
# existence of env row at t0
print("\nt0 pair members missing env row:",int(pd.Series(E[VARS[0]].reindex(pd.MultiIndex.from_arrays([P.track.values,P.time.values])).values).isna().sum()),'(NaN jet250 incl.)', 'env index dup:',E.index.duplicated().sum())
for b in ['atl','pac']:
    sel=(P.basin==b).values
    for v in VARS:
        for ln in lagH:
            x=vals[(v,ln)][sel]
            t=test(x,P.grp.values[sel],P.mon.values[sel],P.season.values[sel])
            rows.append(dict(basin=b,var=v,lag=ln)|t)
            tl=test(x,P.grp.values[sel],P.latband.values[sel],P.season.values[sel])
            lat.append(dict(basin=b,var=v,lag=ln)|tl)
        d=vals[(v,'t0')][sel]-vals[(v,'t0-24h')][sel]
        t=test(d,P.grp.values[sel],P.mon.values[sel],P.season.values[sel])
        chg.append(dict(basin=b,var=v,lag='t0 minus t0-24h')|t)
def finish(rows):
    df=pd.DataFrame(rows); df['q']=bh(df.p.values); return df
A=finish(rows); L=finish(lat); C=finish(chg)
cmp(A,pd.read_csv(H+'stage_a_tests.csv'),['basin','var','lag'],'stage A primary (48)')
cmp(C,pd.read_csv(H+'stage_a_change.csv'),['basin','var','lag'],'stage A change (16)')
cmp(L,pd.read_csv(H+'stage_a_latmatched.csv'),['basin','var','lag'],'stage A latmatched (48)')
# confounds
fxi2=fx.set_index(['track','time'])
cf=[]
for b in ['atl','pac']:
    sel=(P.basin==b).values
    for v in ['lat','age','msl','speed','dp12']:
        x=fxi2[v].reindex(pd.MultiIndex.from_arrays([P.track.values[sel],P.time.values[sel]])).values
        cf.append(dict(basin=b,var=v,lag='t0')|test(x,P.grp.values[sel],P.mon.values[sel],P.season.values[sel]))
Cf=pd.DataFrame(cf); Cf['q']=bh(Cf.p.values)  # BH over the 10 descriptive tests (assumed)
cmp(Cf,pd.read_csv(H+'stage_a_confounds.csv'),['basin','var','lag'],'confounds (10)')
# also q within confounds if ref uses a different BH
ref=pd.read_csv(H+'stage_a_confounds.csv'); print(" ref confound q vs p:\n",ref[['basin','var','p','q']].to_string()); print(" my q\n",Cf[['basin','var','p','q']].to_string())

# ---------- Stage B
cs=pd.read_csv(H+'conv_sample.csv'); cvv=pd.read_csv(H+'conv_values.csv.gz')
print("\nconv_sample rows",len(cs),"conv_values rows",len(cvv),"lag values",sorted(cvv.lag.unique()))
# sample pairs a subset of pairs.csv?
mm=cs.merge(pairs,on=['pair','grp','track','time'],how='left',indicator=True); print("conv_sample rows not in pairs.csv:",(mm._merge!='both').sum(),"pairs per basin/grp",cs.groupby(['basin','grp']).size().to_dict(), "pairs complete:",cs.groupby('pair').size().value_counts().to_dict())
print("conv_values lag0 rows == conv_sample rows? ",len(cvv[cvv.lag==0]),len(cs))
DD=1  # stage B reference SD matches ddof=1
cov=[];rowsB=[]
SC=['jet_dist','jet_cos','thetae850','stab','baroc']
for b in ['atl','pac']:
    for g in (0,1):
        for lg in (0,12,24):
            s=cvv[(cvv.basin==b)&(cvv.grp==g)&(cvv.lag==lg)]
            u=s[s.snapped.astype(str).str.lower()=='true']
            if lg==0: u=u[u.match_km<=25]
            cov.append(dict(basin=b,grp=g,lag=lg,n=len(s),usable=len(u)))
covd=pd.DataFrame(cov); refc=pd.read_csv(H+'stage_b_coverage.csv')
mc=covd.merge(refc,on=['basin','grp','lag'],suffixes=('_m','_r'))
print("\n=== coverage rows",len(mc),"mismatches:",((mc.n_m!=mc.n_r)|(mc.usable_m!=mc.usable_r)).sum())
print(mc[(mc.n_m!=mc.n_r)|(mc.usable_m!=mc.usable_r)])
print("lag0 snapped false / match_km>25 counts:",((cvv.lag==0)&(cvv.snapped.astype(str).str.lower()!='true')).sum(),((cvv.lag==0)&(cvv.match_km>25)).sum())
for b in ['atl','pac']:
    for v in SC:
        for lg in (0,12,24):
            s=cvv[(cvv.basin==b)&(cvv.lag==lg)&(cvv.snapped.astype(str).str.lower()=='true')]
            if lg==0: s=s[s.match_km<=25]
            t=test(s[v].values,s.grp.values,s.mon.values,s.season.values)
            rowsB.append(dict(basin=b,var=v,lag=-lg)|t)
B=finish(rowsB)
cmp(B,pd.read_csv(H+'stage_b_tests.csv'),['basin','var','lag'],'stage B (30)')

# ---- noise diagnostic: z of (p_ref - p_mine) given 2000-rep MC error
import sys
print("\n### MC noise diagnostic and borderline q rows")
for name,mine,fn,keys in [('A',A,'stage_a_tests.csv',['basin','var','lag']),('chg',C,'stage_a_change.csv',['basin','var','lag']),('lat',L,'stage_a_latmatched.csv',['basin','var','lag']),('B',B,'stage_b_tests.csv',['basin','var','lag'])]:
    r=pd.read_csv(H+fn); m=mine.merge(r,on=keys,suffixes=('_m','_r'))
    t=np.minimum(m.p_m/2,1-m.p_m/2).clip(0.0005,0.5)
    sdp=2*np.sqrt(t*(1-t)/2000)
    z=(m.p_r-m.p_m)/sdp
    print(name,"n",len(m),"|z|>2:",(z.abs()>2).sum(),"|z|>3:",(z.abs()>3).sum(),"max|z|=%.2f"%z.abs().max(), "SE ratio ref/mine mean %.3f sd %.3f"%((m.se_r/m.se_m).mean(),(m.se_r/m.se_m).std()))
    bl=m[((m.q_m>0.03)&(m.q_m<0.08))|((m.q_r>0.03)&(m.q_r<0.08))]
    print(bl[keys+['diff_m','diff_r','p_m','p_r','q_m','q_r','se_m','se_r']].to_string())
