import numpy as np, pandas as pd
R='/home/claude/hf-low/research/era5/hf_pattern_phase/results/'
h=pd.read_csv(R+'hf_index.csv'); v=pd.read_csv(R+'pp_values.csv.gz')
VARS=['pc','msl_grad','gmax','wsmax','a_g48','gmax_r','g48_rmax','g48_right','d2m_500','gust_factor','a_hf']
SEAS=np.arange(2004,2026); NB=2000
rng=np.random.default_rng(20261012)
CNT=np.zeros((NB,len(SEAS)))
for b in range(NB):
    CNT[b]=np.bincount(rng.integers(0,len(SEAS),len(SEAS)),minlength=len(SEAS))
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n)
    r=p[o]*n/np.arange(1,n+1); r=np.minimum.accumulate(r[::-1])[::-1]; q[o]=np.minimum(r,1); return q
def pval(d):
    n=len(d); return min(1,2*min(((d>0).sum()+1)/(n+1),((d<0).sum()+1)/(n+1)))
def strat_test(df,var,strat):
    d=df[df[var].notna()&df.tercile.isin(['top','bottom'])]
    ks=sorted(d[strat].unique()); ki={k:i for i,k in enumerate(ks)}
    T=np.zeros((len(SEAS),len(ks),2,2))  # sum,n
    for s,k,t,y in zip(d.season,d[strat],d.tercile,d[var]):
        T[s-2004,ki[k],0 if t=='top' else 1]+=[y,1]
    def stat(C):
        X=np.tensordot(C,T,axes=(-1,0))  # (...,k,g,2)
        nt,nb=X[...,0,1],X[...,1,1]; ok=(nt>0)&(nb>0)
        with np.errstate(all='ignore'):
            dk=X[...,0,0]/nt-X[...,1,0]/nb
        w=np.where(ok,nt,0); return np.where(w.sum(-1)>0,(np.where(ok,dk,0)*w).sum(-1)/w.sum(-1),np.nan)
    obs=stat(np.ones(len(SEAS))); reps=stat(CNT); reps=reps[~np.isnan(reps)]
    tt=d[d.tercile=='top'][var]; bb=d[d.tercile=='bottom'][var]
    se=reps.std(ddof=1); sd=d[var].std(ddof=1)
    Xo=np.tensordot(np.ones(len(SEAS)),T,axes=(0,0)); okk=(Xo[:,0,1]>0)&(Xo[:,1,1]>0); wk=Xo[okk,0,1]
    topw=(Xo[okk,0,0]/Xo[okk,0,1]*wk).sum()/wk.sum()  # top mean reweighted to usable strata; bottom = topw - diff (stratum-weighted bottom mean)
    return dict(var=var,top=topw,bottom=topw-obs,diff=obs,se=se,p=pval(reps),mde80=2.8*se,sd=sd,
                d_std=obs/sd,mde_std=2.8*se/sd,n_top=len(tt),n_bot=len(bb))
def slope_test(df,var):
    d=df[df[var].notna()&df.z.notna()]
    ks=sorted(d.mon.unique()); ki={k:i for i,k in enumerate(ks)}
    T=np.zeros((len(SEAS),len(ks),5)) # n,sz,sy,szy,szz
    for s,k,z,y in zip(d.season,d.mon,d.z,d[var]):
        T[s-2004,ki[k]]+=[1,z,y,z*y,z*z]
    def stat(C):
        X=np.tensordot(C,T,axes=(-1,0)); n=X[...,0]
        with np.errstate(all='ignore'):
            num=X[...,3]-X[...,1]*X[...,2]/n; den=X[...,4]-X[...,1]**2/n
        num=np.where(n>0,num,0);den=np.where(n>0,den,0); return num.sum(-1)/den.sum(-1)
    obs=stat(np.ones(len(SEAS))); reps=stat(CNT)
    se=reps.std(ddof=1); sd=d[var].std(ddof=1)
    return dict(var=var,slope_per_SD=obs,se=se,p=pval(reps),mde80=2.8*se,sd=sd,d_std=obs/sd,mde_std=2.8*se/sd,n=len(d))
# join from hf_index
hh=h[['track','mon','season','z','tercile','onset_lat','in_scope','basin']].rename(columns=lambda c:c if c=='track' else c+'_h')
vv=v.merge(hh,on='track',how='left')
print('join check mismatches: mon',(vv.mon!=vv.mon_h).sum(),'season',(vv.season!=vv.season_h).sum(),'tercile',(vv.tercile!=vv.tercile_h).sum(),
      'z',(abs(vv.z-vv.z_h)>1e-9).sum(),'lat',(abs(vv.onset_lat-vv.onset_lat_h)>1e-9).sum(),'basin',(vv.basin!=vv.basin_h).sum(),'notinscope',(~vv.in_scope_h.astype(bool)).sum())
print('dup (track,anchor):',vv.duplicated(['track','anchor']).sum())
vv['band']=np.clip((vv.onset_lat-20)//5,0,11).astype(int)
out=[]
def run(name,df,basins,strat,csv,slope=False):
    rows=[]
    for b in basins:
        for var in VARS:
            r=slope_test(df[df.basin==b],var) if slope else strat_test(df[df.basin==b],var,strat); r['basin']=b; rows.append(r)
    res=pd.DataFrame(rows); res['q']=bh(res.p.values)
    ref=pd.read_csv(R+csv); m=res.merge(ref,on=['basin','var'],suffixes=('','_ref'))
    cols=(['slope_per_SD'] if slope else ['top','bottom','diff'])+['se','p','sd','mde80','d_std','mde_std','q']+(['n'] if slope else ['n_top','n_bot'])
    print('\n=== %s (%d rows vs %d ref) ==='%(name,len(res),len(ref)))
    bad=0
    for _,r in m.iterrows():
        sdv=r.sd_ref; flags=[]
        nk=['n'] if slope else ['n_top','n_bot']
        for c in nk:
            if r[c]!=r[c+'_ref']: flags.append(c)
        if not slope:
            for c in ['top','bottom']:
                if abs(r[c]-r[c+'_ref'])>max(0.01 if False else 0,1e-4*abs(r[c+'_ref'])) and abs(r[c]-r[c+'_ref'])>1e-6*max(1,abs(r[c+'_ref'])): flags.append(c)
        dc='slope_per_SD' if slope else 'diff'
        if abs(r[dc]-r[dc+'_ref'])>0.02*sdv: flags.append(dc)
        pa,pb=r.p,r.p_ref
        ok=abs(pa-pb)<=0.02 or ((pa<.05)==(pb<.05) and pa<.01 and pb<.01)
        if not ok: flags.append('p')
        if (r.q<.05)!=(r.q_ref<.05): flags.append('q')
        if abs(r.sd-sdv)>1e-6*sdv: flags.append('sd')
        st='OK ' if not flags else 'MISMATCH'+str(flags)
        bad+=bool(flags)
        n=(f"n={int(r.n)}/{int(r.n_ref)}" if slope else f"n={int(r.n_top)}/{int(r.n_top_ref)},{int(r.n_bot)}/{int(r.n_bot_ref)}")
        if slope: vals=f"slope {r[dc]:.5g}/{r[dc+'_ref']:.5g}"
        else: vals=f"top {r.top:.5g}/{r.top_ref:.5g} bot {r.bottom:.5g}/{r.bottom_ref:.5g} diff {r['diff']:.5g}/{r['diff_ref']:.5g}"
        print(f"{st:10s} {r.basin} {r['var']:12s} {n} {vals} se {r.se:.4g}/{r.se_ref:.4g} p {pa:.4f}/{pb:.4f} q {r.q:.4f}/{r.q_ref:.4f} mde {r.mde80:.4g}/{r.mde80_ref:.4g}")
    print('rows compared',len(m),'rows with mismatch',bad)
on=vv[(vv.anchor=='HF_onset')&(vv.match_km<=25)]
print('onset rows total',(vv.anchor=='HF_onset').sum(),'with match<=25',len(on))
run('P primary',on,['atl','pac'],'mon_h','pp_primary.csv')
run('S1 lat',on,['atl','pac'],'band','pp_S1.csv')
run('S2 slope',vv[(vv.anchor=='HF_onset')&(vv.match_km<=25)],['atl','pac'],None,'pp_S2.csv',slope=True)
run('S3 no grn/ice',on[~on.gmax_near_grn_ice],['atl'],'mon_h','pp_S3.csv')
pk=vv[(vv.anchor=='HF_peak')&(vv.match_km<=25)]
print('peak rows',(vv.anchor=='HF_peak').sum(),'<=25',len(pk))
run('S4 peak',pk,['atl','pac'],'mon_h','pp_S4.csv')
# counts
print('\n=== counts ===')
print(h.groupby(['basin','in_scope']).size()); print(h[h.in_scope].groupby(['basin','tercile']).size()); print('out of scope',(~h.in_scope).sum())
for a in ['HF_onset','HF_peak']:
    d=v[(v.anchor==a)&(v.g800_cat>0)]; print(a,'g800 vs g800_cat max abs diff',(d.g800-d.g800_cat).abs().max(),'n>0.5:',((d.g800-d.g800_cat).abs()>0.5).sum(),'n',len(d),'g800_cat<=0:',(v[v.anchor==a].g800_cat<=0).sum())
