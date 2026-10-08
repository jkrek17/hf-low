import json, numpy as np, pandas as pd
from scipy import optimize, stats
R0='/home/claude/hf-low/'
fx=pd.read_csv(R0+'research/era5/hf_structure/results/fixes.csv')
ev=pd.read_csv(R0+'research/era5/tele_intensity/results/event_table.csv')
assert ev.track.is_unique
tj=json.load(open(R0+'docs/data/teleconnections.json'))
oni=np.array(tj['oni']['values'],float); pna=np.array([np.nan if v is None else v for v in tj['pna']['values']],float)
fx=fx.merge(ev[['track','ring']],on='track',how='left')
fx['Dp']=fx.ring-fx.msl
t=fx.time.astype(str)
yr=t.str[:4].astype(int); mo=t.str[4:6].astype(int); dy=t.str[6:8].astype(int)
mi=(yr-2001)*12+(mo-5)-1   # month before
fx['oni']=[oni[i] for i in mi]
d0=pd.Timestamp('2001-05-01')
di=((pd.to_datetime(yr.astype(str)+'-'+mo.astype(str)+'-'+dy.astype(str))-d0).dt.days).values
fx['pna']=[np.nanmean(pna[i-10:i-3]) for i in di]
fx['pna_nmiss']=[np.isnan(pna[i-10:i-3]).sum() for i in di]
LA=np.log(71.7)
def tn_fit(X,y,starts=None):
    # Burridge reparam: gamma=beta/s, tau=1/s ; concave
    n,k=X.shape
    def nll(p):
        g=p[:k]; tau=np.exp(p[k]); e=tau*y-X@g; c=tau*LA-X@g
        ll=np.log(tau)+stats.norm.logpdf(e)-stats.norm.logsf(c)
        # grad
        lam=np.exp(stats.norm.logpdf(c)-stats.norm.logsf(c))
        dg=-(-e)[:,None]*0  # placeholder
        # d ll/d e = -e ; d/d c: +lam
        de_dg=-X; dc_dg=-X
        G=(-e)[:,None]*de_dg + lam[:,None]*dc_dg
        gt=(1/tau + (-e)*y + lam*LA)*tau  # wrt log tau
        grad=np.concatenate([G.sum(0),[gt.sum()]])
        return -ll.sum(), -grad
    b0=np.linalg.lstsq(X,y,rcond=None)[0]; s0=np.std(y-X@b0)
    cands=[np.concatenate([b0/s0,[np.log(1/s0)]])]
    for s in (0.05,0.1,0.2,0.4):
        cands.append(np.concatenate([b0/s,[np.log(1/s)]]))
    if starts is not None: cands.insert(0,starts)
    best=None
    for c in cands:
        r=optimize.minimize(nll,c,jac=True,method='L-BFGS-B',options=dict(maxiter=5000,ftol=1e-15,gtol=1e-10))
        r=optimize.minimize(nll,r.x,method='BFGS',jac=True,options=dict(gtol=1e-9,maxiter=5000))
        if best is None or r.fun<best.fun-1e-10: best=r
    return best
def tn_est(X,y,starts=None):
    r=tn_fit(X,y,starts); k=X.shape[1]; tau=np.exp(r.x[k]); return r.x[:k]/tau,1/tau,r
def design(df,dep='lnDp'):
    return np.column_stack([np.ones(len(df)),df.lnDp,df.lnR,df.lnS,df.oni,df.pna])
def prep(b):
    d=fx[(fx.basin==b)&(fx.Dp>=5)].copy()
    dropped=((fx.basin==b)&(fx.Dp<5)).sum()
    d['oni']=(d.oni-d.oni.mean())/d.oni.std(ddof=0); d['pna']=(d.pna-d.pna.mean())/d.pna.std(ddof=0)
    d['lnDp']=np.log(d.Dp); d['R']=np.maximum(d.gmax_r,25); d['lnR']=np.log(d.R); d['lnS']=np.log(np.sin(np.radians(d.lat)))
    d['y']=np.log(d.g800)
    return d,dropped
def centre(X):
    X=X.copy(); m=X[:,1:].mean(0); X[:,1:]-=m; return X,m
def M2(d):
    X,_=centre(design(d)); b,s,r=tn_est(X,d.y.values); return b,s,r
def ols(X,y): return np.linalg.lstsq(X,y,rcond=None)[0]
def M3(d):
    X=np.column_stack([np.ones(len(d)),d.lnDp,d.lnS,d.oni,d.pna]); return ols(X,d.lnR.values)
out=[]; P=lambda s='':(print(s),out.append(s))
rng=np.random.default_rng(12345)
allp={}
for b in ('atl','pac'):
    d,dr=prep(b)
    P(f'==== {b} ===='); P(f'kept fixes {len(d)}  storms {d.track.nunique()}  seasons {d.season.nunique()}  dropped Dp<5: {dr}  fixes w/ any missing PNA day: {(d.pna_nmiss>0).sum()}  g800 min {d.g800.min()}')
    X,m=centre(design(d)); y=d.y.values
    b2,s2,r=tn_est(X,y)
    # intercept correct for uncentred? report slopes only
    P('M2 slopes lnDp lnR lnSinLat ONI PNA: '+' '.join(f'{v:.4f}' for v in b2[1:])+f'   sigma {s2:.4f}  (const {b2[0]:.4f}, loglik {-r.fun:.3f})')
    # profile over sigma (concave in gamma=beta/sigma at fixed sigma), analytic gradient
    ss=np.linspace(0.07,0.5,44); prof=[]
    for s in ss:
        def f(g,s=s):
            e=y/s-X@g; cc=LA/s-X@g
            lam=np.exp(stats.norm.logpdf(cc)-stats.norm.logsf(cc))
            ll=stats.norm.logpdf(e)-stats.norm.logsf(cc)
            return -ll.sum(), -(X.T@(e-lam))
        rr=optimize.minimize(f,ols(X,y)/s,jac=True,method='BFGS',options=dict(gtol=1e-8,maxiter=2000))
        prof.append(-rr.fun-len(y)*np.log(s))
    j=int(np.argmax(prof)); P(f'profile max at sigma grid {ss[j]:.4f} (grid step 0.01), profile loglik {max(prof):.3f} vs joint {-r.fun:.3f}')
    b3=M3(d); P('M3 slopes lnDp lnSinLat ONI PNA: '+' '.join(f'{v:.4f}' for v in b3[1:]))
    bo=ols(design(d),y); P('OLS ln g800 slopes lnDp lnR lnSinLat ONI PNA: '+' '.join(f'{v:.4f}' for v in bo[1:]))
    c=lambda a,bb:np.corrcoef(a,bb)[0,1]
    P('corr lnDp-lnR %.4f  lnDp-lng %.4f  lnR-lng %.4f  lnSin-lng %.4f  lnDp-lnSin %.4f'%(c(d.lnDp,d.lnR),c(d.lnDp,y),c(d.lnR,y),c(d.lnS,y),c(d.lnDp,d.lnS)))
    # PIT KS
    mu=X@b2*0
    mu=design(d)@np.concatenate([[b2[0]-(m*b2[1:]).sum()],b2[1:]])
    z=(y-mu)/s2; a=(LA-mu)/s2; u=(stats.norm.cdf(z)-stats.norm.cdf(a))/(1-stats.norm.cdf(a))
    P('PIT KS distance %.4f'%stats.kstest(u,'uniform').statistic)
    # bootstrap
    seas=np.array(sorted(d.season.unique())); groups={s:d[d.season==s] for s in seas}
    est=[]; st0=np.concatenate([b2/ s2*0 + 0,[0]])
    for it in range(500):
        pick=rng.choice(seas,len(seas),replace=True)
        db=pd.concat([groups[s] for s in pick])
        Xb,_=centre(design(db)); bb,sb,_=tn_est(Xb,db.y.values)
        m3=M3(db); est.append(np.concatenate([bb[1:4],bb[4:6],m3[3:5],[sb]]))
    est=np.array(est); se=est.std(0,ddof=1)
    names=['M2 lnDp','M2 lnR','M2 lnSinLat','M2 ONI','M2 PNA','M3 ONI','M3 PNA']
    point=np.concatenate([b2[1:],b3[3:5]])
    p=2*stats.norm.sf(np.abs(point/se[:7]))
    P('bootstrap (500 season-block draws, 22 seasons): name est se z p'); 
    for n_,e_,s_,p_ in zip(names,point,se[:7],p): P(f'  {n_:12s} est {e_:8.4f} se {s_:.4f} z {e_/s_:7.3f} p {p_:.4f}')
    P(f'  boot se of sigma {se[7]:.4f}')
    allp[b]=(names,p)
# BH within 14 tests
pa=np.concatenate([allp['atl'][1],allp['pac'][1]]); lab=[f'{b} {n}' for b in ('atl','pac') for n in allp[b][0]]
o=np.argsort(pa); q=np.empty(14); m=14; prev=1
for rank in range(m,0,-1):
    i=o[rank-1]; prev=min(prev,pa[i]*m/rank); q[i]=prev
P('==== BH q (family of 14) ====')
for l,pp,qq in zip(lab,pa,q): P(f'  {l:18s} p {pp:.4f} q {qq:.4f}')
open(R0+'research/era5/gust_scaling/verify/verify_output.txt','w').write('\n'.join(out)+'\n')
