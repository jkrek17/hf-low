"""Independent recomputation of RA-10 numbers (verifier implementation)."""
import numpy as np, pandas as pd, sys
D='/home/claude/hf-low/research/era5/enso_kuroshio/data/'
rng=np.random.default_rng(1)
out=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s); out.append(s)

fx=pd.read_csv(D+'fixes.csv.gz'); fx['t']=pd.to_datetime(fx.time.astype(str),format='%Y%m%d%H')
fx=fx.sort_values(['tid','t']).reset_index(drop=True)
nf=fx.groupby('tid').size()
nxt=fx[['tid','t','lat','lon','msl']].copy(); nxt['t']=nxt.t-pd.Timedelta(hours=24)
m=fx.merge(nxt,on=['tid','t'],suffixes=('','2'))
m['dP']=m.msl-m.msl2
m['B']=m.dP*np.sin(np.radians(60))/np.sin(np.radians((m.lat+m.lat2)/2))
m['dl']=((m.lon2-m.lon+180)%360)-180
m['plon']=(m.lon+m.dl/2)%360; m['plat']=(m.lat+m.lat2)/2
m['tm']=m.t+pd.Timedelta(hours=12)
m=m.reset_index(drop=True)
best=m.loc[m.groupby('tid').B.idxmax()].copy()   # idxmax returns first on ties
P('tracks total',fx.tid.nunique(),'with deepening window',len(best))
best=best.merge(fx.groupby('tid').winter.first().rename('w'),left_on='tid',right_index=True)
# lowest pressure position
lp=fx.loc[fx.groupby('tid').msl.idxmin(),['tid','lon','lat']].rename(columns={'lon':'lplon','lat':'lplat'})
best=best.merge(lp,on='tid')
box=(best.plat>=25)&(best.plat<=67)&(best.plon>=120)&(best.plon<=240)
def inwin(df,kind):
    y=df.w
    if kind=='DJF':
        a=pd.to_datetime(y.astype(str)+'-12-03 00:00'); b=pd.to_datetime((y+1).astype(str)+'-03-01 00:00')
    else:
        a=pd.to_datetime((y+1).astype(str)+'-01-01 00:00'); b=pd.to_datetime((y+1).astype(str)+'-03-28 18:00')+pd.Timedelta(minutes=1)
    return (df.tm>=a)&(df.tm<(b if kind=='DJF' else b))
best['djf']=inwin(best,'DJF'); best['jfm']=inwin(best,'JFM')
W=np.arange(1979,2026)
def pop(thr,kind='djf'):
    return best[box&(best.B>=thr)&best[kind]]
bomb=pop(24); rap=pop(12)
cb=bomb.groupby('w').size().reindex(W,fill_value=0); cr=rap.groupby('w').size().reindex(W,fill_value=0)
P('bombs/winter mean,min,max',round(cb.mean(),2),cb.min(),cb.max(),'| rapid mean',round(cr.mean(),2),cr.min(),cr.max(),'| total bombs',cb.sum())
am=pd.read_csv(D+'pipelineA_match.csv'); hf=set(am[am.gust800_kt>=71.7].tid)
hfp=best[box&best.djf&best.tid.isin(hf)]
chf=hfp.groupby('w').size().reindex(W,fill_value=0)
P('HF-reaching (DJF, Pacific box, deepening pt) per winter 2004-2025 mean',round(chf.loc[2004:2025].mean(),2),'min',chf.loc[2004:2025].min(),'max',chf.loc[2004:2025].max(),'| HF tracks total in match',len(hf))

# indices
wi=pd.read_csv(D+'winter_indices.csv').set_index('winter').loc[W]
sb=pd.read_csv(D+'sst_boxes_raw.csv')
for c in ['nino34','emi_a','emi_b','emi_c']:
    clim=sb[(sb.year>=1991)&(sb.year<=2020)].groupby('month')[c].mean()
    sb[c+'_an']=sb[c]-sb.month.map(clim)
sb['emi']=sb.emi_a_an-0.5*sb.emi_b_an-0.5*sb.emi_c_an
def son(c): return pd.Series({y:sb[(sb.year==y)&sb.month.isin([9,10,11])][c].mean() for y in W})
son_emi=son('emi'); son_n34=son('nino34_an')
def djf(c): return pd.Series({y:np.mean([sb[(sb.year==y)&(sb.month==12)][c].iloc[0],sb[(sb.year==y+1)&(sb.month==1)][c].iloc[0],sb[(sb.year==y+1)&(sb.month==2)][c].iloc[0]]) for y in W})
d_emi=djf('emi')
P('check: DJF EMI from raw vs file max abs diff',round(float(np.abs(d_emi.values-wi.emi.values).max()),4))
z=lambda s:(np.asarray(s,float)-np.mean(s))/np.std(s,ddof=1)
Z=dict(emi=z(wi.emi),n34=z(wi.n34),n4mn3=z(wi.n4mn3),pna=z(wi.pna),yr=z(W),semi=z(son_emi),sn34=z(son_n34))

def winter_stat(df,col,fn=np.mean):
    g=df.groupby('w')[col]; return g.apply(lambda v:fn(v.values)).reindex(W), df.groupby('w').size().reindex(W,fill_value=0)

def design(cols,mask=None):
    X=np.column_stack([np.ones(len(W))]+[Z[c] for c in cols]); return X
def wls(y,X,w):
    sw=np.sqrt(w); b,*_=np.linalg.lstsq(X*sw[:,None],y*sw,rcond=None); return b
def fit(y,w,cols,idx):
    X=design(cols)[idx]; return wls(y[idx],X,w[idx])[1]
def tstat(Xt,yt):
    b,*_=np.linalg.lstsq(Xt,yt,rcond=None); r=yt-Xt@b; s2=r@r/(len(yt)-Xt.shape[1])
    return b[1]/np.sqrt(s2*np.linalg.inv(Xt.T@Xt)[1,1])
def fl_p(y,w,cols,idx,nperm,rg,delta_vec=None):
    X=design(cols)[idx]; yy=y[idx]; sw=np.sqrt(w[idx]); Xt=X*sw[:,None]; yt=yy*sw
    Xr=np.delete(Xt,1,axis=1); br,*_=np.linalg.lstsq(Xr,yt,rcond=None); fr=Xr@br; er=yt-fr
    t0=tstat(Xt,yt); cnt=0
    # vectorised
    n=len(yt); perms=np.array([rg.permutation(n) for _ in range(nperm)])
    Ys=fr[None,:]+er[perms]
    XtX=np.linalg.inv(Xt.T@Xt); H=XtX@Xt.T; B=Ys@H.T; R=Ys-B@Xt.T
    s2=(R**2).sum(1)/(n-Xt.shape[1]); t=B[:,1]/np.sqrt(s2*XtX[1,1])
    return (1+np.sum(np.abs(t)>=abs(t0)-1e-12))/(nperm+1),t0
def boot(y,w,cols,idx,nb,rg):
    ii=np.where(idx)[0] if idx.dtype==bool else idx; X=design(cols); r=[]
    for _ in range(nb):
        s=rg.choice(ii,len(ii)); 
        try: r.append(wls(y[s],X[s],w[s])[1])
        except Exception: pass
    return np.percentile(r,[2.5,97.5])

allidx=np.ones(len(W),bool)
def outcome(df,col,fn=np.mean):
    y,n=winter_stat(df,col,fn); return y.values.astype(float),n.values.astype(float)
res={}
def run(name,y,w,cols,idx=allidx):
    ok=idx&(w>0)&np.isfinite(y); b=fit(y,w,cols,ok); res[name]=b; return b
yl,wb=outcome(bomb,'plon'); ya,_=outcome(bomb,'plat')
yl1,wr=outcome(rap,'plon'); ya1,_=outcome(rap,'plat')
ylp,_=outcome(bomb,'lplon')
def bc(df):
    g=df.groupby('w'); return (g.apply(lambda d:np.average(d.plon,weights=d.B))).reindex(W).values.astype(float)
b10=best[box&(best.B>=24)&best.jfm]; y10,w10=outcome(b10,'plon')
b13,_=outcome(bomb,'plon',lambda v:np.percentile(v,75)); b14,_=outcome(bomb,'plon',lambda v:np.percentile(v,25))
yhf,whf=outcome(hfp,'plon')
e=W>=2004
run('P1',yl,wb,['emi','n34','yr']); run('P2',ya,wb,['emi','n34','yr'])
run('S1',yl1,wr,['emi','n34','yr']); run('S2',ya1,wr,['emi','n34','yr'])
run('S3',yl,wb,['semi','sn34','yr']); run('S4',ya,wb,['semi','sn34','yr'])
run('S5',ylp,wb,['emi','n34','yr']); run('S6',yl,wb,['n4mn3','n34','yr'])
run('S7',yl,wb,['emi','n34']); run('S8',yl,wb,['emi','yr'])
run('S9',bc(rap),wr,['emi','n34','yr']); run('S10',y10,w10,['emi','n34','yr'])
run('S11',yhf,whf,['emi','n34'],e); 
run('S13',b13,wb,['emi','n34','yr']); run('S14',b14,wb,['emi','n34','yr'])
run('S15',yl,wb,['emi','n34','yr'],W<=2003); run('S16',yl,wb,['emi','n34','yr'],W>=2004)
P('\n## Point estimates (deg per SD)')
for k,v in res.items(): P(k,round(float(v),3))
# S12
el=wi.n34.values>=0.5; cp=el&(wi.n4.values>wi.n3.values); ep=el&~(wi.n4.values>wi.n3.values)
P('S12 El Nino winters',el.sum(),'CP',cp.sum(),'EP',ep.sum())
P('S12 CP-EP winter-mean (unweighted) %.3f ; bomb-weighted %.3f ; pooled storms %.3f'%(yl[cp].mean()-yl[ep].mean(),
  np.average(yl[cp],weights=wb[cp])-np.average(yl[ep],weights=wb[ep]),
  bomb[bomb.w.isin(W[cp])].plon.mean()-bomb[bomb.w.isin(W[ep])].plon.mean()))
# P1,P2 inference
NP=20000; NB=5000
for nm,y in [('P1',yl),('P2',ya)]:
    p,t=fl_p(y,wb,['emi','n34','yr'],allidx,NP,rng); ci=boot(y,wb,['emi','n34','yr'],allidx,NB,rng)
    p2,_=fl_p(y,wb,['emi','n34','yr'],allidx,NP,np.random.default_rng(99))
    P(nm,'FL perm p=%.4f (rerun other seed %.4f) t=%.3f; boot 95%% CI [%.3f, %.3f]'%(p,p2,t,*ci))
# power
P('\n## Planted-effect power for P1 (null surrogate = reduced-model fit + permuted residuals)')
X=design(['emi','n34','yr']); sw=np.sqrt(wb); Xt=X*sw[:,None]; yt=yl*sw
Xr=np.delete(Xt,1,axis=1); br,*_=np.linalg.lstsq(Xr,yt,rcond=None); fr=Xr@br; er=yt-fr
pg=np.random.default_rng(5); NPL=300
for d in [0,1.5,2.0,2.5,3.0]:
    hit=0
    for _ in range(NPL):
        ys=(fr+er[pg.permutation(len(er))])/sw + d*Z['emi']
        p,_=fl_p(ys,wb,['emi','n34','yr'],allidx,299,pg); hit+=p<0.05
    P('delta',d,'power',round(hit/NPL,3))
# item 5
P('\n## Item 5')
def wsd(v,w): m=np.average(v,weights=w); return np.sqrt(np.average((v-m)**2,weights=w))
P('SD winter-mean bomb lon: unweighted %.3f weighted %.3f | lat unweighted %.3f weighted %.3f'%(yl.std(ddof=1),wsd(yl,wb),ya.std(ddof=1),wsd(ya,wb)))
P('weighted ddof-corrected lon %.3f lat %.3f'%(wsd(yl,wb)*np.sqrt(47/46),wsd(ya,wb)*np.sqrt(47/46)))
P('corr EMI-N34 %.3f ; EMI-PNA %.3f ; EMI-SON EMI %.3f ; EMI-N4mN3 %.3f'%tuple(np.corrcoef(wi.emi,v)[0,1] for v in [wi.n34,wi.pna,son_emi,wi.n4mn3]))
open('/home/claude/hf-low/research/era5/enso_deepening/verify/VERIFICATION.md','w').write('# Verification output (raw, verify_numbers.py)\n\n```\n'+'\n'.join(out)+'\n```\n')
