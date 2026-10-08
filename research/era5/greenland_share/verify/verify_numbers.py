import sys, numpy as np, pandas as pd, statsmodels.api as sm
sys.path.insert(0,'/home/claude/hf-low/research/era5/freq_split')
import split as S
R='/home/claude/hf-low/'
T=S.load_tracks(R+'research/era5/hf_history/results/all_tracks.csv.gz'); T=T[T.basin=='atl'].reset_index(drop=True)
I=S.indices('/mnt/project-files/teleconnection-test/cpc_indices',R)
gh=pd.read_csv(R+'research/era5/nao_share_barrier/gh_daily.csv',parse_dates=['date']).set_index('date').GH.asfreq('D')
I['GH_lag']=gh.rolling(7,min_periods=7).mean().shift(4)
D=S.window_days('octapr',2004,2025)
raw=I.reindex(D.date)[['NAO_lag','GH_lag']].reset_index(drop=True)
ok=raw.notna().all(axis=1).values
D,raw=D[ok].reset_index(drop=True),raw[ok].reset_index(drop=True)
mu,sd=raw.mean(),raw.std(ddof=0); z=(raw-mu)/sd
Tw=S.assign(T,'octapr',2004,2025)
zt=(I.reindex(Tw.gen)[['NAO_lag','GH_lag']].reset_index(drop=True)-mu)/sd
okt=zt.notna().all(axis=1).values
Tw=Tw[okt].copy(); zt=zt[okt]; Tw['NAO']=zt.NAO_lag.values; Tw['GH']=zt.GH_lag.values
print('1: tracks',len(Tw),'HF',Tw.hf.sum())
fx=pd.read_csv(R+'research/era5/intensity/results/fixes_2004.csv.gz'); fx=fx[fx.basin=='atl']
env=pd.read_csv(R+'research/era5/intensity/results/env_2004.csv.gz')
print('  with fix rows',Tw.track.isin(set(fx.track)).sum())
print('2: daily r',np.corrcoef(z.NAO_lag,z.GH_lag)[0,1],'storm r',np.corrcoef(Tw.NAO,Tw.GH)[0,1])
sm_=z.groupby(S.window_days('octapr',2004,2025).season.values[ok]).mean()
print('  season-mean r',np.corrcoef(sm_.NAO_lag,sm_.GH_lag)[0,1])
C=S.daily_counts(Tw,D,'atl',['all','hf'])
def share(cols,C=C,X=None):
    X=S.design(C,z,cols) if X is None else X
    b={o:S.fast_pois(C[o].values.astype(float),np.asarray(X,float)) for o in ('all','hf')}
    return b
for cols in (['GH_lag'],['NAO_lag'],['GH_lag','NAO_lag']):
    b=share(cols)
    print('3/4',cols,{c:(round(np.exp(b['hf'][i]-b['all'][i]),3),round(np.exp(b['hf'][i]),3),round(np.exp(b['all'][i]),3)) for i,c in enumerate(cols)})
# 5
Tw['gt']=pd.qcut(Tw.GH,3,labels=[1,2,3]); Tw['nt']=pd.qcut(Tw.NAO,3,labels=[1,2,3])
print('5',Tw.groupby('gt').hf.agg(['mean','size']).round(3).T.to_string())
print(Tw.groupby(['nt','gt']).hf.agg(['mean','size']).round(3).unstack().to_string())
# 6
fx=fx.sort_values(['track','time']); first=fx.groupby('track').first()
Tw=Tw.merge(first[['lat','lon']].rename(columns={'lat':'flat','lon':'flon'}),left_on='track',right_index=True,how='inner')
Tw['flon']=np.where(Tw.flon>180,Tw.flon-360,Tw.flon)
env['key']=env.track.astype(str)+'_'+env.time.astype(str); 
fe=first.reset_index()[['track','time']]
e=fe.merge(env[['track','time','eady']],on=['track','time'],how='left')
Tw=Tw.merge(e[['track','eady']],on='track',how='left')
print('storms with fix rows',len(Tw),'eady na',Tw.eady.isna().sum())
def covs(T):
    X=pd.DataFrame({'GH':T.GH.values,'NAO':T.NAO.values},index=T.index)
    for mo in sorted(T.gen.dt.month.unique())[1:]: X[f'm{mo}']=(T.gen.dt.month.values==mo).astype(float)
    X['trend']=(T.season.values-T.season.mean())/10.; X['const']=1.
    return X
X=covs(Tw)
for nm in ('flon','minp','eady'):
    m=Tw[nm].notna()
    print('6',nm,round(sm.OLS(Tw[nm][m],X[m]).fit().params['GH'],3))
def lg(Xm): 
    return sm.GLM(Tw.hf.values,Xm,family=sm.families.Binomial()).fit().params['GH']
print('7 base',lg(X))
zs=lambda s:(s-s.mean())/s.std(ddof=0)
X2=X.copy()
for a,b_ in (('flat','flat'),('flon','flon')): X2[a]=zs(Tw[a]).values
X2['flat2']=zs(Tw.flat**2).values; X2['flon2']=zs(Tw.flon**2).values
g2=lg(X2); print('  +geo',g2, g2/lg(X))
X3=X2.copy(); X3['minp']=zs(Tw.minp).values
nd=fx.groupby('track').ndr24.max().reindex(Tw.track).values
X3['ndr']=np.nan_to_num(nd,nan=0.0); X3['ndrmiss']=np.isnan(nd).astype(float)
g3=lg(X3); print('  +minp,ndr',g3,g3/lg(X))
# 8
idx=fx.loc[fx.groupby('track').msl.idxmin()].set_index('track')
la=idx.lat.reindex(Tw.track).values; lo=idx.lon.reindex(Tw.track).values; lo=np.where(lo>180,lo-360,lo)
ind=((la>=55)&(la<=67)&(lo>=-50)&(lo<=-15)).astype(float)
X4=X.copy(); X4['box']=ind; print('8',lg(X4),'n in box',ind.sum())
# 9
rng=np.random.default_rng(1); seasons=np.arange(2004,2026)
sidx={s:np.where(C.season.values==s)[0] for s in seasons}
Xd=np.asarray(S.design(C,z,['GH_lag','NAO_lag']),float)
bs=[]
for _ in range(500):
    pick=rng.choice(seasons,22,replace=True); rows=np.concatenate([sidx[s] for s in pick])
    Cb=C.iloc[rows].reset_index(drop=True); Xb=Xd[rows]
    bs.append(S.fast_pois(Cb.hf.values.astype(float),Xb)[0]-S.fast_pois(Cb['all'].values.astype(float),Xb)[0])
print('9',np.percentile(bs,[2.5,97.5]))
print('--- extra')
b=share(['GH_lag','NAO_lag']); print('HF|NAO exp',np.exp(b['hf'][0]),np.exp(b['all'][0]))
Tw['gt']=pd.qcut(Tw.GH,3,labels=[1,2,3]); Tw['nt']=pd.qcut(Tw.NAO,3,labels=[1,2,3])
print('5b (9622)',Tw.groupby('gt').hf.agg(['mean','size']).round(3).T.to_string())
print(Tw.groupby(['nt','gt']).hf.agg(['mean','size']).round(3).unstack().to_string())
for name,(a,b_,c,d) in {'lon-50..-15 lat55-67':(55,67,-50,-15),'wide':(50,70,-60,-10),'strict':(55,67,-50,-15)}.items():
    i2=((la>=a)&(la<=b_)&(lo>=c)&(lo<=d)).astype(float); X5=X.copy(); X5['box']=i2; print(name,lg(X5),i2.sum())
# alt: box on 0..360 lon? and strict inequality; or peak-gust fix
i3=((la>55)&(la<67)&(lo>-50)&(lo<-15)).astype(float); X5=X.copy(); X5['box']=i3; print('strict ineq',lg(X5))
# box using track peak_lat/peak_lon
pl=Tw.peak_lon.values; pl=np.where(pl>180,pl-360,pl); i4=((Tw.peak_lat.values>=55)&(Tw.peak_lat.values<=67)&(pl>=-50)&(pl<=-15)).astype(float)
X5=X.copy(); X5['box']=i4; print('peak pos',lg(X5),i4.sum())
