import sys, numpy as np, pandas as pd, datetime as dt
sys.path.insert(0,'/home/claude/hf-low/research/era5/freq_split')
from split import read_cpc
R='/home/claude/hf-low/research/era5/'
T=pd.read_csv(R+'hf_history/results/all_tracks.csv.gz',dtype={'start':str})
T=T[T.basin=='atl'].copy()
T['gen']=pd.to_datetime(T.start.str[:8],format='%Y%m%d')
yr=np.where((T.gen.dt.month>=10),T.gen.dt.year,T.gen.dt.year-1)
T['season']=yr; T['day']=(T.gen.values-pd.to_datetime(dict(year=yr,month=10,day=1)).values).astype('timedelta64[D]').astype(int)
T=T[(T.season>=2004)&(T.season<=2025)&(T.day>=0)&(T.day<210)].copy()
print('dup track ids',T.track.duplicated().sum())
T['hf']=T.gust800_kt>=71.7
F=pd.read_csv(R+'highlat/gustloc_fixes.csv')
F['surv']=~((F.max_dkm>400)&(F.max_dgl<300))
g=F.groupby('track')
fl={'T':g.surv.any(),'G100':g.g800_gl100.max()>=71.7,'G300':g.g800_gl300.max()>=71.7,'R400':g.g800_r400.max()>=71.7}
for k,v in fl.items():
    T[k]=T.hf & np.where(T.track.isin(v.index),T.track.map(v).fillna(False).astype(bool),True)
T['barr']=T.hf&~T['T']
print(len(T),T.hf.sum(),{k:int(T[k].sum()) for k in fl},T.barr.sum())
print('hf tracks without fix rows',(T.hf&~T.track.isin(F.track)).sum())
n60=T.peak_lat<60
print('N60',n60.sum(),(T.hf&n60).sum()); dd=~T.barr
print('D',dd.sum(),(T.hf&dd).sum())
nao=read_cpc('/mnt/project-files/teleconnection-test/cpc_indices/norm.daily.nao.index.b500101.current.ascii')
gh=pd.read_csv(R+'nao_share_barrier/gh_daily.csv',parse_dates=['date']).set_index('date').GH.asfreq('D')
lag=lambda s:s.rolling(7,min_periods=7).mean().shift(4)
rows=[(s,k,pd.Timestamp(dt.date(s,10,1)+dt.timedelta(days=k))) for s in range(2004,2026) for k in range(210)]
D=pd.DataFrame(rows,columns=['season','day','date']);D['month']=D.date.dt.month
D['nao']=lag(nao).reindex(D.date).values; D['gh']=lag(gh).reindex(D.date).values
D=D[D.nao.notna()].reset_index(drop=True); print('days',len(D),'gh nan',D.gh.isna().sum())
D['zn']=(D.nao-D.nao.mean())/D.nao.std(ddof=0)
D['zg']=(D.gh-D.gh.mean())/D.gh.std(ddof=0)
print('corr',D.zn.corr(D.zg), 'gh-days',D.gh.notna().sum())
def cnt(mask):
    x=T[mask]; c=x.groupby(['season','day']).size()
    return c.reindex(pd.MultiIndex.from_frame(D[['season','day']]),fill_value=0).values.astype(float)
def X(cols,idx=None):
    d=D if idx is None else D.iloc[idx]
    M=[d[c].values for c in cols]+[(d.month.values==m).astype(float) for m in (11,12,1,2,3,4)]
    M+=[(d.season.values-D.season.mean())/10, np.ones(len(d))]
    return np.column_stack(M)
def fit(y,Xm):
    b=np.zeros(Xm.shape[1]);b[-1]=np.log(y.mean())
    for _ in range(50):
        mu=np.exp(Xm@b);st=np.linalg.solve(Xm.T@(Xm*mu[:,None]),Xm.T@(y-mu));b+=st
        if abs(st).max()<1e-11:break
    return b
Y={'all':cnt(np.ones(len(T),bool)),'hf':cnt(T.hf.values),'T':cnt(T['T'].values),'barr':cnt(T.barr.values)}
for k in ('G100','G300','R400'):Y[k]=cnt(T[k].values)
Y['all_n60']=cnt(n60.values);Y['hf_n60']=cnt((T.hf&n60).values)
Y['all_d']=cnt(dd.values);Y['hf_d']=cnt((T.hf&dd).values)
b=lambda y,cols,k=0,idx=None:fit(y if idx is None else y[idx],X(cols,idx))[k]
e=np.exp
a=b(Y['all'],['zn']);h=b(Y['hf'],['zn']);print('2',e(h),e(a),e(h-a))
for k in ('T','G100','G300','R400'):print('3',k,e(b(Y[k],['zn'])-a))
print('3 N60',e(b(Y['hf_n60'],['zn'])-b(Y['all_n60'],['zn'])),'D',e(b(Y['hf_d'],['zn'])-b(Y['all_d'],['zn'])))
print('4',e(b(Y['T'],['zn'])),e(b(Y['barr'],['zn'])))
# GH models: restrict to rows with GH
ok=D.gh.notna().values
print('6 NAO|GH',e(b(Y['hf'],['zn','zg'],0)-b(Y['all'],['zn','zg'],0)))
print('6 GH only',e(b(Y['hf'],['zg'])-b(Y['all'],['zg'])),e(b(Y['T'],['zg'])-b(Y['all'],['zg'])))
# bootstrap
rng=np.random.default_rng(1);seas=np.arange(2004,2026);byS={s:np.where(D.season.values==s)[0] for s in seas}
out=[]
for _ in range(1500):
    idx=np.concatenate([byS[s] for s in rng.choice(seas,22)])
    A=b(Y['all'],['zn'],0,idx);H=b(Y['hf'],['zn'],0,idx);Tt=b(Y['T'],['zn'],0,idx)
    Ag=b(Y['all'],['zn','zg'],0,idx);Hg=b(Y['hf'],['zn','zg'],0,idx)
    out.append((H-A,Tt-A,Hg-Ag))
o=np.array(out);print('boot SE T-orig',(o[:,1]-o[:,0]).std(ddof=1),'NAO|GH change',(o[:,2]-o[:,0]).std(ddof=1),'CI',e(np.percentile(o[:,0],[2.5,97.5])))
print({k:int(v.sum()) for k,v in Y.items()})
