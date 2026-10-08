import pandas as pd, numpy as np, json, statsmodels.api as sm
R='/home/claude/hf-low/research/era5/'
at=pd.read_csv(R+'hf_history/results/all_tracks.csv.gz')
fx=pd.read_csv(R+'intensity/results/fixes_2004.csv.gz')
fp=pd.read_csv(R+'share_p99/results/fix_p99.csv.gz')
cuts=json.load(open(R+'share_p99/results/cuts.json'))['cuts_kt']
at=at[at.basin.isin(['atl','pac'])].copy()
d=pd.to_datetime(at.start.astype(str),format='%Y%m%d%H')
s=np.where(d.dt.month>=10,d.dt.year,d.dt.year-1)
at['s']=s; at['day']=(d.dt.normalize()-pd.to_datetime(pd.Series(s).astype(str)+'-10-01').values).dt.days.values
at['d0']=d.dt.normalize().values
allt=at.copy()
pop=at[(at.s>=2004)&(at.s<=2025)&(at.day>=0)&(at.day<=209)].copy()
pop['week']=pop.day//7
print('pop tracks',pop.groupby('basin').size().to_dict(), 'dup track ids',pop.track.duplicated().sum())
fx['time']=fx.time.astype('int64'); fp['time']=fp.time.astype('int64')
fx=fx.merge(fp[['track','time','p98','p99']],on=['track','time'],how='left')
fx[['p98','p99']]=fx[['p98','p99']].fillna(0)
print('fix_p99 rows unmatched in fixes_2004:',len(fp)-fx.p99.gt(0).sum())
D=set(fx[fx.dp12<=-3.6].track)
idxs={'M':fx.groupby('track').g800.max(),'P99':fx.groupby('track').p99.max(),'P98':fx.groupby('track').p98.max()}
def hf(df,lab):
    if lab=='ref': return (df.gust800_kt>=71.7).values
    ix=df.track.map(idxs[lab]); return (ix>=cuts[lab]).fillna(False).values
for df in (pop,allt):
    df['D']=df.track.isin(D).values
    for lab in ['ref','M','P99','P98']: df['HF_'+lab]=hf(df,lab)
print('counts per basin/label: HF, D, HD')
for b in ['atl','pac']:
    p=pop[pop.basin==b]
    for lab in ['ref','M','P99','P98']:
        print(b,lab,'HF',p['HF_'+lab].sum(),'D',p.D.sum(),'HD',(p.D&p['HF_'+lab]).sum())
print('(3) D tracks with M label != P99 label')
for b in ['atl','pac']:
    p=pop[(pop.basin==b)&pop.D]; print(b,(p.HF_M!=p.HF_P99).sum(),'M-only',(p.HF_M&~p.HF_P99).sum(),'P99-only',(~p.HF_M&p.HF_P99).sum())
def rr(b,lab,restrict):
    src=pop if restrict else allt
    src=src[src.basin==b]
    idx=pd.read_csv(R+f'hemispheric/results/oos_index_{"atl" if b=="atl" else "pac"}.csv').sort_values(['season','week'])
    z=((idx.idx-idx.idx.mean())/idx.idx.std(ddof=1)).values
    ss=idx.season.values; wk=idx.week.values
    coefs={}
    for out in ['D','HD']:
        y=src.D.values if out=='D' else (src.D&src['HF_'+lab]).values
        dd=src.d0.values
        yy=[];lag=[];X=[]
        for sea,w in zip(ss,wk):
            st=pd.Timestamp(f'{sea}-10-01')+pd.Timedelta(days=7*int(w))
            e=st+pd.Timedelta(days=7)
            m=(dd>=np.datetime64(st))&(dd<np.datetime64(e))
            l=(dd>=np.datetime64(st-pd.Timedelta(days=7)))&(dd<np.datetime64(st))
            yy.append(y[m].sum()); lag.append(y[l].sum())
            mo=(st+pd.Timedelta(days=3)).month
            X.append([1.0*(mo==k) for k in [10,11,12,1,2,3,4]])
        X=np.column_stack([np.array(X),np.log1p(lag),z])
        r=sm.GLM(np.array(yy),X,family=sm.families.Poisson()).fit()
        coefs[out]=r.params[-1]
    return np.exp(coefs['HD']-coefs['D'])
claim={('ref','atl'):1.252,('ref','pac'):1.137,('M','atl'):1.203,('M','pac'):1.137,('P99','atl'):1.180,('P99','pac'):1.140,('P98','atl'):1.180,('P98','pac'):1.151}
print('(2) RR: label basin claimed recomputed(pop-only lag) recomputed(lag incl pre-Oct1 tracks)')
for (lab,b),c in claim.items():
    a=rr(b,lab,True); a2=rr(b,lab,False)
    print(lab,b,c,round(a,4),round(a2,4),'MATCH' if abs(a-c)<=.002 else 'MISMATCH')
