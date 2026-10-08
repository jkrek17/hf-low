import numpy as np, pandas as pd, statsmodels.api as sm, warnings
warnings.filterwarnings('ignore')
R='/home/claude/hf-low/research/era5/'
def dt(x): return pd.to_datetime(x.astype('int64').astype(str),format='%Y%m%d%H')
def basin_of(lat,lon):
    a=(lat>=30)&(lat<=67)&((lon>=262)|(lon<=10))
    p=(lat>=27)&(lat<=67)&(lon>=135)&(lon<=240)
    return np.where(a,'atl',np.where(p,'pac',''))
def load_tracker(name):
    f=pd.read_csv(f'/mnt/project-files/second-tracker/{name}_fixes.csv.gz')
    f['time']=f.time.astype('int64'); print(name,'tracks on read',f.track.nunique())
    f['dt']=dt(f.time)
    return f
def tracker_table(name,thr=-3.6):
    f=load_tracker(name)
    if name=='V':
        m=pd.read_csv('/home/claude/hf-low/research/era5/second_tracker_conversion/work/v_mslp500.csv.gz')
        f=f.merge(m,on=['track','time'],how='left'); f['p']=f.mslp500
        idx=f.groupby('track').val.idxmax()
    else:
        f['p']=f.val; idx=f.groupby('track').val.idxmin()
    pk=f.loc[idx.values,['track','lat','lon']]; pk['basin']=basin_of(pk.lat.values,pk.lon.values)
    first=f.groupby('track').dt.min().rename('start')
    # deepening
    g=f[f.dt.dt.hour.isin([0,12])&f.p.notna()][['track','dt','p']].copy()
    g2=g.copy(); g2['dt']=g2.dt+pd.Timedelta(hours=12); g2=g2.rename(columns={'p':'pprev'})
    j=g.merge(g2,on=['track','dt'])
    j['dp']=j.p-j.pprev
    D=set(j.track[j.dp<=thr])
    ev=pd.read_csv(R+f'second_tracker/results/events_{name}.csv')
    hf=set(ev.tr_track[ev.matched==True].astype(int))
    t=pk.merge(first,on='track'); t['D']=t.track.isin(D); t['HF']=t.track.isin(hf)
    return t[t.basin!=''][['track','basin','start','D','HF']]
def A_table():
    a=pd.read_csv(R+'hf_history/results/all_tracks.csv.gz')
    fx=pd.read_csv(R+'intensity/results/fixes_2004.csv.gz',usecols=['track','dp12'])
    D=set(fx.track[fx.dp12<=-3.6])
    a['start']=dt(a.start); a['D']=a.track.isin(D); a['HF']=a.gust800_kt>=71.7
    return a[a.basin.isin(['atl','pac'])][['track','basin','start','D','HF']]
def weekly(t,basin):
    t=t[t.basin==basin]
    rows=[]
    for s in range(2004,2022):
        o=pd.Timestamp(s,10,1)
        for w in range(30):
            ws=o+pd.Timedelta(days=7*w); we=ws+pd.Timedelta(days=7)
            cur=t[(t.start>=ws)&(t.start<we)]; prv=t[(t.start>=ws-pd.Timedelta(days=7))&(t.start<ws)]
            r=dict(season=s,week=w,month=ws.month)
            for k,fn in [('N',lambda x:len(x)),('D',lambda x:x.D.sum()),('HD',lambda x:(x.D&x.HF).sum())]:
                r[k]=fn(cur); r['p'+k]=fn(prv)
            rows.append(r)
    return pd.DataFrame(rows)
def fit(df,out):
    X=pd.DataFrame({f'm{m}':(df.month==m).astype(float) for m in [10,11,12,1,2,3,4]})
    X['lag']=np.log1p(df['p'+out]); X['idx']=df.z
    r=sm.GLM(df[out],X,family=sm.families.Poisson()).fit()
    return r.params['idx']
res=[]
tabs={'A':A_table(),'V':tracker_table('V'),'M':tracker_table('M')}
for tn,t in tabs.items():
    for b in ['atl','pac']:
        w=weekly(t,b)
        ix=pd.read_csv(R+f'hemispheric/results/oos_index_{b}.csv').sort_values(['season','week'])
        ix=ix[ix.season.between(2004,2021)]; assert len(ix)==540
        w=w.merge(ix[['season','week','idx']],on=['season','week'])
        w['z']=(w.idx-w.idx.mean())/w.idx.std(ddof=0)
        bN,bD,bHD=[fit(w,o) for o in ['N','D','HD']]
        res.append((tn,b,w.N.sum(),w.D.sum(),w.HD.sum(),np.exp(bHD-bD),np.exp(bD-bN)))
r=pd.DataFrame(res,columns=['tracker','basin','N','D','HD','RR_A2','RR_A1'])
print(r.round(4).to_string(index=False))
# spot check
import random
random.seed(1)
m=pd.read_csv('/home/claude/hf-low/research/era5/second_tracker_conversion/work/v_mslp500.csv.gz')
V=load_tracker('V').set_index(['track','time'])
lats=19.5+1.5*np.arange(38); lons=1.5*np.arange(240)
ok=0;n=0
for _,row in m.sample(30,random_state=7).iterrows():
    k=(int(row.track),int(row.time)); lat,lon=V.loc[k,'lat'],V.loc[k,'lon']
    d=pd.Timestamp(str(int(row.time)),tz=None) if False else pd.to_datetime(str(int(row.time)),format='%Y%m%d%H')
    h=int((d-pd.Timestamp(1959,1,1)).total_seconds()//3600)
    arr=np.load(f'/home/claude/hf-low/research/era5/second_tracker_conversion/work/mslp/{h//48}.npy')[(h%48)//6]
    LA,LO=np.meshgrid(lats,lons,indexing='ij')
    dl=np.radians(LO-lon); p1=np.radians(LA); p2=np.radians(lat)
    c=np.sin(p1)*np.sin(p2)+np.cos(p1)*np.cos(p2)*np.cos(dl)
    dist=6371*np.arccos(np.clip(c,-1,1))
    v=arr[dist<=500].min(); n+=1; ok+=abs(v-row.mslp500)<=0.01
    print(k,round(v,3),round(row.mslp500,3))
print('spot-check matched',ok,'of',n)
