import pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
import os
R=os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','..','..'))+'/'  # repo root
def hav(la1,lo1,la2,lo2):
    la1,lo1,la2,lo2=map(np.radians,(la1,lo1,la2,lo2))
    a=np.sin((la2-la1)/2)**2+np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*6371*np.arcsin(np.sqrt(a))
def dt(s): return pd.to_datetime(s.astype(str),format='%Y%m%d%H')
arc={}
for b,f in (('atl','Atl'),('pac','Pac')):
    d=pd.read_csv(R+f'data/hf_lows/HF_Data_-_{f}.csv',dtype={'date':str})
    d=d[d.Category=='HF']
    d=d[d.date.astype(str).str.fullmatch(r'\d{10}')].copy()
    d['t']=pd.to_datetime(d.date,format='%Y%m%d%H',errors='coerce'); d=d.dropna(subset=['t'])
    d['lon']=d.Longitude%360; d['lat']=d.Latitude
    arc[b]=d.sort_values('t').reset_index(drop=True)
g=pd.concat([pd.read_csv(R+f'research/era5/gfs_test/results/gfs_tracks_{y}.csv.gz') for y in range(2021,2026)])
g['t']=dt(g.valid_time)
g=g[g.basin.isin(['atl','pac'])&g.t.dt.hour.isin([0,12])&g.lead.isin([0,12,24,36,48])].copy()
g['season']=np.where(g.t.dt.month>=6,g.t.dt.year,g.t.dt.year-1)
g=g[g.season.between(2021,2025)].reset_index(drop=True)
g['lon']=g.lon%360
def label(df,lo,hi,tcol='t'):
    y=np.zeros(len(df),int)
    for b in ('atl','pac'):
        a=arc[b]; at=a.t.values.astype('datetime64[s]').astype(np.int64)
        idx=np.where(df.basin.values==b)[0]
        tt=df[tcol].values.astype('datetime64[s]').astype(np.int64)
        for i in idx:
            l=np.searchsorted(at,tt[i]+lo,'right') if lo_open else np.searchsorted(at,tt[i]+lo,'left')
            r=np.searchsorted(at,tt[i]+hi,'right')
            if r>l:
                s=a.iloc[l:r]
                if (hav(df.lat.values[i],df.lon.values[i],s.lat.values,s.lon.values)<=600).any(): y[i]=1
    return y
H=3600
lo_open=False; g['Y0']=label(g,-3*H,3*H)
lo_open=True; g['Y24']=label(g,0,24*H)
print('total',len(g),g.Y0.sum(),g.Y24.sum())
def sc(f,y):
    a=int((f&y).sum());b=int((f&~y).sum());c=int((~f&y).sum());d=int((~f&~y).sum());n=a+b+c+d
    pod=a/(a+c);far=b/(a+b) if a+b else np.nan;bias=(a+b)/(a+c)
    den=(a+c)*(c+d)+(a+b)*(b+d);hss=2*(a*d-b*c)/den
    return round(hss,3),round(pod,3),round(far,3),round(bias,3)
for b in ('atl','pac'):
    s=g[(g.basin==b)&(g.lead==0)]; print(b,'lead0 n',len(s),s.Y0.sum())
    print(' nominal L0',sc((s.g800>64).values,s.Y0.values==1))
    s24=g[(g.basin==b)&(g.lead==24)]
    print(' nominal L24 Y0',sc((s24.g800>64).values,s24.Y0.values==1))
    print(' AUC L0',round(roc_auc_score(s.Y0,s.g800),3))
    T={}
    for ss in range(2021,2026):
        o=s[s.season!=ss]; k=int(o.Y0.sum()); T[ss]=np.sort(o.g800.values)[::-1][k]
    print(' T',[round(T[x],1) for x in T])
    for L in (0,24):
        x=g[(g.basin==b)&(g.lead==L)]; thr=x.season.map(T)
        print(f' matched L{L}',sc((x.g800>thr).values,x.Y0.values==1))
    globals()['T_'+b]=T
# ERA5
e=pd.read_csv(R+'research/era5/intensity/results/fixes_2004.csv.gz')
e['t']=dt(e.time); e=e[e.basin.isin(['atl','pac'])&e.t.dt.hour.isin([0,12])&e.t.dt.month.isin([10,11,12,1,2,3,4])&e.season.between(2021,2025)].reset_index(drop=True)
e['lon']=e.lon%360
lo_open=False; e['Y0']=label(e,-3*H,3*H)
print('ERA5',len(e),e.Y0.sum())
for b in ('atl','pac'):
    s=e[e.basin==b]; print(b,sc((s.g800>71.7).values,s.Y0.values==1))
# events
rows=[]
for b in ('atl','pac'):
    a=arc[b]
    f=a.sort_values('t').groupby('ID').head(1)
    f=f[f.t.dt.month.isin([10,11,12,1,2,3,4])].copy()
    f['season']=np.where(f.t.dt.month>=6,f.t.dt.year,f.t.dt.year-1)
    f=f[f.season.between(2021,2025)]; f['b']=b; rows.append(f)
ev=pd.concat(rows); print('events',len(ev),ev.groupby('b').size().to_dict())
def vt(t):
    fl=t.floor('12h'); return fl if t.hour%12<6 else fl+pd.Timedelta(hours=12)
n_cov=n0=nany=0; n0b=nanyb=0
for _,r in ev.iterrows():
    v=vt(r.t); c=g[(g.t==v)]
    c=c[hav(r.lat,r.lon,c.lat.values,c.lon.values)<=600]
    if len(c)==0: continue
    n_cov+=1
    T=(T_atl if r.b=='atl' else T_pac)[r.season]
    if ((c.lead==0)&(c.g800>T)).any(): n0+=1
    if (c.g800>T).any(): nany+=1
    cb=c[c.basin==r.b]
    if ((cb.lead==0)&(cb.g800>T)).any(): n0b+=1
    if (cb.g800>T).any(): nanyb+=1
print('cov',n_cov,'lead0',n0,'any',nany,'| own-basin fixes only',n0b,nanyb)
