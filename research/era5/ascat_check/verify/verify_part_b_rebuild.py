import sys,glob,numpy as np,pandas as pd,datetime as dt
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from blosc_nolib import decode
d=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','results','scat_footprints.csv.gz'));d=d[d.cover>=0.40].reset_index(drop=True)
rng=np.random.default_rng(1);p=rng.permutation(len(d));s=d.iloc[p[:40]]
if (s.idx>=71.7).sum()<10:
    hi=d.iloc[[i for i in p if d.idx[i]>=71.7]]; lo=d.iloc[[i for i in p if d.idx[i]<71.7]]
    s=pd.concat([hi.iloc[:10],lo.iloc[:30]])
print('hi in sample',(s.idx>=71.7).sum())
T=pd.read_csv('research/era5/hf_history/results/all_tracks.csv.gz').set_index('track')
lsm=np.load('/mnt/project-files/era5-cache/lsm_nh_0p25.npy')
H={}
for f in sorted(glob.glob('/mnt/project-files/era5-cache/gust_nh_12h/gust_*.npz')):
    z=np.load(f)
    for k,h in enumerate(z['hours_since_1900']):H[int(h)]=(f,k)
cache={};cells={}
def field(h):
    f,k=H[h]
    if f not in cache:
        cache.clear();z=np.load(f);cache[f]=(decode(z['data'].tobytes()).reshape(tuple(z['shape'])),)
    return cache[f][0][k]/2.0
def gc(la1,lo1,la2,lo2):
    la1,lo1,la2,lo2=map(np.radians,(la1,lo1,la2,lo2))
    a=np.sin((la2-la1)/2)**2+np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*6371*np.arcsin(np.sqrt(a))
def ocean(lat,lon):
    r=np.floor((90-lat)/.25).astype(int);c=np.floor(lon/.25).astype(int)
    return (lsm[r,c]+lsm[r+1,c]+lsm[r,(c+1)%1440]+lsm[r+1,(c+1)%1440])/4<0.5
def bil(F,lat,lon):
    y=(90-lat)/.25;x=lon/.25;r=np.floor(y).astype(int);c=np.floor(x).astype(int);fy=y-r;fx=x-c
    c1=(c+1)%1440;c=c%1440
    return F[r,c]*(1-fy)*(1-fx)+F[r+1,c]*fy*(1-fx)+F[r,c1]*(1-fy)*fx+F[r+1,c1]*fy*fx
rows=[];ok=0;idxok=0
for _,r in s.iterrows():
    t=T.loc[int(r.track)];pt=str(int(t.peak_time));assert pt==str(int(r.peak_time)),(pt,r.peak_time)
    D=dt.datetime.strptime(pt,'%Y%m%d%H');assert D.hour in(0,12)
    mins=int((D-dt.datetime(1970,1,1)).total_seconds()//60)
    yrs=[D.year-1,D.year,D.year+1] ; arr=[]
    for y in(D.year,D.year-1,D.year+1):
        fn=f'/mnt/project-files/ascat-check/cells/cells_{r.sensor}_{y}.npz'
        if (fn) not in cells:
            try:cells[fn]=np.load(fn)['rows']
            except Exception:cells[fn]=np.zeros((0,5),np.int64)
        a=cells[fn];arr.append(a[np.abs(a[:,2]-mins)<=90])
    a=np.unique(np.concatenate(arr),axis=0)
    a=a[(a[:,4]&1)==0]
    lat=-89.875+.25*a[:,0];lon=.125+.25*a[:,1]
    m=gc(lat,lon,t.peak_lat,t.peak_lon)<=800
    a,lat,lon=a[m],lat[m],lon[m]
    o=ocean(lat,lon);a,lat,lon=a[o],lat[o],lon[o]
    h=int((D-dt.datetime(1900,1,1)).total_seconds()//3600);F=field(h)
    n=len(a);sm=(a[:,3]*0.2*1.94384449).max() if n else np.nan
    em=bil(F,lat,lon).max() if n else np.nan
    # idx check
    R=np.arange(361)[:,None];C=np.arange(1440)[None,:]
    la=90-.25*R;lo=.25*C
    dd=gc(la,lo,t.peak_lat,t.peak_lon)
    ocn=lsm<0.5
    ix=F[(dd<=800)&ocn].max()
    good=abs(n-r.n_cells)<=1 and abs(sm-r.scat_max)<=.1 and abs(em-r.era5_max_fp)<=.1
    ok+=good;idxok+=abs(ix-r.idx)<=1
    rows.append((int(r.track),r.sensor,pt,r.idx,r.n_cells,n,round(r.scat_max,2),round(sm,2),round(r.era5_max_fp,2),round(em,2),round(ix,2),good))
print(pd.DataFrame(rows,columns='track sens pt idx n_file n_mine sc_file sc_mine e_file e_mine idx_mine ok'.split()).to_string())
print('match',ok,'/40; idx within1',idxok/40)
