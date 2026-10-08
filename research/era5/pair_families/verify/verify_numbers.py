import numpy as np, pandas as pd
R='/home/claude/hf-low/research/era5/'
tr=pd.read_csv(R+'hf_history/results/all_tracks.csv.gz')
fx=pd.read_csv(R+'intensity/results/fixes_2004.csv.gz')
lc=pd.read_csv(R+'hf_history/results/lifecycle_events.csv')
def hrs(x):
    d=pd.to_datetime(x.astype(str),format='%Y%m%d%H'); return ((d-pd.Timestamp('2000-01-01')).dt.total_seconds()/3600).values
def bear(la1,lo1,la2,lo2):
    la1,lo1,la2,lo2=map(np.radians,(la1,lo1,la2,lo2)); dl=lo2-lo1
    y=np.sin(dl)*np.cos(la2); x=np.cos(la1)*np.sin(la2)-np.sin(la1)*np.cos(la2)*np.cos(dl)
    return np.degrees(np.arctan2(y,x))%360
def dist(la1,lo1,la2,lo2):
    la1,lo1,la2,lo2=map(np.radians,(la1,lo1,la2,lo2))
    a=np.sin((la2-la1)/2)**2+np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 6371.0*2*np.arcsin(np.sqrt(a))
def pairs(T,shift_fn=None):
    """T: dict with arrays. returns (n_pairs,n_rear,n_og,n_og_rear)"""
    pass
def run(basin,nperm=300,seed=20261008):
    t=tr[(tr.basin==basin)&(tr.gust800_kt>=71.7)&(tr.season>=2004)].reset_index(drop=True)
    f=fx[fx.track.isin(t.track)&(fx.basin==basin)].copy(); f['h']=hrs(f.time); f=f.sort_values(['track','h'])
    lcb=lc[lc.basin==basin].drop_duplicates('track').set_index('track')
    n=len(t); F=[]
    for k in t.track:
        g=f[f.track==k]; lon=np.unwrap(np.radians(g.lon.values)); F.append((g.h.values,g.lat.values,np.degrees(lon)))
    ok=np.array([len(x[0])>0 for x in F]); 
    peak=hrs(t.peak_time); start=hrs(t.start)
    first=np.array([x[0][0] if len(x[0]) else np.nan for x in F]); last=np.array([x[0][-1] if len(x[0]) else np.nan for x in F])
    fm=f.groupby('track').msl.first().reindex(t.track).values
    go=lcb.gen_obs.reindex(t.track).fillna(False).astype(bool).values
    og=go&((first-start)<=12)
    flat=np.array([x[1][0] if len(x[0]) else np.nan for x in F]); flon=np.array([x[2][0] if len(x[0]) else np.nan for x in F])
    def count(sh):
        npair=nr=ng=ngr=0
        for d in range(n):
            if not ok[d]: continue
            td=first[d]+sh[d]
            for p in range(n):
                if p==d or not ok[p]: continue
                h,la,lo=F[p]; hp=h+sh[p]
                if not (hp[0]<td<=hp[-1]): continue
                if abs(td-(peak[p]+sh[p]))>48: continue
                if len(h)<2: continue
                i=min(max(np.searchsorted(hp,td,side='left')-1,0),len(h)-2)  # bracket (t_i, t_i+1]; D first fix sits on a fix time
                w=(td-hp[i])/(hp[i+1]-hp[i]); pla=la[i]+w*(la[i+1]-la[i]); plo=lo[i]+w*(lo[i+1]-lo[i])
                dlo=((flon[d]-plo+180)%360)-180
                if dist(pla,plo,flat[d],plo+dlo)>1500: continue
                head=bear(la[i],lo[i],la[i+1],lo[i+1]); b=bear(pla,plo,flat[d],plo+dlo)
                rear=120<=((b-head)%360)<=240
                npair+=1; nr+=rear
                if og[d]: ng+=1; ngr+=rear
        return npair,nr,ng,ngr
    obs=count(np.zeros(n)); print(basin,'tracks',n,'pairs,rear,og,ogrear',obs,'f=%.4f'%(obs[1]/obs[0]),flush=True)
    rng=np.random.default_rng(seed)
    mon=((t.peak_time//10000)%100).values; grp=pd.Series(range(n)).groupby([t.season.values,mon]).apply(list)
    fs=[]
    for _ in range(nperm):
        sh=np.zeros(n)
        for idx in grp:
            idx=np.array(idx); newp=rng.permutation(peak[idx]); sh[idx]=newp-peak[idx]
        c=count(sh); fs.append(c[1]/c[0] if c[0] else np.nan)
    fs=np.array(fs); fo=obs[1]/obs[0]
    print(basin,'null f mean %.4f'%np.nanmean(fs),'p_one_sided %.4f'%((1+np.nansum(fs>=fo))/(1+np.sum(~np.isnan(fs)))))
import sys
run(sys.argv[1],int(sys.argv[2]) if len(sys.argv)>2 else 300)
