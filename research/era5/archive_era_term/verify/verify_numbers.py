import json, numpy as np, pandas as pd, statsmodels.api as sm
from datetime import datetime, timedelta
R='/home/claude/hf-low/'
ERA=datetime(2009,11,23)
out=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s); out.append(s)
def ts(v): return datetime.strptime(str(int(v)),'%Y%m%d%H')
d=json.load(open(R+'docs/data/hf-lows.json'))
F=d['lowFields']; lows=[dict(zip(F,r)) for r in d['lows']]
# ---------- A
P('== SET A')
for b in ['atl','pac']:
    idx=pd.read_csv(R+f'research/era5/hemispheric/results/oos_index_{b}.csv')
    z=((idx.idx-idx.idx.mean())/idx.idx.std(ddof=0)).values
    ws=[datetime(s,10,1)+timedelta(days=7*k) for s in range(2004,2026) for k in range(30)]
    assert len(ws)==len(idx)==660
    bl=[l for l in lows if l['basin']==b and l['cls']=='low']
    st=np.array([ts(l['start']) for l in bl],dtype='datetime64[h]')
    hf=np.array([l['hfN'] for l in bl])
    mon=np.array([(w+timedelta(days=3)).month for w in ws])
    era=np.array([1.0 if w>=ERA else 0.0 for w in ws])
    P(b,'corr(idx,era)=%.3f'%np.corrcoef(z,era)[0,1])
    for name,mask in [('h1',hf>=0),('h2',hf>=2),('h3',hf>=3),('b',hf==1),('m',hf==2)]:
        s=st[mask]
        cnt=np.array([((s>=np.datetime64(w,'h'))&(s<np.datetime64(w+timedelta(days=7),'h'))).sum() for w in ws])
        prev=np.array([((s>=np.datetime64(w-timedelta(days=7),'h'))&(s<np.datetime64(w,'h'))).sum() for w in ws])
        D=np.column_stack([(mon==m).astype(float) for m in [10,11,12,1,2,3,4]]+[np.log1p(prev),z])
        r1=sm.GLM(cnt,D,family=sm.families.Poisson()).fit()
        D2=np.column_stack([D,era]); r2=sm.GLM(cnt,D2,family=sm.families.Poisson()).fit()
        P(b,name,'total',cnt.sum(),'M1 RR %.3f'%np.exp(r1.params[-1]),'M2 RR %.3f'%np.exp(r2.params[-2]),'era coef %+.3f'%r2.params[-1])
# ---------- B
P('== SET B')
EXB=(ERA-datetime(2004,6,1)).days/365.25; EXA=(datetime(2026,6,1)-ERA).days/365.25
P('exposure',round(EXB,2),round(EXA,2))
for b in ['atl','pac']:
    ev=[]
    for l in lows:
        if l['basin']!=b or not(2004<=l['season']<=2025): continue
        seen={};
        for f in l['fixes']:
            seen.setdefault(f[0],f)
        fx=sorted(seen.values(),key=lambda f:f[0])
        h=[f for f in fx if f[3]=='HF']
        if h: ev.append((ts(h[0][0]),len(h)))
    for k in (1,2,3):
        e=[t for t,n in ev if n>=k]
        eb=sum(t<ERA for t in e); ea=len(e)-eb
        P(b,'>=%d'%k,'per season %.1f'%(len(e)/22),'before %.1f after %.1f'%(eb/EXB,ea/EXA))
    for k in (2,3):
        sh=[]
        for lo in (True,False):
            tot=[n for t,n in ev if (t<ERA)==lo]
            sh.append(sum(n>=k for n in tot)/len(tot))
        P(b,'share >=%d before %.3f after %.3f'%(k,*sh))
# ---------- C
P('== SET C')
for src in ['archive','proxy']:
    e=pd.read_csv(R+f'research/era5/climatology_atlas/results/events_{src}.csv')
    e['t']=pd.to_datetime(e.t_on)
    for b in ['atl','pac']:
        x=e[e.basin==b]; bf=x[x.t<ERA]; af=x[x.t>=ERA]
        P(src,b,'mean %.2f -> %.2f'%(bf.hf_h.mean(),af.hf_h.mean()),'n',len(bf),len(af),'share6 %.3f -> %.3f'%((bf.hf_h==6).mean(),(af.hf_h==6).mean()))
# ---------- D
P('== SET D')
t=pd.read_csv(R+'research/era5/hf_history/results/all_tracks.csv.gz')
m=pd.read_csv(R+'research/era5/gust_depth/results/matches.csv')
t=t[(t.season>=2004)&(t.season<=2025)&(t.gust800_kt>=71.7)].copy()
P('n tracks',len(t))
t['listed']=t.track.isin(set(m.track.dropna().astype(int)))
t['pt']=t.peak_time.map(ts); t['bef']=t.pt<ERA
for b in ['atl','pac']:
    x=t[t.basin==b]
    P(b,'before %.3f (n %d) after %.3f (n %d)'%(x[x.bef].listed.mean(),x.bef.sum(),x[~x.bef].listed.mean(),(~x.bef).sum()))
P('overall %.3f'%t.listed.mean())
# ---------- E
P('== SET E')
fp=pd.read_csv(R+'research/era5/hf_probability_tracks/results/fix_probs_2004.csv.gz')
P('rows',len(fp))
fp['t']=pd.to_datetime(fp.time.astype(str),format='%Y%m%d%H')
frames=[]
for b in ['Atl','Pac']:
    a=pd.read_csv(R+f'data/hf_lows/HF_Data_-_{b}.csv'); a=a[a.Category=='HF'].copy()
    a['t']=pd.to_datetime(a.date.astype(str),format='%Y%m%d%H',errors='coerce')
    P(b,'malformed date rows dropped:',int(a.t.isna().sum())); a=a.dropna(subset=['t']); frames.append(a)
a=pd.concat(frames,ignore_index=True)
a['season']=a.t.dt.year-(a.t.dt.month<6)
a=a[(a.season>=2004)&(a.season<=2025)].reset_index(drop=True)
P('archive HF fixes',len(a))
def hav(la1,lo1,la2,lo2):
    la1,lo1,la2,lo2=map(np.radians,(la1,lo1,la2,lo2))
    h=np.sin((la2-la1)/2)**2+np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*6371*np.arcsin(np.sqrt(h))
g={k:v for k,v in fp.groupby('t')}
matched=[]  # (track, archive time)
unm=0
for r in a.itertuples():
    lon=r.Longitude%360; hr=r.t.hour
    if hr in (0,12): cands=[r.t]; lim=600
    elif hr in (6,18): cands=[r.t-timedelta(hours=6),r.t+timedelta(hours=6)]; lim=800
    else: unm+=1; continue
    best=None
    for c in cands:
        if c not in g: continue
        x=g[c]; dd=hav(r.Latitude,lon,x.lat.values,x.lon.values%360)
        i=dd.argmin()
        if dd[i]<=lim and (best is None or dd[i]<best[0]): best=(dd[i],x.track.values[i])
    if best: matched.append((best[1],r.t))
    else: unm+=1
P('matched',len(matched),'unmatched',unm)
mt={}
for tr,tt in matched: mt.setdefault(tr,[]).append(tt)
lab=np.zeros(len(fp),int)
trk=fp.track.values; tm=fp.t.values
for i in range(len(fp)):
    L=mt.get(trk[i])
    if L:
        t0=fp.t.iloc[i]
        lab[i]=int(any(t0<x<=t0+timedelta(hours=24) for x in L))
fp['arch24']=lab
def hss(f,o):
    f=np.asarray(f).astype(bool);o=np.asarray(o).astype(bool)
    a_=(f&o).sum();b_=(f&~o).sum();c_=(~f&o).sum();d_=(~f&~o).sum()
    return 2*(a_*d_-b_*c_)/((a_+c_)*(c_+d_)+(a_+b_)*(b_+d_))
fp['fc']=fp.P_F>=fp.cut_F
bef=fp.t<ERA
P('positives hf24',int(fp.hf24.sum()),'arch24',int(fp.arch24.sum()))
P('n before',int(bef.sum()),'arch24 pos before',int(fp.arch24[bef].sum()),'n after',int((~bef).sum()),'pos after',int(fp.arch24[~bef].sum()))
for nm,col in [('hf24','hf24'),('fc','fc')]:
    P(nm,'HSS overall %.3f before %.3f after %.3f'%(hss(fp[col],fp.arch24),hss(fp[col][bef],fp.arch24[bef]),hss(fp[col][~bef],fp.arch24[~bef])))
open(R+'research/era5/archive_era_term/verify/VERIFICATION_raw.txt','w').write('\n'.join(out)+'\n')
