import numpy as np, pandas as pd
from scipy.optimize import minimize
R='/home/claude/hf-low/research/era5/'
ev=pd.read_csv(R+'hf_history/results/lifecycle_events.csv')
ev=ev[(ev.era=='2004+')&(ev.tc==False)].copy()
cat=pd.read_csv(R+'hf_history/results/era5_hf_catalog.csv')
cat=cat[cat.role=='event'].drop_duplicates('track')[['track','gust800_kt','peak_time']]
tr=pd.read_csv(R+'hf_history/results/era5_hf_catalog_tracks.csv')
on=tr[tr.basin.notna()&(tr.g800>=71.7)].sort_values('time').groupby('track').first()[['time','msl']].rename(columns={'time':'t_on','msl':'msl_on_tr'})
fx=pd.read_csv(R+'hf_structure/results/fixes.csv').sort_values('time')
nfx=fx.groupby('track').size().rename('nfx')
f1=fx.groupby('track').first()
d=ev.merge(cat,on='track',how='left').merge(on,left_on='track',right_index=True,how='left').merge(f1[['time','gmax_r','gmax_rel','gmax_lat','speed_kt','msl','terrain']],left_on='track',right_index=True,how='left').merge(nfx,left_on='track',right_index=True,how='left')
print('n',len(d),d.basin.value_counts().to_dict(),'missing fixes',d.gmax_r.isna().sum(),'time mismatch',(d.time!=d.t_on).sum(),'no t_on',d.t_on.isna().sum())
d['pac']=(d.basin=='pac').astype(int)
d['late']=(d.h_on_minp>=0).astype(int)
m=d.mon_on
d['mg']=np.where(m.isin([12,1,2]),'DJF',np.where(m.isin([3,4,5]),'MAM','SON'))
print('\n== 2 LATE shares')
print('all %.3f atl %.3f pac %.3f'%(d.late.mean(),d[d.pac==0].late.mean(),d[d.pac==1].late.mean()))
print('>=6 %.3f >=12 %.3f'%((d.h_on_minp>=6).mean(),(d.h_on_minp>=12).mean()))
d['T1']=(d.gmax_r>400).astype(float)
d['T2']=((d.gmax_lat>=60)&(d.pac==0)).astype(float)
d['T3']=(d.gust800_kt<76.7).astype(float)
d['T4']=((d.gmax_rel>=90)&(d.gmax_rel<270)).astype(float)
d['T5']=np.where(d.maxdeep.isna(),np.nan,(d.maxdeep>=1).astype(float))
d['T6']=(d.speed_kt-d.speed_kt.mean())/d.speed_kt.std()
d['T7']=d.pac.astype(float)
print('\n== 3 prevalences')
print('T1 n',int(d.T1.sum()),'atl %.3f pac %.3f'%(d[d.pac==0].T1.mean(),d[d.pac==1].T1.mean()))
print('T2 n',int(d.T2.sum()),'share of atl %.3f'%d[d.pac==0].T2.mean())
for t in['T3','T4','T5']: print(t,'atl %.3f pac %.3f'%(d[d.pac==0][t].mean(),d[d.pac==1][t].mean()))
print('maxdeep missing',d.maxdeep.isna().sum(),'; T6 mean %.1f sd %.1f'%(d.speed_kt.mean(),d.speed_kt.std()),'speed missing',d.speed_kt.isna().sum())
print('\n== 4 raw late share with/without')
for t,sub in [('T1',d),('T2',d[d.pac==0]),('T3',d),('T4',d),('T5',d),('T7',d)]:
    s=sub.dropna(subset=[t]); print(t,'with %.3f (n=%d) without %.3f (n=%d)'%(s[s[t]==1].late.mean(),(s[t]==1).sum(),s[s[t]==0].late.mean(),(s[t]==0).sum()))
# standardised (basin-month cells, pooled weights)
def std(t,sub):
    s=sub.dropna(subset=[t]); s=s.assign(cell=s.basin+s.mg); w=s.cell.value_counts(normalize=True); out=[]
    for v in(1,0):
        r=s[s[t]==v].groupby('cell').late.mean(); out.append((r*w.reindex(r.index)).sum()/w.reindex(r.index).sum())
    return out
print('standardised (my own, basin x month cells pooled weights):')
for t,sub in [('T1',d),('T2',d[d.pac==0]),('T3',d),('T4',d),('T5',d),('T7',d)]: print(t,'%.3f vs %.3f'%tuple(std(t,sub)))
# logistic
def fit(X,y,l2=1e-8):
    X=np.asarray(X,float);y=np.asarray(y,float)
    f=lambda b:(np.logaddexp(0,X@b)-y*(X@b)).sum()+l2*(b@b)
    g=lambda b:X.T@(1/(1+np.exp(-X@b))-y)+2*l2*b
    return minimize(f,np.zeros(X.shape[1]),jac=g,method='BFGS').x
def design(s,preds,basin=True):
    cols=[np.ones(len(s))]+[s[p].values for p in preds]
    if basin: cols.append(s.pac.values)
    cols+= [(s.mg=='DJF').values.astype(float),(s.mg=='MAM').values.astype(float)]
    return np.column_stack(cols)
print('\n== 5 odds ratios')
for t,sub,bas in [('T1',d,1),('T2',d[d.pac==0],0),('T3',d,1),('T4',d,1),('T5',d,1),('T6',d,1),('T7',d,0)]:
    s=sub.dropna(subset=[t]); b=fit(design(s,[t],bool(bas)),s.late); print(t,'n=%d OR %.2f'%(len(s),np.exp(b[1])))
print('\n== 6 late storms')
L=d[d.late==1]; print('n late',len(L))
print('single HF fix (lifecycle n_hf==1) %.3f; fixes.csv count==1 %.3f'%((L.n_hf==1).mean(),(L.nfx==1).mean()))
print('onset==peak-gust fix (catalog peak_time==t_on) %.3f'%((L.peak_time==L.t_on).mean()))
print('h_on_minp==0 share of late %.3f ; h_pg_minp==h_on_minp (peak at onset, lifecycle) %.3f'%((L.h_on_minp==0).mean(),(L.h_pg_minp==L.h_on_minp).mean()))
dm=L.msl-L.minp
print('onset msl within 2 hPa of minp (fixes.csv msl) %.3f ; tracks msl %.3f'%((dm.abs()<=2).mean(),((L.msl_on_tr-L.minp).abs()<=2).mean()))
dd=d.msl-d.minp
print('median onset msl-minp: early %.1f, late>=6h %.1f, late (all) %.1f, late h==0 %.1f'%(dd[d.h_on_minp<0].median(),dd[d.h_on_minp>=6].median(),dd[d.late==1].median(),dd[d.h_on_minp==0].median()))
print('\n== 7 joint LOSO')
from scipy.stats import rankdata
s=d.dropna(subset=['T1','T3','T4','T5','T6']).reset_index(drop=True); print('n',len(s))
def loso(preds):
    p=np.zeros(len(s)); X=design(s,preds)
    for se in s.season.unique():
        te=(s.season==se).values; b=fit(X[~te],s.late[~te]); p[te]=1/(1+np.exp(-X[te]@b))
    return p
def auc(y,p): r=rankdata(p); n1=y.sum(); n0=len(y)-n1; return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
y=s.late.values
pj=loso(['T1','T3','T4','T5','T6']); pb=loso([])
print('AUC joint %.3f basin+month %.3f'%(auc(y,pj),auc(y,pb)))
br=lambda p:((p-y)**2).mean()
# base-rate reference: LOSO basin-month cell rate
pc=np.zeros(len(s)); cell=s.basin+s.mg
for se in s.season.unique():
    te=(s.season==se).values; r=s[~te].groupby(cell[~te]).late.mean(); pc[te]=cell[te].map(r).values
print('Brier joint %.4f basin-month LOSO base %.4f skill %.3f ; vs logistic basin+month LOSO: %.3f ; vs overall mean %.3f'%(br(pj),br(pc),1-br(pj)/br(pc),1-br(pj)/br(pb),1-br(pj)/(y.mean()*(1-y.mean()))))
print('\n== T5 prevalence, denominator incl. 38 missing as 0:')
for b in(0,1):
    s=d[d.pac==b]; print('pac' if b else 'atl','%.3f'%((s.maxdeep>=1).sum()/len(s)),'defined-only %.3f'%s.T5.mean(),'missing',s.maxdeep.isna().sum())
