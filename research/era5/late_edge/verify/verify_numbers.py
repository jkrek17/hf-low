import importlib.util, numpy as np, pandas as pd, statsmodels.formula.api as smf
ERA="/home/claude/hf-low/research/era5/"
sp=importlib.util.spec_from_file_location("lw",ERA+"late_wind/late_wind.py"); lw=importlib.util.module_from_spec(sp); sp.loader.exec_module(lw)
M=lw.load(); M=M[M.basin=="atl"].copy()
T=pd.read_csv(ERA+"hf_history/results/era5_hf_catalog_tracks.csv"); T=T[T.track.isin(M.track)]
print("basin values in tracks of atl events:",T.basin.value_counts(dropna=False).to_dict())
out=T.basin.isna()&(T.lat>=50)
touch=set(T[out].track)
buf=touch|set(T[T.basin.notna()&((T.lat>64)|((T.lon>=7)&(T.lon<=10)))].track)
M["touch"]=M.track.isin(touch).astype(int); M["tb"]=M.track.isin(buf).astype(int)
M["season"]=np.where(pd.to_datetime(M.t_on).dt.month>=6,pd.to_datetime(M.t_on).dt.year,pd.to_datetime(M.t_on).dt.year-1)
print("seasons:",M.season.nunique(), "NaN T2:",M.T2.isna().sum())
def n(d):return len(d)
print("\n== counts ==")
for lab,d in [("Atlantic",M),("T2==1",M[M.T2==1]),("TOUCH all",M[M.touch==1]),("TOUCH & T2==1",M[(M.touch==1)&(M.T2==1)]),
              ("bufTOUCH all",M[M.tb==1]),("bufTOUCH & T2==1",M[(M.tb==1)&(M.T2==1)])]: print(f"{lab}: n={n(d)}")
def lr(d,f,term="T2",lab=""):
    d=d.dropna(subset=[c for c in ["late","late6","T2","T3","T6","mg","touch"] if c in d])
    m=smf.logit(f,d).fit(disp=0,cov_type="cluster",cov_kwds={"groups":d.season.values})
    b,s=m.params[term],m.bse[term]
    print(f"{lab}: n={len(d)} clusters={d.season.nunique()} OR={np.exp(b):.3f} CI=[{np.exp(b-1.96*s):.3f},{np.exp(b+1.96*s):.3f}] p={m.pvalues[term]:.4g} se={s:.4f}")
print("\n== logistic ==")
nt=M[M.touch==0]
lr(M,"late~T2",lab="Full Atlantic")
lr(nt,"late~T2",lab="E1 nonTOUCH")
lr(M[M.touch==1],"late~T2",lab="E2 TOUCH")
lr(nt,"late6~T2",lab="E5 nonTOUCH LATE6")
lr(nt[nt.terrain==0],"late~T2",lab="E6 nonTOUCH no terrain")
lr(M[M.tb==0],"late~T2",lab="E7 nonTOUCH buffered")
lr(nt,"late~T2+T3+T6+C(mg)",lab="E4 nonTOUCH adj")
lr(M,"late~T2+touch+T2:touch",term="T2:touch",lab="E3 interaction OR")
print("terrain value counts:",M.terrain.value_counts(dropna=False).to_dict())
print("\n== late shares in nonTOUCH ==")
for v in (1,0):
    d=nt[nt.T2==v]; print(f"T2=={v}: n={len(d)} late={d.late.sum()} share={d.late.mean():.4f}")
for v in (1,0):
    d=M[M.T2==v]; print(f"(full) T2=={v}: n={len(d)} late share={d.late.mean():.4f}")
