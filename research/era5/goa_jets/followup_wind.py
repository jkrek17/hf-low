"""Lead 1 (PREREGISTRATION_followup.md): is the Cook Inlet sea-minus-air result a gust-formula artefact? ERA5 PROXY, 10 m MEAN wind instead of gust.

Same Z3 GAP episodes as results/episodes.csv. WMAX = largest 10 m mean wind over Z3 ocean cells at the episode peak time; GF = zone-max gust / WMAX.
Writes results/followup_wind_tests.csv, followup_wind_episodes.csv, followup_wind.txt. usage: followup_wind.py
"""
import os, sys
import numpy as np, pandas as pd
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import features as F

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
rng = np.random.default_rng(16)
out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
E = pd.read_csv(os.path.join(RES, "episodes.csv"), parse_dates=["time"])
E = E[(E.zone == "Z3") & (E.regime == "GAP")].dropna(subset=["STAB", "DEPTH", "INTERIOR_DP"]).reset_index(drop=True)
W = F.WORK; m = F.ZM["Z3"]
wm = []
for t in E.time:
    z = np.load(os.path.join(W, t.strftime("%Y%m%d%H") + ".npz"))
    wm.append(float(np.hypot(z["u10"], z["v10"])[m].max() * F.KT))
E["WMAX"] = wm; E["GF"] = E.gust_kt / E.WMAX
order = np.lexsort((-E.gust_kt.values, -E.WMAX.values))        # WMAX descending, ties by gust
E["wind_strong"] = 0; E.loc[order[:13], "wind_strong"] = 1
for c in ("STAB", "INTERIOR_DP", "DEPTH"): E[c + "_z"] = (E[c] - E[c].mean()) / E[c].std()
E.to_csv(os.path.join(RES, "followup_wind_episodes.csv"), index=False)
P(f"Z3 GAP episodes: {len(E)}; gust-strong {int(E.strong_ep.sum())}; wind-strong (top 13 by WMAX) {int(E.wind_strong.sum())}; overlap {int(((E.strong_ep == 1) & (E.wind_strong == 1)).sum())}")
P(f"WMAX kt: median {E.WMAX.median():.1f}, range {E.WMAX.min():.1f}-{E.WMAX.max():.1f}; wind-strong cut (13th largest) {np.sort(E.WMAX.values)[-13]:.1f} kt; GF median {E.GF.median():.3f} IQR {E.GF.quantile(.25):.3f}-{E.GF.quantile(.75):.3f}")


def logit(X, y, lam=0.5):
    X1 = np.column_stack([np.ones(len(y)), X])
    fn = lambda b: np.sum(np.logaddexp(0, X1 @ b) - y * (X1 @ b)) + lam * np.sum(b[1:] ** 2)
    return minimize(fn, np.zeros(X1.shape[1]), method="BFGS").x


def ols(X, y):
    X1 = np.column_stack([np.ones(len(y)), X]); return np.linalg.lstsq(X1, y, rcond=None)[0]


def run(name, fac, kind, ycol):
    cols = [fac + "_z", "DEPTH_z"] if fac != "DEPTH" else ["DEPTH_z"]
    def est(d):
        X = d[cols].values.astype(float); y = d[ycol].values.astype(float)
        return (logit(X, y)[1] if kind == "logit" else ols(X, y)[1])
    e0 = est(E); g = {s: x for s, x in E.groupby("season")}; ss = list(g); b = []
    for _ in range(2000):
        d = pd.concat([g[s] for s in rng.choice(ss, len(ss))])
        if kind == "logit" and d[ycol].nunique() < 2: continue
        b.append(est(d))
    b = np.array(b); p = max(2 * min((b <= 0).mean(), (b >= 0).mean()), 1 / 2000)
    if kind == "logit":
        r = dict(test=name, factor=fac, outcome=ycol, effect="OR per SD", est=np.exp(e0), lo=np.exp(np.percentile(b, 2.5)), hi=np.exp(np.percentile(b, 97.5)), p=p)
    else:
        sdy = E[ycol].std()
        r = dict(test=name, factor=fac, outcome=ycol, effect="change in outcome per SD", est=e0, lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5), p=p, outcome_sd=sdy)
    return r


R = [run("1", "STAB", "logit", "wind_strong"), run("2", "INTERIOR_DP", "logit", "wind_strong"), run("3", "STAB", "ols", "WMAX"),
     run("4", "INTERIOR_DP", "ols", "WMAX"), run("5", "STAB", "ols", "GF")]
T = pd.DataFrame(R)
p = T.p.values; o = np.argsort(p); q = np.empty(len(p)); mm = len(p)
q[o] = np.minimum.accumulate((p[o] * mm / (np.arange(mm) + 1))[::-1])[::-1]; T["q"] = np.minimum(q, 1)
T.to_csv(os.path.join(RES, "followup_wind_tests.csv"), index=False)
P(T.round(3).to_string(index=False))
# reference: the same STAB test on gust-strong, replicated here so the comparison sits in one place
ref = run("ref", "STAB", "logit", "strong_ep"); P("\nreference (gust-strong, same code):", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in ref.items()})
t1, t3, t5 = T.iloc[0], T.iloc[2], T.iloc[4]
if t1.est >= 1.5 and t1.lo > 1 and t3.lo > 0: verdict = "NOT a gust artefact (test 1 OR>=1.5 excluding 1 and test 3 slope>0 excluding 0)"
elif t5.lo > 0 or t5.hi < 0: verdict = "PARTLY the gust formula (test 5 interval excludes 0)"
else: verdict = "CANNOT TELL"
P("\nPRE-REGISTERED READING:", verdict)
open(os.path.join(RES, "followup_wind.txt"), "w").write("\n".join(out)); print("\n".join(out))
