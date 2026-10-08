"""RA-26 follow-up: does PR 86's high-latitude lateness odds ratio (T2) survive removing tracks that touch the domain edge?
ERA5 PROXY, pipeline A, 2004-05..2025-26, non-TC, Atlantic. Plan: PREREGISTRATION.md. Committed files only.
usage: python3 late_edge.py -> results/late_edge.txt, results/tests.csv
"""
import os, importlib.util, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); ERA = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(HERE, "results"); os.makedirs(RES, exist_ok=True)
spec = importlib.util.spec_from_file_location("lw", os.path.join(ERA, "late_wind/late_wind.py"))
lw = importlib.util.module_from_spec(spec); spec.loader.exec_module(lw)
M = lw.load(); A = M[M.basin == "atl"].copy()
T = pd.read_csv(os.path.join(ERA, "hf_history/results/era5_hf_catalog_tracks.csv")); T = T[T.track.isin(A.track)].copy()
T["lon10"] = np.where(T.lon > 180, T.lon - 360, T.lon)
T["out"] = T.basin.isna()
T["touch"] = T.out & (T.lat >= 50)
T["buf"] = T.touch | ((~T.out) & ((T.lat > 64) | ((T.lon10 >= 7) & (T.lon10 <= 10))))
g = T.groupby("track")[["touch", "buf"]].max()
A = A.merge(g, left_on="track", right_index=True)
A["T6"] = (A.speed_kt - A.speed_kt.mean()) / A.speed_kt.std()
rng = np.random.default_rng(7); NB = 2000
out = []
def say(s=""): out.append(s); print(s)

def fit(d, y, cols, boot=False):
    X = sm.add_constant(pd.concat([d[[c for c in cols if c != "mg"]].astype(float)] + ([pd.get_dummies(d.mg, drop_first=True, dtype=float)] if "mg" in cols else []), axis=1))
    m = sm.Logit(d[y], X).fit(disp=0, cov_type="cluster", cov_kwds={"groups": d.season})
    return m
def row(name, d, y="late", cols=("T2",)):
    cols = list(cols); d = d.dropna(subset=[y] + [c for c in cols if c != "mg"])
    m = fit(d, y, cols); b, se, p = m.params["T2"], m.bse["T2"], m.pvalues["T2"]
    seas = d.season.unique(); idx = {s: d[d.season == s] for s in seas}; bs = []
    for _ in range(NB):
        sm_ = pd.concat([idx[s] for s in rng.choice(seas, len(seas))])
        try: bs.append(fit(sm_, y, cols).params["T2"])
        except Exception: pass
    lo, hi = np.percentile(bs, [2.5, 97.5])
    hi_sh = d[(d.T2 == 1)][y].mean(); lo_sh = d[d.T2 == 0][y].mean()
    return dict(test=name, n=len(d), n_hi=int((d.T2 == 1).sum()), late_hi=hi_sh, late_lo=lo_sh, OR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se),
                boot_lo=np.exp(lo), boot_hi=np.exp(hi), p=p, se=se, mde80=np.exp(2.8 * se), power_or2=stats.norm.cdf(np.log(2) / se - 1.96))
full = row("check: full Atlantic sample (RA-7 T2, expect 3.10)", A)
R = [
 row("E1 T2, not touching (primary)", A[~A.touch]),
 row("E2 T2, touching", A[A.touch]),
]
# E3 interaction
A["TX"] = A.T2 * A.touch.astype(float); A["touchf"] = A.touch.astype(float)
X = sm.add_constant(A[["T2", "touchf", "TX"]]); m3 = sm.Logit(A.late, X).fit(disp=0, cov_type="cluster", cov_kwds={"groups": A.season})
e3 = dict(test="E3 interaction T2 x touch (log-OR difference, touch minus not)", n=len(A), n_hi=int(A.T2.sum()), late_hi=np.nan, late_lo=np.nan, OR=np.exp(m3.params["TX"]), lo=np.exp(m3.params["TX"] - 1.96 * m3.bse["TX"]), hi=np.exp(m3.params["TX"] + 1.96 * m3.bse["TX"]),
          boot_lo=np.nan, boot_hi=np.nan, p=m3.pvalues["TX"], se=m3.bse["TX"], mde80=np.exp(2.8 * m3.bse["TX"]), power_or2=stats.norm.cdf(np.log(2) / m3.bse["TX"] - 1.96))
R.append(e3)
N = A[~A.touch]
R += [row("E4 T2, not touching, adjusted (marginal, speed, month)", N, cols=("T2", "T3", "T6", "mg")),
      row("E5 T2, not touching, LATE6", N, y="late6"),
      row("E6 T2, not touching, terrain-flagged removed", N[~N.terrain.astype(bool)]),
      row("E7 T2, not touching (buffered definition)", A[~A.buf])]
D = pd.DataFrame(R)
ps = D.p.values; o = np.argsort(ps); m = len(ps); q = np.empty(m); run = 1.0
for r, i in list(enumerate(o, 1))[::-1]:
    run = min(run, ps[i] * m / r); q[i] = run
D["q"] = q
D.to_csv(os.path.join(RES, "tests.csv"), index=False)
say("RA-26 follow-up: domain-edge artefact check. ERA5 PROXY, pipeline A, Atlantic non-TC events 2004-05..2025-26. Terrain-flagged fixes kept unless stated.")
say(f"Atlantic events {len(A)}; T2 (onset max >=60N) {int(A.T2.sum())}; touching north/east edge {int(A.touch.sum())} (T2: {int(A[A.T2==1].touch.sum())}); buffered-touching {int(A.buf.sum())} (T2: {int(A[A.T2==1].buf.sum())})")
f = full; say(f"CHECK full sample: OR {f['OR']:.2f} [{f['lo']:.2f},{f['hi']:.2f}], late share hi/lo {f['late_hi']:.3f}/{f['late_lo']:.3f}, n {f['n']}")
say("")
say(f"{'test':62s} {'n':>5s} {'nhi':>4s} {'late hi':>8s} {'late lo':>8s} {'OR':>5s} {'95% CI (CR1)':>14s} {'boot CI':>13s} {'p':>7s} {'q':>6s} {'MDE80':>6s} {'pow2':>5s}")
for _, r in D.iterrows():
    say(f"{r.test:62s} {int(r.n):5d} {int(r.n_hi):4d} {r.late_hi:8.3f} {r.late_lo:8.3f} {r.OR:5.2f} [{r.lo:5.2f},{r.hi:6.2f}] [{r.boot_lo:5.2f},{r.boot_hi:5.2f}] {r.p:7.4f} {r.q:6.3f} {r.mde80:6.2f} {r.power_or2:5.2f}")
say(f"tests passing BH q<0.05: {(q<0.05).sum()} of {len(q)}")
e1 = D.iloc[0]
rule = ("PERSISTS" if (e1.OR >= 1.76 and e1.lo > 1) else
        "SHRINKS" if e1.lo > 1 else
        "VANISHES (edge artefact)" if (e1.lo <= 1 and e1.hi < 3.10 and e1.power_or2 >= 0.8) else "CAN'T TELL")
say(f"decision rule on E1: {rule}")
open(os.path.join(RES, "late_edge.txt"), "w").write("\n".join(out) + "\n")
