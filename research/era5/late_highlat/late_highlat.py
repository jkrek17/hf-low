"""RA-26: is late HF onset in the high-latitude Atlantic barrier flow behind the low?
ERA5 PROXY, pipeline A (research/era5/hf_history), seasons 2004-05..2025-26, non-TC events. Plan: PREREGISTRATION.md.
Reads committed files only (no pull).   usage: python3 late_highlat.py  -> results/late_highlat.txt, results/tests.csv
"""
import os, importlib.util
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(HERE, "results"); os.makedirs(RES, exist_ok=True)
SEED, NB = 7, 2000
spec = importlib.util.spec_from_file_location("lw", os.path.join(ERA, "late_wind/late_wind.py"))
lw = importlib.util.module_from_spec(spec); spec.loader.exec_module(lw)

M = lw.load()
M["barrier"] = None
F = pd.read_csv(os.path.join(ERA, "hf_structure/results/fixes.csv")).sort_values(["track", "time"])
on = F.groupby("track").first()
M = M.merge(on[["gmax_coast_km"]].rename(columns={"gmax_coast_km": "coast"}), left_on="track", right_index=True, how="left")
# peak-gust HF fix (first of the maximum g800)
pk = F.loc[F.groupby("track").g800.idxmax()].set_index("track")[["gmax_r", "gmax_coast_km", "gmax_lat", "terrain"]]
pk.columns = ["pk_r", "pk_coast", "pk_lat", "pk_terrain"]
M = M.merge(pk, left_on="track", right_index=True, how="left")
# wind direction at the maximum
G = pd.read_csv(os.path.join(ERA, "highlat/gustloc_fixes.csv"))[["track", "time", "max_wdir"]]
M = M.merge(G.rename(columns={"time_x": "time"}), on=["track", "time"], how="left")
M["barrier"] = ((M.coast <= 300) & (M.gmax_r > 400)).astype(float)
M["near"] = (M.coast <= 300).astype(float)
M["far"] = (M.gmax_r > 400).astype(float)
M["ctr"] = (M.gmax_r <= 400).astype(float)
M["north"] = np.where(M.max_wdir.isna(), np.nan, ((M.max_wdir >= 315) | (M.max_wdir <= 45)).astype(float))
M["pkbar"] = ((M.pk_coast <= 300) & (M.pk_r > 400)).astype(float)
M["mon"] = M.mg
H = M[(M.basin == "atl") & (M.T2 == 1)].copy()
S = M[(M.basin == "atl") & (M.T2 == 0)].copy()
rng = np.random.default_rng(SEED)

out = []
def say(s=""): out.append(s); print(s)

def diff_boot(d, y, lab):
    d = d.dropna(subset=[y]); a = d[d[lab] == 1]; b = d[d[lab] == 0]
    sh = lambda x: x[y].mean()
    obs = sh(a) - sh(b)
    seas = d.season.unique(); idx = {s: d[d.season == s] for s in seas}
    bd = []
    for _ in range(NB):
        smp = pd.concat([idx[s] for s in rng.choice(seas, len(seas))])
        x, z = smp[smp[lab] == 1][y], smp[smp[lab] == 0][y]
        if len(x) and len(z): bd.append(x.mean() - z.mean())
    bd = np.array(bd)
    p = min(1.0, 2 * (min((bd <= 0).sum(), (bd >= 0).sum()) + 1) / (len(bd) + 1))
    lo, hi = np.percentile(bd, [2.5, 97.5]); se = bd.std()
    pw = stats.norm.cdf(0.25 / se - 1.96) if se > 0 else np.nan
    return dict(n_late=len(a), n_non=len(b), late=sh(a), non=sh(b), diff=obs, lo=lo, hi=hi, p=p, se=se, power_025=pw, mde80=2.8 * se)

tests = []
def add(name, d, y, lab="late"):
    r = diff_boot(d, y, lab); r["test"] = name; tests.append(r)

H["late6"] = H.late6.astype(int)
add("P1 barrier share, onset fix, LATE", H, "barrier")
add("S1 coast <=300 km", H, "near")
add("S2 centre distance >400 km", H, "far")
add("S3 northerly wind at max", H, "north")
add("S4 barrier, LATE6", H, "barrier", "late6")
add("S5 barrier, terrain-flagged removed", H[~H.terrain.astype(bool)], "barrier")
add("S6 barrier at peak-gust fix", H, "pkbar")
add("S7 maximum within 400 km of centre", H, "ctr")
add("S8 barrier, Atlantic south of 60N (contrast)", S, "barrier")
T = pd.DataFrame(tests)
# S9 adjusted logistic, CR1 by season
import statsmodels.api as sm
X = H[["barrier", "T3", "T6"]].copy()
X = pd.concat([X, pd.get_dummies(H.mon, drop_first=True, dtype=float)], axis=1); X = sm.add_constant(X)
ok = X.notna().all(axis=1) & H.late.notna()
fit = sm.Logit(H.late[ok], X[ok]).fit(disp=0, cov_type="cluster", cov_kwds={"groups": H.season[ok]})
b, se_, pv = fit.params["barrier"], fit.bse["barrier"], fit.pvalues["barrier"]
orr = (np.exp(b), np.exp(b - 1.96 * se_), np.exp(b + 1.96 * se_))
# BH
ps = list(T.p) + [pv]; order = np.argsort(ps); m = len(ps); q = np.empty(m)
run = 1.0
for rnk, i in list(enumerate(order, 1))[::-1]:
    run = min(run, ps[i] * m / rnk); q[i] = run
T["q"] = q[:len(T)]
T.to_csv(os.path.join(RES, "tests.csv"), index=False)

say("RA-26 late HF onset, high-latitude Atlantic. ERA5 PROXY, pipeline A, 2004-05..2025-26, non-TC. Terrain-flagged fixes kept unless stated.")
say(f"events total {len(M)}, Atlantic {int((M.basin=='atl').sum())}, group (Atlantic, onset max >=60N) {len(H)}: late {int(H.late.sum())}, non-late {int((1-H.late).sum())}, LATE6 {int(H.late6.sum())}, terrain-flagged {int(H.terrain.sum())}; dropped {lw.INFO['dropped']}")
say(f"wind direction available for {int(H.north.notna().sum())} of {len(H)}; late with wdir {int(H[H.late==1].north.notna().sum())}")
say("")
say(f"{'test':48s} {'nL':>4s} {'nN':>4s} {'late':>6s} {'non':>6s} {'diff':>6s} {'95% CI':>15s} {'p':>7s} {'q':>7s} {'pow.25':>6s} {'MDE80':>6s}")
for _, r in T.iterrows():
    say(f"{r.test:48s} {int(r.n_late):4d} {int(r.n_non):4d} {r.late:6.3f} {r.non:6.3f} {r['diff']:6.3f} [{r.lo:6.3f},{r.hi:6.3f}] {r.p:7.4f} {r.q:7.4f} {r.power_025:6.2f} {r.mde80:6.3f}")
say(f"{'S9 adjusted OR barrier (late vs non-late)':48s} OR {orr[0]:.2f} [{orr[1]:.2f},{orr[2]:.2f}] p {pv:.4f} q {q[-1]:.4f}  n {int(ok.sum())}")
say(f"tests passing BH q<0.05: {(q<0.05).sum()} of {len(q)}")
p1 = T.iloc[0]
rule = ("YES" if (p1.late >= 0.60 and p1.non <= 0.35 and p1.lo > 0) else
        "PARTLY" if p1.lo > 0 else
        "NO (well powered)" if (p1.hi < 0.15 and p1.power_025 >= 0.8) else "CAN'T TELL")
say(f"decision rule on P1: {rule}")
open(os.path.join(RES, "late_highlat.txt"), "w").write("\n".join(out) + "\n")
