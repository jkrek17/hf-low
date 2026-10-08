"""Amendment 1: do moisture/latent-heating proxies explain the shortfall? Pipeline A, ERA5 proxy. usage: moisture.py REPO OUTDIR"""
import os, sys, warnings
import numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
warnings.filterwarnings("ignore")
root, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(root, "research/era5/intensity"))
import model as M
res = os.path.join(root, "research/era5/intensity/results")
D = M.load(os.path.join(res, "fixes_2004.csv.gz"), os.path.join(res, "env_2004.csv.gz"))
D["hf24"] = D.hf24.astype(int)
T1 = pd.read_csv(os.path.join(root, "research/era5/intensity_extra/results/features_tier1.csv.gz"), dtype={"time": str})
D = D.merge(T1, on=["track", "time"], how="inner")
full = M.SETS["full"]; moist = ["tcwv", "flux", "airsea", "ivt500", "ivtmax", "precip6"]
base = [c for c in full]; withm = base + [c for c in moist if c not in base]
seasons = np.sort(D.season.unique()); rng = np.random.default_rng(20261011)
B = rng.integers(0, len(seasons), (1000, len(seasons)))
L = []
def w(s): L.append(s); print(s)
def loso(df, cols, y, kind):
    p = np.zeros(len(df))
    for s in seasons:
        te = (df.season == s).values
        X = df[cols].values.astype(float); pr = M.Prep(X[~te])
        if kind == "lin": p[te] = LinearRegression().fit(pr(X[~te]), y[~te]).predict(pr(X[te]))
        else: p[te] = M.fit(pr(X[~te]), y[~te]).predict_proba(pr(X[te]))[:, 1]
    return p
def seas_sum(v): return np.array([v[(sub.season == s).values].sum() for s in seasons])
w(f"Pipeline A, ERA5 proxy; tier-1 sample {len(D)} fixes, {len(seasons)} seasons")
# H4a: deepening rate
sub = D[(D.cls >= 0)].reset_index(drop=True); y = sub.ndr24.values
e0 = (y - loso(sub, base, y, "lin")) ** 2; e1 = (y - loso(sub, withm, y, "lin")) ** 2
tot = ((y - y.mean()) ** 2); s0, s1, st = seas_sum(e0), seas_sum(e1), seas_sum(tot)
g = (s0.sum() - s1.sum()) / st.sum(); gb = (s0[B].sum(1) - s1[B].sum(1)) / st[B].sum(1)
w(f"H4a R2 for 24 h deepening rate: base {1 - s0.sum()/st.sum():.4f}, with moisture {1 - s1.sum()/st.sum():.4f}, gain {g:+.4f} [{np.percentile(gb,5):+.4f}, {np.percentile(gb,95):+.4f}] (n {len(sub)})")
# H4b: hf24 given realised deepening
sub["ndr2"] = sub.ndr24 ** 2; yy = sub.hf24.values
pc = np.zeros(len(sub))
for s in seasons:
    te = (sub.season == s).values; trn = sub[~te]
    for (b, mo), gg in sub[te].groupby(["basin", "month"]):
        t_ = trn[(trn.basin == b) & (trn.month == mo)].hf24; pc[gg.index] = (t_.sum() + .5) / (len(t_) + 1)
ob = ((loso(sub, base + ["ndr24", "ndr2"], yy, "log") - yy) ** 2); om = ((loso(sub, withm + ["ndr24", "ndr2"], yy, "log") - yy) ** 2)
cs = seas_sum((pc - yy) ** 2); sb, sm = seas_sum(ob), seas_sum(om)
g2 = (sb.sum() - sm.sum()) / cs.sum(); g2b = (sb[B].sum(1) - sm[B].sum(1)) / cs[B].sum(1)
w(f"H4b hf24 BSS given realised deepening (oracle rows): base {1 - sb.sum()/cs.sum():+.4f}, with moisture {1 - sm.sum()/cs.sum():+.4f}, gain {g2:+.4f} [{np.percentile(g2b,5):+.4f}, {np.percentile(g2b,95):+.4f}]")
# plain forecast gain, same rows
fb = ((loso(sub, base, yy, "log") - yy) ** 2); fm = ((loso(sub, withm, yy, "log") - yy) ** 2)
sfb, sfm = seas_sum(fb), seas_sum(fm)
g3 = (sfb.sum() - sfm.sum()) / cs.sum(); g3b = (sfb[B].sum(1) - sfm[B].sum(1)) / cs[B].sum(1)
w(f"     forecast (no oracle) BSS gain from moisture, same rows: {g3:+.4f} [{np.percentile(g3b,5):+.4f}, {np.percentile(g3b,95):+.4f}]")
w(f"H4 verdict: deepening gain >= 0.01 {'yes' if g >= 0.01 else 'no'}; wind-given-deepening gain < 0.005 {'yes' if g2 < 0.005 else 'no'}")
# H5: misses vs hits among rapid deepeners (cut = count-matched via training quantile)
pf = np.zeros(len(sub)); cut = np.zeros(len(sub))
for s in seasons:
    te = (sub.season == s).values; X = sub[base].values.astype(float); pr = M.Prep(X[~te])
    m = M.fit(pr(X[~te]), yy[~te]); pf[te] = m.predict_proba(pr(X[te]))[:, 1]
    cut[te] = np.quantile(m.predict_proba(pr(X[~te]))[:, 1], 1 - yy[~te].mean())
rd = (sub.ndr24 >= 1.0).values & (yy == 1)          # rapid deepeners that did reach HF within 24 h
hit = (pf >= cut)[rd]; r = sub[rd].reset_index(drop=True); r["miss"] = (~hit).astype(int)
w(f"H5 rapid deepeners reaching HF: {len(r)} fixes, misses {r.miss.sum()} ({r.miss.mean():.3f})")
pa = np.zeros(len(r)); 
for s in seasons:
    te = (r.season == s).values
    if te.sum() == 0 or r.miss[~te].nunique() < 2: pa[te] = r.miss[~te].mean(); continue
    X = r[moist].values.astype(float); pr = M.Prep(X[~te]); pa[te] = M.fit(pr(X[~te]), r.miss.values[~te]).predict_proba(pr(X[te]))[:, 1]
auc = M.roc_auc_score(r.miss, pa)
ab = []
for b in B:
    idx = np.concatenate([np.where(r.season.values == seasons[i])[0] for i in b])
    if r.miss.values[idx].min() != r.miss.values[idx].max(): ab.append(M.roc_auc_score(r.miss.values[idx], pa[idx]))
w(f"H5 LOSO AUC of the six moisture proxies for miss vs hit: {auc:.3f} [{np.percentile(ab,5):.3f}, {np.percentile(ab,95):.3f}]")
for c in moist:
    a_, b_ = r[r.miss == 1][c], r[r.miss == 0][c]
    w(f"   {c:7s} mean misses {a_.mean():8.3f} hits {b_.mean():8.3f} std diff {(a_.mean()-b_.mean())/r[c].std():+.2f}")
open(os.path.join(out, "moisture.txt"), "w").write("\n".join(L) + "\n")
