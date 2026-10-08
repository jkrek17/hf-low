"""RA-21 POST HOC diagnostics (not in PREREG.md; logged there). Pipeline A, ERA5 proxy.
PH1 where the gain of trees beyond smooth additive splines (M1 - M3) sits, by current g800.
PH2 g800 partial dependence extended to the 99.9th percentile (the pre-registered 2-98 grid stops at 71.1 kt, before the HF fixes).
PH3 M3b = M3 + the five pre-specified product terms, LOSO: how much of M1 - M3 do those pairs explain.
PH4 splines for the three Hart phase-space parameters together (B, VTL, VTU), and for g800 plus the three, LOSO, as a share of M3 - M0.
usage: posthoc.py REPO_ROOT WORKDIR OUTDIR
"""
import os, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
root, work, out = sys.argv[1:4]
sys.path.insert(0, os.path.join(root, "research/era5/intensity"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model as M
import run_loso as R
from joblib import Parallel, delayed
from sklearn.ensemble import HistGradientBoostingClassifier
D, X, y, seas, COLS, CONT, PAIRS, ci = R.D, R.X, R.y, R.seas, R.COLS, R.CONT, R.PAIRS, R.ci
us = np.unique(seas)
Z = np.load(os.path.join(work, "probs.npz")); P = {k: Z[k].astype(float) for k in Z.files}
B = np.random.default_rng(20261010).integers(0, len(us), (1000, len(us)))
lines = []
def w(s=""):
    lines.append(s); print(s)
def ssum(v, m=None):
    m = np.ones(len(v), bool) if m is None else m
    return np.array([v[m & (seas == s)].sum() for s in us])
bs = lambda p: (p - y) ** 2

# PH1
cur = D.g800.values
st = {"<55 kt": cur < 55, "55-71.7 kt": (cur >= 55) & (cur < 71.7), ">=71.7 kt": cur >= 71.7}
tot = ssum(bs(P["M3"])) - ssum(bs(P["M1"]))
w("PH1 (post hoc) share of the M3 -> M1 Brier reduction by current g800 (M3 = smooth additive splines, M1 = boosted trees):")
for k, m in st.items():
    d = ssum(bs(P["M3"]), m) - ssum(bs(P["M1"]), m)
    b = d[B].sum(1) / tot[B].sum(1)
    w(f"  {k:11s} share {d.sum() / tot.sum() * 100:5.1f}% [{np.percentile(b, 5) * 100:.1f}, {np.percentile(b, 95) * 100:.1f}]")
tot3 = ssum(bs(P["M0"])) - ssum(bs(P["M3"]))
w("PH1b share of the M0 -> M3 reduction (what smooth curves capture), same strata:")
for k, m in st.items():
    d = ssum(bs(P["M0"]), m) - ssum(bs(P["M3"]), m)
    w(f"  {k:11s} share {d.sum() / tot3.sum() * 100:5.1f}%")

# PH3
def fold4(s):
    tr, te = np.where(seas != s)[0], np.where(seas == s)[0]
    prep = M.Prep(X[tr]); Ztr, Zte = prep(X[tr]), prep(X[te])
    o = []
    for feats in (["B", "VTL", "VTU"], ["B", "VTL", "VTU", "g800"], ["sst", "flux", "lat", "eady", "vadv500", "msl", "jet250"]):
        idx = [ci[c] for c in feats]
        Str, Ste, _ = R.spline_block(Ztr, Zte, idx)
        m = R.lr(np.hstack([Ztr, Str]), y[tr]); o.append(m.predict_proba(np.hstack([Zte, Ste]))[:, 1])
    return te, o
def fold(s):
    tr, te = np.where(seas != s)[0], np.where(seas == s)[0]
    prep = M.Prep(X[tr]); Ztr, Zte = prep(X[tr]), prep(X[te])
    cidx = [ci[c] for c in CONT]
    Str, Ste, _ = R.spline_block(Ztr, Zte, cidx)
    ptr = np.column_stack([Ztr[:, ci[a]] * Ztr[:, ci[b]] for a, b in PAIRS]); pte = np.column_stack([Zte[:, ci[a]] * Zte[:, ci[b]] for a, b in PAIRS])
    m = R.lr(np.hstack([Ztr, Str, ptr]), y[tr])
    return te, m.predict_proba(np.hstack([Zte, Ste, pte]))[:, 1]
parts = Parallel(n_jobs=4)(delayed(fold)(s) for s in us)
p3b = np.zeros(len(y))
for te, p in parts: p3b[te] = p
P["M3b"] = p3b
clim = ssum(bs(P["clim"]))
def g(a, b):
    d = (ssum(bs(P[b])) - ssum(bs(P[a]))) / clim.sum()
    bt = (ssum(bs(P[b]))[B].sum(1) - ssum(bs(P[a]))[B].sum(1)) / clim[B].sum(1)
    return d.sum(), np.percentile(bt, 5), np.percentile(bt, 95), int((d > 0).sum())
w("")
w(f"PH3 (post hoc) M3b = M3 + the five pre-specified products. BSS: M3 {1 - ssum(bs(P['M3'])).sum() / clim.sum():+.4f}, M3b {1 - ssum(bs(P['M3b'])).sum() / clim.sum():+.4f}, M1 {1 - ssum(bs(P['M1'])).sum() / clim.sum():+.4f}")
for a, b, lab in (("M3b", "M3", "M3b - M3 (the five pairs)"), ("M1", "M3b", "M1 - M3b (what remains for trees)")):
    e, lo, hi, n = g(a, b)
    w(f"  {lab}: {e:+.4f} [{lo:+.4f}, {hi:+.4f}], seasons better {n}/22")

# PH4
clim = ssum(bs(P["clim"]))
def g(a, b):
    d = (ssum(bs(P[b])) - ssum(bs(P[a]))) / clim.sum()
    bt = (ssum(bs(P[b]))[B].sum(1) - ssum(bs(P[a]))[B].sum(1)) / clim[B].sum(1)
    return d.sum(), np.percentile(bt, 5), np.percentile(bt, 95), int((d > 0).sum())
parts4 = Parallel(n_jobs=4)(delayed(fold4)(s) for s in us)
for k, lab in enumerate(["Hart B, VTL, VTU", "Hart B, VTL, VTU + g800", "seven other predictors picked post hoc from the S9 table (sst, flux, lat, eady, vadv500, msl, jet250)"]):
    p4 = np.zeros(len(y))
    for te, o in parts4: p4[te] = o[k]
    P["PH4_%d" % k] = p4
    e, lo, hi, n = g("PH4_%d" % k, "M0")
    w(f"PH4 (post hoc) M0 + splines for {lab}: {e:+.4f} [{lo:+.4f}, {hi:+.4f}], seasons better {n}/22, share of M3 - M0 {e / g('M3', 'M0')[0] * 100:.0f}%")

# PH2
m1 = HistGradientBoostingClassifier(**R.HGB).fit(X, y)
prep = M.Prep(X); m0 = R.lr(prep(X), y)
Xs = X[np.random.default_rng(20261010).choice(len(X), 5000, replace=False)]
j = COLS.index("g800")
grid = np.percentile(X[:, j], np.linspace(0.1, 99.9, 30)); grid = np.unique(grid)
w("")
w("PH2 (post hoc) g800 partial dependence, logit, 0.1-99.9th percentile (M0 logistic, M1 trees):")
rows = []
for gv in grid:
    A = Xs.copy(); A[:, j] = gv
    rows.append((gv, m0.decision_function(prep(A)).mean(), m1.decision_function(A).mean()))
for gv, a, b in rows:
    w(f"  g800 {gv:6.1f} kt  M0 {a:+.2f}  M1 {b:+.2f}")
pd.DataFrame(rows, columns=["g800", "M0_logit", "M1_logit"]).round(3).to_csv(os.path.join(out, "posthoc_g800_pd.csv"), index=False)
open(os.path.join(out, "posthoc.txt"), "w").write("\n".join(lines) + "\n")
