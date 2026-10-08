"""RA-21 partial dependence and Friedman H^2 for M0 (logistic) and M1 (boosted trees), fitted on all 22 seasons (in-sample, descriptive; no p values).
Run only because P1 gave 'real, as predicted' (PREREG.md). Pipeline A, ERA5 proxy.

usage: pdp.py REPO_ROOT OUTDIR
"""
import os, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
root, out = sys.argv[1:3]
sys.path.insert(0, os.path.join(root, "research/era5/intensity"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "4")
import model as M
from sklearn.ensemble import HistGradientBoostingClassifier
import run_loso as R   # loads D and the constants only (main guarded)

COLS, CONT, PAIRS, HGB = R.COLS, R.CONT, R.PAIRS, R.HGB
X, y = R.X, R.y
prep = M.Prep(X)
m0 = R.lr(prep(X), y)
m1 = HistGradientBoostingClassifier(**HGB).fit(X, y)
rng = np.random.default_rng(20261010)
S = rng.choice(len(X), 5000, replace=False)
Xs = X[S]

def f0(A): return m0.decision_function(prep(A))
def f1(A): return m1.decision_function(A)

def pd1(f, j, grid):
    v = []
    for g in grid:
        A = Xs.copy(); A[:, j] = g; v.append(f(A).mean())
    return np.array(v)

rows, curves = [], {}
for c in CONT:
    j = COLS.index(c)
    col = X[:, j]; col = col[~np.isnan(col)]
    grid = np.percentile(col, np.linspace(2, 98, 25))
    a, b = pd1(f0, j, grid), pd1(f1, j, grid)
    curves[c] = (grid, a, b)
    A = np.vstack([np.ones_like(grid), grid]).T
    coef = np.linalg.lstsq(A, b, rcond=None)[0]
    r2 = 1 - ((b - A @ coef) ** 2).sum() / ((b - b.mean()) ** 2).sum()
    step = np.diff(b); k = int(np.argmax(np.abs(step)))
    rows.append(dict(feature=c, grid_lo=grid[0], grid_hi=grid[-1], M0_range_logit=a.max() - a.min(), M1_range_logit=b.max() - b.min(),
                     M1_linear_R2=r2, M1_steepest_step_at=(grid[k] + grid[k + 1]) / 2, M1_steepest_step_logit=step[k],
                     M1_minus_M0_max_gap_logit=np.max(np.abs((b - b.mean()) - (a - a.mean())))))
    pd.DataFrame(dict(feature=c, x=grid, M0_logit=a, M1_logit=b)).to_csv(os.path.join(out, "pd_curves.csv"), mode="a", header=not os.path.exists(os.path.join(out, "pd_curves.csv")), index=False)
T = pd.DataFrame(rows)
T.round(4).to_csv(os.path.join(out, "pd_summary.csv"), index=False)

# Friedman H^2 on the logit scale, 800 fixes, for the five pre-specified pairs
Hs = Xs[:800]
def H2(f, a, b):
    ia, ib = COLS.index(a), COLS.index(b)
    fab = np.empty(len(Hs)); fa = np.empty(len(Hs)); fb = np.empty(len(Hs))
    for i in range(len(Hs)):
        A = Hs.copy(); A[:, ia] = Hs[i, ia]; A[:, ib] = Hs[i, ib]; fab[i] = f(A).mean()
        A = Hs.copy(); A[:, ia] = Hs[i, ia]; fa[i] = f(A).mean()
        A = Hs.copy(); A[:, ib] = Hs[i, ib]; fb[i] = f(A).mean()
    fab, fa, fb = fab - fab.mean(), fa - fa.mean(), fb - fb.mean()
    return float(((fab - fa - fb) ** 2).sum() / (fab ** 2).sum())
hrows = [dict(pair=f"{a} x {b}", H2_M0=H2(f0, a, b), H2_M1=H2(f1, a, b)) for a, b in PAIRS]
Hd = pd.DataFrame(hrows); Hd.round(4).to_csv(os.path.join(out, "pd_interactions.csv"), index=False)

lines = ["Partial dependence of M0 (logistic) and M1 (boosted trees), logit scale, fitted on all 22 seasons, 5,000 random fixes, 25 grid points (2nd-98th percentile). Descriptive.",
         T.round(3).to_string(index=False), "", "Friedman H^2 (logit scale, 800 fixes); M0 is additive on the logit scale so its H^2 is 0 up to rounding:", Hd.round(4).to_string(index=False)]
open(os.path.join(out, "pd.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
show = ["g800", "B", "VTL", "VTU", "sst", "flux", "lat", "msl", "dp12"]
fig, ax = plt.subplots(3, 3, figsize=(11, 8.5))
for a_, c in zip(ax.ravel(), show):
    g, a, b = curves[c]
    a_.plot(g, a - a.mean(), color="#888", label="logistic (PR 12)"); a_.plot(g, b - b.mean(), color="#c0392b", label="boosted trees")
    a_.set_title(c, fontsize=10); a_.set_ylabel("partial dependence (logit)", fontsize=8)
ax[0, 0].legend(fontsize=8)
fig.suptitle("Pipeline A, ERA5 proxy: P(HF within 24 h), partial dependence, all 22 seasons", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(out, "pd.png"), dpi=110)
