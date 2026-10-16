"""Stage 4 audit (POST HOC, logged in PREREGISTRATION.md): are the Q2 leads carried by calendar month, a few seasons, or chance?

For each of the 15 Q2 tests: (a) the same fit with a late-winter term (January-April = 1) added, (b) even-season and odd-season fits,
(c) leave-one-season-out range, (d) a placebo: the factor shuffled among episodes of the same season (1,000 times), p = share of
shuffles with an estimate at least as large in size. ERA5 PROXY. Reads results/episodes.csv.
"""
import os
import numpy as np, pandas as pd
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
rng = np.random.default_rng(17)
E = pd.read_csv(os.path.join(RES, "episodes.csv"), parse_dates=["time"])
E["late"] = E.time.dt.month.isin([1, 2, 3, 4]).astype(float)
FACT = {"Z1": ["DEPTH", "LOW_DIST", "GRAD", "COASTPILE", "STAB"], "Z2": ["DEPTH", "LOW_DIST", "GRAD", "COASTPILE", "STAB"],
        "Z3": ["DEPTH", "LOW_DIST", "GRAD", "INTERIOR_DP", "STAB"]}


def fit(X, y, lam=0.5):
    X1 = np.column_stack([np.ones(len(y)), X])
    fn = lambda b: np.sum(np.logaddexp(0, X1 @ b) - y * (X1 @ b)) + lam * np.sum(b[1:] ** 2)
    return minimize(fn, np.zeros(X1.shape[1]), method="BFGS").x


def est(d, fac, extra=()):
    z = lambda c: (d[c] - d[c].mean()) / d[c].std()
    cols = [z(fac)] + ([] if fac == "DEPTH" else [z("DEPTH")]) + [d[c] for c in extra]
    return fit(np.column_stack(cols), d.strong_ep.values.astype(float))[1]


out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
rows = []
for zn, rg in [("Z1", "BARRIER"), ("Z2", "BARRIER"), ("Z3", "GAP")]:
    d = E[(E.zone == zn) & (E.regime == rg)].reset_index(drop=True)
    for fac in FACT[zn]:
        o = np.exp(est(d, fac)); om = np.exp(est(d, fac, ["late"]))
        ev = np.exp(est(d[d.season % 2 == 0], fac)); od = np.exp(est(d[d.season % 2 == 1], fac))
        loo = [np.exp(est(d[d.season != s], fac)) for s in sorted(d.season.unique())]
        sh = []
        for _ in range(1000):
            dd = d.copy()
            dd[fac] = dd.groupby("season")[fac].transform(lambda x: rng.permutation(x.values))
            sh.append(est(dd, fac))
        sh = np.array(sh); e0 = np.log(o)
        rows.append(dict(zone=zn, regime=rg, factor=fac, n=len(d), strong=int(d.strong_ep.sum()), OR=o, OR_month_adj=om, OR_even=ev, OR_odd=od,
                         loo_min=min(loo), loo_max=max(loo), placebo_p=float((np.abs(sh) >= abs(e0)).mean())))
A = pd.DataFrame(rows); A.to_csv(os.path.join(RES, "audit.csv"), index=False)
P("POST HOC audit of the 15 Q2 estimates (OR per SD, strong vs 54-64 kt episodes)")
P(A.round(2).to_string(index=False))
c = E.groupby("zone").apply(lambda g: g.groupby(g.time.dt.month).size()).unstack(0).fillna(0).astype(int)
P("\nepisodes per calendar month by zone"); P(c.to_string())
open(os.path.join(RES, "audit.txt"), "w").write("\n".join(out)); print("\n".join(out))
