"""RA-7 post hoc diagnostics (not in PREREGISTRATION.md; labelled post hoc in the README).
(a) joint model with season-clustered intervals; (b) how far above the minimum is the pressure at onset for late storms.
usage: python3 posthoc.py -> results/posthoc.txt
"""
import os
import numpy as np, pandas as pd
from scipy import stats
import late_wind as lw

HERE = lw.HERE


def joint_cr(d, preds, basin=True):
    cols = [np.ones(len(d))] + ([d.pac.values] if basin else []) + [d[k].values.astype(float) for k in preds] + \
           [(d.mg == "DJF").values.astype(float), (d.mg == "MAM").values.astype(float)]
    X = np.column_stack(cols)
    y = d.late.values.astype(float)
    b, _ = lw.irls(X, y)
    p = 1 / (1 + np.exp(-(X @ b)))
    A = np.linalg.inv((X.T * (p * (1 - p))) @ X)
    s = X * (y - p)[:, None]
    g = d.season.values
    meat = sum(np.outer(s[g == k].sum(0), s[g == k].sum(0)) for k in np.unique(g))
    G, n, kk = len(np.unique(g)), len(d), X.shape[1]
    V = G / (G - 1) * (n - 1) / (n - kk) * A @ meat @ A
    se = np.sqrt(np.diag(V))
    tc = stats.t.ppf(0.975, G - 1)
    names = ["const"] + (["pac"] if basin else []) + preds + ["DJF", "MAM"]
    return pd.DataFrame(dict(term=names, OR=np.exp(b), lo=np.exp(b - tc * se), hi=np.exp(b + tc * se),
                             p=2 * stats.t.sf(abs(b / se), G - 1)))


def main():
    M = lw.load()
    o = []
    w = o.append
    w("RA-7 post hoc diagnostics (ERA5 PROXY, pipeline A). Not pre-registered.")
    d = M.dropna(subset=["T5"])
    jp = ["T1", "T3", "T4", "T5", "T6"]
    w("== (a) joint model, all events, season-clustered 95% intervals ==")
    w(joint_cr(d, jp).round(3).to_string(index=False))
    da = d[d.basin == "atl"]
    w("== (a) Atlantic only, T2 added ==")
    w(joint_cr(da, jp + ["T2"], basin=False).round(3).to_string(index=False))
    w("")
    w("== (b) pressure at onset minus minimum pressure, hPa ==")
    M["gap"] = M.msl - M.minp
    for name, s in (("early (onset before minimum)", M[M.h_on_minp < 0]), ("late, onset at the minimum fix (h=0)", M[M.h_on_minp == 0]),
                    ("late by >= 6 h", M[M.h_on_minp >= 6]), ("late by >= 12 h", M[M.h_on_minp >= 12])):
        w(f"{name}: n={len(s)}, median {s.gap.median():.1f}, quartiles {s.gap.quantile(.25):.1f}/{s.gap.quantile(.75):.1f}, share within 2 hPa of the minimum {np.mean(s.gap <= 2):.3f}")
    L = M[M.late == 1]
    w(f"all late storms: share whose pressure at onset is within 2 hPa of the minimum {np.mean(L.gap <= 2):.3f}; within 5 hPa {np.mean(L.gap <= 5):.3f}")
    w(f"minimum depth (hPa) median: late {L.minp.median():.1f}, early {M[M.late == 0].minp.median():.1f}")
    w("late share by minimum-pressure tercile: " + ", ".join(f"{k} {v:.3f}" for k, v in M.groupby(pd.qcut(M.minp, 3), observed=True).late.mean().items()))
    open(os.path.join(lw.RES, "posthoc.txt"), "w").write("\n".join(o) + "\n")
    print("\n".join(o))


if __name__ == "__main__":
    main()
