"""Planted-effect check of the RA-18 decision rule. ERA5 PROXY, pipeline A. Plan: PREREGISTRATION.md (Power and limits).

usage: python3 -I power_m.py REPO OUTDIR [NSIM] [NBOOT]
Outcomes are simulated from the fitted primary outcome model with the direct pattern coefficient set so that the true
mediated fraction M is 0, 0.3, 0.5 or 0.8 (ingredient coefficients and the pattern-to-ingredient shifts kept); the
verdict is read from a season-block bootstrap interval of beta_tot (an approximation to the permutation test) and of M.
"""
import os, sys, json
from multiprocessing import Pool

import numpy as np
import pandas as pd

ROOT, OUT = sys.argv[1], sys.argv[2]
NSIM = int(sys.argv[3]) if len(sys.argv) > 3 else 60
NBOOT = int(sys.argv[4]) if len(sys.argv) > 4 else 400
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import part2 as P  # noqa: E402  (its module-level code only defines functions)

TARGETS = [0.0, 0.3, 0.5, 0.8]


def planted(d, target):
    M, x, Z, y = d["M"], d["x"], d["Z"], d["y"]
    k = Z.shape[1]
    a = np.linalg.lstsq(np.column_stack([M, x]), Z, rcond=None)[0][-1]
    Xo = np.column_stack([M, Z, x])
    bo = P.logit_fit(Xo, y)
    tot0 = P.logit_fit(np.column_stack([M, x]), y)[-1]
    nm = M.shape[1]

    def truth(b):
        p0 = P.sig(Xo @ b)
        Xi = Xo.copy(); Xi[:, nm:nm + k] += a[None, :]
        Xd = Xo.copy(); Xd[:, -1] += 1
        return (P.sig(Xi @ b) - p0).mean(), (P.sig(Xd @ b) - p0).mean()
    b = bo.copy()
    if target == 0:
        b[nm:nm + k] = 0.0; b[-1] = tot0
    else:
        nie, _ = truth(b)
        lo, hi = 0.0, 4.0
        for _ in range(60):
            b[-1] = (lo + hi) / 2
            nie, nde = truth(b)
            m = nie / (nie + nde)
            lo, hi = ((b[-1], hi) if m > target else (lo, b[-1]))  # NDE rises with the direct coefficient, so M falls as b grows
    nie, nde = truth(b)
    return b, nie / (nie + nde), Xo


def one(args):
    basin, d, target, s = args
    rng = np.random.default_rng(900 + s)
    b, mtrue, Xo = planted(d, target)
    y = (rng.random(len(d["y"])) < P.sig(Xo @ b)).astype(float)
    sea = d["sea"]; seas = np.unique(sea); idx = {q: np.where(sea == q)[0] for q in seas}
    M, x, Z = d["M"], d["x"], d["Z"]
    pt = P.one_pass(M, x, Z, y, parts=False)
    bt, bm = [], []
    for _ in range(NBOOT):
        ii = np.concatenate([idx[q] for q in rng.choice(seas, len(seas))])
        r = P.one_pass(M[ii], x[ii], Z[ii], y[ii], parts=False)
        bt.append(r["beta_tot"]); bm.append(r["M"])
    bt, bm = np.array(bt), np.array(bm)
    ok = bt >= 0.01
    sig_tot = np.percentile(bt, 2.5) > 0
    lo, hi = np.nanpercentile(np.where(ok, bm, np.nan), [2.5, 97.5])
    m = pt["M"]
    if not sig_tot or ok.mean() < 0.95:
        v = "unresolved_total"
    elif m >= 0.5 and lo > 1 / 3:
        v = "mostly"
    elif m < 1 / 3 and hi < 0.5:
        v = "little"
    else:
        v = "partly"
    return dict(basin=basin, target=target, true_M=mtrue, est=m, lo=lo, hi=hi, verdict=v, sig_tot=bool(sig_tot))


if __name__ == "__main__":
    D, _ = P.build(ROOT)
    jobs = [(b, D[b], t, s) for b in D for t in TARGETS for s in range(NSIM)]
    with Pool(4) as pool:
        res = pool.map(one, jobs, chunksize=4)
    R = pd.DataFrame(res)
    os.makedirs(OUT, exist_ok=True)
    R.to_csv(os.path.join(OUT, "power_m_sims.csv"), index=False)
    S = (R.groupby(["basin", "target", "verdict"]).size() / R.groupby(["basin", "target"]).size()).rename("share").reset_index()
    S = S.merge(R.groupby(["basin", "target"]).agg(true_M=("true_M", "mean"), mean_est=("est", "mean"), mean_width=("hi", lambda v: np.nan), n=("est", "size")).reset_index()[["basin", "target", "true_M", "mean_est", "n"]], on=["basin", "target"])
    S.to_csv(os.path.join(OUT, "power_m.csv"), index=False)
    print(S.pivot_table(index=["basin", "target", "true_M"], columns="verdict", values="share", fill_value=0).round(2))
