"""Part B: conversion effect, depth reached versus gust at a given depth. usage: python3 -I run_b.py OUTDIR [NBOOT]. Plan: PREREGISTRATION.md (23f1b5d)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import placelib as P

OUT = sys.argv[1]
NBOOT = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
os.makedirs(OUT, exist_ok=True)
sys.argv = ["part1", P.REPO, OUT]
sys.path.insert(0, os.path.join(P.REPO, "research", "era5", "share_mechanism"))
import part1 as S1  # noqa: E402
C = P.C
Tbase = S1.load(-3.6)
MD = C.month_dummies()


def logit(X, y, pen, iters=50):
    b = np.zeros(X.shape[1]); b[0] = 0.0
    for _ in range(iters):
        eta = np.clip(X @ b, -30, 30); mu = 1 / (1 + np.exp(-eta))
        g = X.T @ (y - mu) - pen * b
        H = (X * (mu * (1 - mu))[:, None]).T @ X + np.diag(pen) + 1e-8 * np.eye(len(b))
        st = np.linalg.solve(H, g); b = b + st
        if np.abs(st).max() < 1e-8:
            break
    return b


def design(d, kind, cells=None):
    cols = [MD[d.row.values], d.x.values[:, None]]           # month dummies, then x (column 7)
    pen = [np.full(7, 1e-6), [0.0]]
    if kind >= 2:
        z = (d.minp.values - 965.0) / 15.0
        cols += [z[:, None], (z ** 2)[:, None]]; pen += [[1e-6], [1e-6]]
    if kind >= 3:
        cell = (np.floor(d.peak_lat.values / 5) * 100 + np.floor(d.lon_u.values / 10)).astype(int)
        U = sorted(set(cell))
        cols.append(np.array([[1.0 if c == u else 0.0 for u in U[1:]] for c in cell]) if len(U) > 1 else np.zeros((len(d), 0)))
        pen.append(np.full(max(len(U) - 1, 0), 1.0))
    return np.hstack(cols), np.concatenate([np.asarray(p, float) for p in pen]), 7


def slopes(d):
    out = []
    for k in (1, 2, 3):
        X, pen, ix = design(d, k)
        out.append(logit(X, d.hf.values.astype(float), pen)[ix])
    return np.array(out)


def main():
    rng = np.random.default_rng(20261011)
    res, rows = [], []
    for b in C.BASINS:
        t = Tbase[Tbase.basin == b]
        win = C.assign_week(t)
        d = win[win.deep == 1].copy()
        x = C.index(b)[0]
        d["x"] = x[d.row.values]
        d = d.reset_index(drop=True)
        obs = slopes(d)
        sea = d.sea.values
        S = sorted(set(sea))
        bs = np.empty((NBOOT, 3))
        for i in range(NBOOT):
            pick = rng.integers(0, len(S), len(S))
            dd = pd.concat([d[sea == S[j]] for j in pick], ignore_index=True)
            bs[i] = slopes(dd)
        att2, att3 = 1 - bs[:, 1] / bs[:, 0], 1 - bs[:, 2] / bs[:, 0]
        for nm, a, o in (("M2_minp", att2, 1 - obs[1] / obs[0]), ("M3_minp_cell", att3, 1 - obs[2] / obs[0])):
            lo, hi = np.percentile(a, [2.5, 97.5])
            p = float(min(1, 2 * min((a <= 0).mean(), (a >= 0).mean())) or 1.0 / NBOOT)
            rows.append(dict(basin=b, measure=nm, attenuation=o, lo=lo, hi=hi, p=p))
        res.append(dict(basin=b, n_deepening=len(d), n_hf=int(d.hf.sum()),
                        b_M1=obs[0], b_M2=obs[1], b_M3=obs[2],
                        b_M1_ci=np.percentile(bs[:, 0], [2.5, 97.5]).tolist(), b_M2_ci=np.percentile(bs[:, 1], [2.5, 97.5]).tolist(), b_M3_ci=np.percentile(bs[:, 2], [2.5, 97.5]).tolist()))
        print(b, res[-1], flush=True)
    A = pd.DataFrame(rows)
    A["q"] = C.bh(A.p.values)
    A.to_csv(os.path.join(OUT, "partB_attenuation.csv"), index=False)
    json.dump(res, open(os.path.join(OUT, "partB_models.json"), "w"), indent=1, default=float)
    print(A)


if __name__ == "__main__":
    main()
