"""Weekly bomb counts and a small Poisson GLM (numpy only). ERA5 proxy, pipeline A. See PREREGISTRATION.md."""
import os, datetime as dt
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.dirname(HERE)
def weekly(basin):
    fx = pd.read_csv(os.path.join(ERA, "intensity/results/fixes_2004.csv.gz"))
    fx = fx[(fx.basin == basin) & (fx.ndr24 >= 1) & (fx.lat >= 30)].copy()
    fx["t"] = pd.to_datetime(fx.time.astype(str), format="%Y%m%d%H")
    first = fx.sort_values("t").groupby("track").first().reset_index()
    cnt = {}
    for t in first.t:
        d = t.date(); s = d.year if d.month >= 10 else d.year - 1
        k = (d - dt.date(s, 10, 1)).days // 7
        cnt[(s, k)] = cnt.get((s, k), 0) + 1
    ix = pd.read_csv(os.path.join(ERA, f"hemispheric/results/oos_index_{basin}.csv"))
    ix = ix[ix.season >= 2004].copy().reset_index(drop=True)
    ix["bomb"] = [cnt.get((s, k), 0) for s, k in zip(ix.season, ix.week)]
    ix["month"] = [(dt.date(s, 10, 1) + dt.timedelta(days=7 * k + 3)).month for s, k in zip(ix.season, ix.week)]
    ix["prev"] = ix.groupby("season").bomb.shift(1).fillna(ix.bomb.mean())
    return ix
def design(df, extra=()):
    M = pd.get_dummies(df.month, dtype=float).values
    cols = [M, df[["prev"]].values.astype(float)] + [np.asarray(e, float).reshape(len(df), -1) for e in extra]
    return np.hstack(cols)
def poisson_fit(X, y, ridge=0.0, pen_from=None, iters=50):
    n, p = X.shape; b = np.zeros(p); b[:] = 0; mu = np.full(n, max(y.mean(), 1e-3)); eta = np.log(mu)
    P = np.zeros(p)
    if pen_from is not None: P[pen_from:] = ridge
    else: P[:] = ridge * 0
    for _ in range(iters):
        mu = np.exp(eta); z = eta + (y - mu) / mu; W = mu
        A = X.T @ (W[:, None] * X) + np.diag(P) + 1e-8 * np.eye(p)
        bn = np.linalg.solve(A, X.T @ (W * z)); 
        if np.max(np.abs(bn - b)) < 1e-8: b = bn; break
        b = bn; eta = np.clip(X @ b, -20, 8)
    return b
def dev(y, mu):
    mu = np.maximum(mu, 1e-9); t = np.where(y > 0, y * np.log(np.maximum(y, 1e-9) / mu), 0.0)
    return 2 * np.sum(t - (y - mu))
