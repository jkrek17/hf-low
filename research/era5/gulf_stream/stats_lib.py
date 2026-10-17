"""Shared pieces: winter-block bootstrap on per-winter sufficient statistics (X'X, X'y)."""
import numpy as np, pandas as pd
def design_wind(d, stations):
    """station dummies + log Vg + dT bins (4 non-reference: stable<0, 3-6, 6-9, >=9; reference 0-3)."""
    X = [(d.station == s).astype(float).values for s in stations]
    X.append(np.log(d.vg.values))
    dT = d.dT.values
    for lo, hi in [(-9, 0), (3, 6), (6, 9), (9, 99)]:
        X.append(((dT >= lo) & (dT < hi)).astype(float))
    return np.column_stack(X)
class Boot:
    """Holds per-winter X'X and X'y for one model; solves any resample of winters."""
    def __init__(self, X, y, winter):
        self.w = np.unique(winter); self.k = X.shape[1]
        self.xtx = np.stack([X[winter == w].T @ X[winter == w] for w in self.w])
        self.xty = np.stack([X[winter == w].T @ y[winter == w] for w in self.w])
        self.n = np.array([(winter == w).sum() for w in self.w])
    def fit(self, counts=None):
        c = np.ones(len(self.w)) if counts is None else counts
        A = np.tensordot(c, self.xtx, 1); b = np.tensordot(c, self.xty, 1)
        return np.linalg.solve(A + 1e-9 * np.eye(self.k), b)
    def boot(self, fn, B, rng):
        out = np.empty(B)
        for i in range(B):
            c = np.bincount(rng.integers(0, len(self.w), len(self.w)), minlength=len(self.w)).astype(float)
            try: out[i] = fn(self.fit(c))
            except np.linalg.LinAlgError: out[i] = np.nan
        return out
def pval(draws):
    d = draws[~np.isnan(draws)]
    return min(1.0, 2 * (min((d <= 0).sum(), (d >= 0).sum()) + 1) / (len(d) + 1))
def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.minimum(q, 1)
