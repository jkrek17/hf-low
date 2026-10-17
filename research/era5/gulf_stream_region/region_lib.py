"""Designs and contrasts for the region-wide tests (see PREREGISTRATION.md). Reuses the winter-block bootstrap of ../gulf_stream/stats_lib.py."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gulf_stream"))
from stats_lib import Boot, pval, bh
BINS = [(-9, 0), (3, 6), (6, 9), (9, 99)]          # reference 0-3 C
def prep_A(A):
    d = A.dropna(subset=["v_ascat", "vg", "dT"]); d = d[(d.vg >= 15) & (d.v_ascat >= 5)].copy()
    d["y"] = np.log(d.v_ascat); d["lvg"] = np.log(d.vg); return d
def prep_B(B):
    d = B.dropna(subset=["v_ascat", "v_gfs", "dT"]); d = d[(d.v_gfs >= 5) & (d.v_ascat >= 5)].copy()
    d["y"] = np.log(d.v_ascat / d.v_gfs); return d
def design(d, with_vg):
    box = (d.lat.astype(str) + "_" + d.lon.astype(str)).values; ub = np.unique(box)
    X = [(box == b).astype(np.float32) for b in ub]
    if with_vg: X.append(d.lvg.values.astype(np.float32))
    for lo, hi in BINS: X.append(((d.dT.values >= lo) & (d.dT.values < hi)).astype(np.float32))
    X = np.column_stack(X); k = len(ub) + (1 if with_vg else 0)
    keep = X[:, :len(ub)].sum(0) > 0
    return X, k
def contrast_high(k): return lambda b: 0.5 * (b[k + 2] + b[k + 3])          # mean of 6-9 and >=9 bins
def fit_boot(d, with_vg, winter_col="winter"):
    X, k = design(d, with_vg); return Boot(X.astype(np.float64), d.y.values.astype(np.float64), d[winter_col].values), k
def boot_diff(B1, f1, B2, f2, nb, rng):
    out = np.empty(nb)
    for i in range(nb):
        c1 = np.bincount(rng.integers(0, len(B1.w), len(B1.w)), minlength=len(B1.w)).astype(float)
        c2 = np.bincount(rng.integers(0, len(B2.w), len(B2.w)), minlength=len(B2.w)).astype(float)
        out[i] = f1(B1.fit(c1)) - f2(B2.fit(c2))
    return out
def pct(x): return 100 * (np.exp(x) - 1)
