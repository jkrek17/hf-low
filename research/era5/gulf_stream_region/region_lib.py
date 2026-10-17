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
def box_id(d): return (d.lat.round(1).astype(str) + "_" + d.lon.round(1).astype(str)).values
def design_raw(d, with_vg):
    """columns: [log Vg if with_vg] + four dT bins (reference 0-3 C); box fixed effects are absorbed by within-box demeaning in fit_boot."""
    X = [d.lvg.values] if with_vg else []
    for lo, hi in BINS: X.append(((d.dT.values >= lo) & (d.dT.values < hi)).astype(np.float64))
    return np.column_stack(X), (1 if with_vg else 0)
def demean(M, ids):
    out = M.astype(np.float64).copy(); codes = pd.factorize(ids)[0]; n = codes.max() + 1; cnt = np.bincount(codes, minlength=n)
    if out.ndim == 1: out -= (np.bincount(codes, out, minlength=n) / cnt)[codes]
    else:
        for c in range(out.shape[1]): out[:, c] -= (np.bincount(codes, out[:, c], minlength=n) / cnt)[codes]
    return out
def contrast_high(k): return lambda b: 0.5 * (b[k + 2] + b[k + 3])          # mean of 6-9 and >=9 bins
def fit_boot(d, with_vg, winter_col="winter"):
    X, k = design_raw(d, with_vg); ids = box_id(d)
    return Boot(demean(X, ids), demean(d.y.values, ids), d[winter_col].values), k
def boot_diff(B1, f1, B2, f2, nb, rng):
    out = np.empty(nb)
    for i in range(nb):
        c1 = np.bincount(rng.integers(0, len(B1.w), len(B1.w)), minlength=len(B1.w)).astype(float)
        c2 = np.bincount(rng.integers(0, len(B2.w), len(B2.w)), minlength=len(B2.w)).astype(float)
        out[i] = f1(B1.fit(c1)) - f2(B2.fit(c2))
    return out
def pct(x): return 100 * (np.exp(x) - 1)
