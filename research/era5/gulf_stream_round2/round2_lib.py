"""Design and contrast helpers for round 2 (generalises ../gulf_stream_region/region_lib.py: any temperature column, optional categorical controls)."""
import sys, os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gulf_stream_region")); sys.path.insert(0, os.path.join(HERE, "..", "gulf_stream"))
from stats_lib import Boot, pval, bh
from region_lib import demean, contrast_high, boot_diff, pct
BINS = [(-9, 0), (3, 6), (6, 9), (9, 99)]                     # reference 0-3 C
def box_id(d): return (d.lat.round(1).astype(str) + "_" + d.lon.round(1).astype(str)).values
def cat_dummies(d, cols):
    """dummies for each categorical control; the largest category in the subset is the reference; empty categories are dropped"""
    X = []
    for c in cols:
        v = d[c].values; cats = pd.Series(v).value_counts()
        for cat in cats.index[1:]:
            if cats[cat] > 0: X.append((v == cat).astype(np.float64))
    return X
def design2(d, with_vg, dTcol="dT", extras=()):
    X = [d.lvg.values] if with_vg else []
    for lo, hi in BINS: X.append(((d[dTcol].values >= lo) & (d[dTcol].values < hi)).astype(np.float64))
    X += cat_dummies(d, extras)
    return np.column_stack(X), (1 if with_vg else 0)
def fit2(d, with_vg, dTcol="dT", extras=(), ycol="y"):
    X, k = design2(d, with_vg, dTcol, extras); ids = box_id(d)
    return Boot(demean(X, ids), demean(d[ycol].values, ids), d["winter"].values), k
def tadv_bin(x): return np.select([x <= -4, x <= -1, x < 1, x < 4], ["a_vcold", "b_cold", "c_ref", "d_warm"], "e_vwarm")
