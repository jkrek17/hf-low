"""Descriptive risk tables for the answer in plain words (added after the primary run; not tests).

Weighted share of candidate times with G_T >= 71.7 kt (ERA5 proxy), pooled over all 22 seasons, by terciles of
pairs of ingredients, to show whether two ingredients together do more or less than their separate effects on the
risk scale. Reads zall.pkl from analyse.py. Writes results/tables.txt.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
RES = os.path.join(HERE, "results")
out = []


def say(s):
    print(s); out.append(s)


def wq(v, w, qs):
    o = np.argsort(v); cw = np.cumsum(w[o]) / w.sum()
    return [v[o][np.searchsorted(cw, q)] for q in qs]


z = pd.read_pickle(f"{WORK}/zall.pkl")
say("Weighted share of low-in-region times with G_T >= 71.7 kt (ERA5 proxy), all 22 seasons, by tercile of two ingredients")
say("(z is in SD of the fit-half weighted distribution; 'low' = bottom tercile, 'high' = top tercile of the weighted sample)")
for a, b in (("GH", "depth"), ("GH", "STAB"), ("GH", "GRAD"), ("depth", "STAB")):
    ta = wq(z[a].values, z.w.values, [1 / 3, 2 / 3]); tb = wq(z[b].values, z.w.values, [1 / 3, 2 / 3])
    ca = np.digitize(z[a], ta); cb = np.digitize(z[b], tb)
    nm = {"depth": "reference-low MSLP (high = shallow)"}
    say(f"\nrows: {nm.get(a, a)} tercile; columns: {nm.get(b, b)} tercile (cells: share, n cases/n rows)")
    lab = ["low", "mid", "high"]
    say("            " + "  ".join(f"{l:>20s}" for l in lab))
    for i in range(3):
        cells = []
        for j in range(3):
            m = (ca == i) & (cb == j)
            p = np.average(z.case[m], weights=z.w[m]) if m.sum() else np.nan
            cells.append(f"{p:6.3f} ({int(z.case[m].sum()):3d}/{int(m.sum()):4d})")
        say(f"  {lab[i]:8s}  " + "  ".join(f"{c:>20s}" for c in cells))
open(f"{RES}/tables.txt", "w").write("\n".join(out) + "\n")
