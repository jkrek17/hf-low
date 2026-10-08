"""POST HOC diagnostics of why the matched sensitivity disagrees with the adjusted logistic on GH.

Written after the primary results and the matched sensitivity were read; none of this was pre-registered and
none of it enters the F1 family. ERA5 proxy. Needs analyse.py to have run. Writes results/posthoc.txt.
Each row refits the GH effect (logit per SD, ingredients GRAD, STAB, NAO, motion kept) under a different way of
adjusting for the low, so the reader can see how much the answer depends on that choice.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
RES = os.path.join(HERE, "results")
sys.path.insert(0, HERE)
import analyse as an
out = []


def say(s):
    print(s); out.append(s)


def gh_effect(z, cols, nboot=500, seed=1):
    y = z.case.values.astype(float)
    X = an.mat(z, cols)
    j = cols.index("GH") + 1
    b = an.logit(X, y)
    sid = an.seasons_idx(z)
    rng = np.random.default_rng(seed)
    bb = np.array([an.logit(X[ix], y[ix])[j] for ix in (an.boot_idx(sid, rng) for _ in range(nboot))])
    return b[j], np.percentile(bb, 2.5), np.percentile(bb, 97.5)


z = pd.read_pickle(f"{WORK}/zall.pkl").reset_index(drop=True)
say("POST HOC diagnostics (not pre-registered). ERA5 proxy, outcome Y_T (G_T >= 71.7 kt), all 22 seasons unless stated")
c = z[["GH", "GRAD", "STAB", "NAO", "depth", "x", "y", "Dcf"]].corr()
say("weighted-sample correlations with GH are below; unweighted over the 1,076 sample rows: " +
    ", ".join(f"{k} {c.loc['GH', k]:+.2f}" for k in ["GRAD", "STAB", "NAO", "depth", "x", "y", "Dcf"]))
cc = z[z.case == 1][["GH", "depth", "x", "y", "Dcf"]].corr()
say("within cases only, GH with: " + ", ".join(f"{k} {cc.loc['GH', k]:+.2f}" for k in ["depth", "x", "y", "Dcf"]))
z["d2"] = z.depth ** 2
z["dx"], z["dy"] = z.depth * z.x, z.depth * z.y
z["d3"] = z.depth ** 3
o = np.argsort(z.depth.values)
cw = np.cumsum(z.w.values[o]) / z.w.sum()
q1, q2 = z.depth.values[o][np.searchsorted(cw, 1 / 3)], z.depth.values[o][np.searchsorted(cw, 2 / 3)]
z["t1"] = (z.depth <= q1).astype(float)
z["t2"] = ((z.depth > q1) & (z.depth <= q2)).astype(float)
base_ing = ["GH", "GRAD", "STAB", "NAO", "u", "v"]
specs = {
    "pre-registered baseline (linear depth, quadratic position)": an.BASE + base_ing,
    "+ depth squared and cubed": an.BASE + ["d2", "d3"] + base_ing,
    "+ depth squared and depth x position": an.BASE + ["d2", "dx", "dy"] + base_ing,
    "depth in weighted terciles (+ linear within)": an.BASE + ["t1", "t2"] + base_ing,
}
for lbl, cols in specs.items():
    for nm, zz in (("all seasons", z), ("even-start seasons", z[z.half == "discovery"].reset_index(drop=True)),
                   ("odd-start seasons", z[z.half == "heldout"].reset_index(drop=True))):
        b, lo, hi = gh_effect(zz, cols)
        say(f"  {lbl:55s} {nm:18s} GH OR per SD {np.exp(b):.2f} [{np.exp(lo):.2f}, {np.exp(hi):.2f}]")
deep = z[z.depth <= q1].reset_index(drop=True)
say(f"\nrestricted to the deepest weighted tercile of reference-low MSLP (<= {q1 * an.wstats(pd.read_csv(f'{RES}/sample_times.csv').assign(season=0), ['ref_msl'])['ref_msl'][1] + an.wstats(pd.read_csv(f'{RES}/sample_times.csv'), ['ref_msl'])['ref_msl'][0]:.0f} hPa): "
    f"{int(deep.case.sum())} cases in {len(deep)} rows")
for lbl, cols in (("baseline", an.BASE + base_ing), ("baseline + depth squared", an.BASE + ["d2"] + base_ing)):
    b, lo, hi = gh_effect(deep, cols)
    say(f"  {lbl:30s} GH OR per SD {np.exp(b):.2f} [{np.exp(lo):.2f}, {np.exp(hi):.2f}]")
for k in ("STAB", "GRAD", "NAO"):
    j = (an.BASE + base_ing).index(k) + 1
    y = deep.case.values.astype(float)
    X = an.mat(deep, an.BASE + base_ing)
    b = an.logit(X, y)[j]
    sid = an.seasons_idx(deep)
    rng = np.random.default_rng(2)
    bb = np.array([an.logit(X[ix], y[ix])[j] for ix in (an.boot_idx(sid, rng) for _ in range(500))])
    say(f"  deepest tercile, {k} OR per SD {np.exp(b):.2f} [{np.exp(np.percentile(bb, 2.5)):.2f}, {np.exp(np.percentile(bb, 97.5)):.2f}]")
open(f"{RES}/posthoc.txt", "w").write("\n".join(out) + "\n")
