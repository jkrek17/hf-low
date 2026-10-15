"""POST HOC (logged in PREREGISTRATION.md): the pre-registered half-split tests could not run (no class had 30 of the
rarer class in each half). All 22 seasons, season-block bootstrap (2,000), BH over the 15 class x factor tests, then the
5 pooled-across-classes tests (class dummies) as a second family. ERA5 PROXY. Also describes the 'other' class.
"""
import os, sys
import numpy as np, pandas as pd
MODE = sys.argv[1] if len(sys.argv) > 1 else "orig"      # "orig" = pre-registered classes; "s2" = S widened to 45-225 deg bearing (post hoc)
sys.argv = sys.argv[:1]
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
from scipy.optimize import minimize
pk = pd.read_csv(os.path.join(RES, "episodes.csv"), parse_dates=["time"])
if MODE == "s2":
    inb = lambda x, a, b: ((x - a) % 360) < ((b - a) % 360)
    pk["cls"] = np.where(inb(pk.bearing, 45, 225) & inb(pk.wdir, 225, 315), "S", pk.cls.replace({"S": "other"}))
    SUF = "_s2"
else:
    SUF = "_orig_classes"
FACT = ["LOW_DEPTH", "LOW_DIST", "GRAD", "GH", "COLD"]
rng = np.random.default_rng(15)
out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))


def fit(X, y, lam=0.5):   # ridge 0.5 on slopes: keeps tiny classes from separating
    X1 = np.column_stack([np.ones(len(y)), X])
    f = lambda b: np.sum(np.logaddexp(0, X1 @ b) - y * (X1 @ b)) + lam * np.sum(b[1:] ** 2)
    return minimize(f, np.zeros(X1.shape[1]), method="BFGS").x


def design(d, fac, dummies):
    cols = [fac + "_z"] + ([] if fac == "LOW_DEPTH" else ["LOW_DEPTH_z"]) + dummies
    return d[cols].values.astype(float)


def one(d, fac, dummies=()):
    return fit(design(d, fac, list(dummies)), d.hf_ep.values.astype(float))[1]


def boot(d, fac, dummies=(), n=2000):
    g = {s: x for s, x in d.groupby("season")}; ss = list(g); r = []
    for _ in range(n):
        dd = pd.concat([g[s] for s in rng.choice(ss, len(ss))])
        if dd.hf_ep.nunique() > 1: r.append(one(dd, fac, dummies))
    return np.array(r)


def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty(len(p)); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]; return np.minimum(q, 1)


def mdor(n, prev):
    return float(np.exp(2.8 / np.sqrt(n * prev * (1 - prev))))


rows = []
Z = pk[pk.cls.isin(["S", "E", "W"])].copy()
for c in ["S", "E", "W"]:
    d = Z[Z.cls == c].copy()
    for f in FACT: d[f + "_z"] = (d[f] - d[f].mean()) / d[f].std()
    for f in FACT:
        b = boot(d, f); est = one(d, f)
        rows.append(dict(family="class", cls=c, factor=f, n=len(d), hf=int(d.hf_ep.sum()), OR=np.exp(est), lo=np.exp(np.percentile(b, 2.5)),
                         hi=np.exp(np.percentile(b, 97.5)), p=2 * min((b <= 0).mean(), (b >= 0).mean()), mdor=mdor(len(d), d.hf_ep.mean())))
d = Z.copy()
for f in FACT: d[f + "_z"] = d.groupby("cls")[f].transform(lambda x: (x - x.mean()) / x.std())
d["cE"] = (d.cls == "E").astype(float); d["cW"] = (d.cls == "W").astype(float)
for f in FACT:
    b = boot(d, f, ["cE", "cW"]); est = one(d, f, ["cE", "cW"])
    rows.append(dict(family="pooled", cls="S+E+W", factor=f, n=len(d), hf=int(d.hf_ep.sum()), OR=np.exp(est), lo=np.exp(np.percentile(b, 2.5)),
                     hi=np.exp(np.percentile(b, 97.5)), p=2 * min((b <= 0).mean(), (b >= 0).mean()), mdor=mdor(len(d), d.hf_ep.mean())))
R = pd.DataFrame(rows); R["q"] = np.nan
for fam in ("class", "pooled"):
    m = R.family == fam; R.loc[m, "q"] = bh(R.loc[m, "p"].values)
R.to_csv(os.path.join(RES, "posthoc_tests" + SUF + ".csv"), index=False)
P("POST HOC, all 22 seasons, episode-level; OR per SD for hurricane-force vs storm-force proxy episodes (ridge 0.5, adjusted for LOW_DEPTH except for itself)")
P(R.round(3).to_string())
o = pk[pk.cls == "other"]
P("\n'other' class (", len(o), "episodes): bearing sector of the gust maximum from the cape and wind-from sector")
P(pd.crosstab(pd.cut(o.bearing, [0, 45, 135, 225, 315, 360], labels=["N", "E", "S", "W", "N2"], right=False).replace("N2", "N"),
              pd.cut(o.wdir, [0, 45, 135, 225, 315, 360], labels=["N", "E", "S", "W", "N2"], right=False).replace("N2", "N"), rownames=["max is"], colnames=["wind from"]).to_string())
P("\nHF share by class:", pk.groupby("cls").hf_ep.mean().round(3).to_dict(), "; season counts with >=1 episode per class:", pk.groupby("cls").season.nunique().to_dict())
open(os.path.join(RES, "posthoc" + SUF + ".txt"), "w").write("\n".join(out)); print("\n".join(out))
