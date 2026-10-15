"""Stage 2: classify events, build episodes and ingredients, run the pre-registered tests and composites.

ERA5 PROXY. Plan: PREREGISTRATION.md (committed before the pull). Fields come from stage1_pull.py ($ERA5_WORK/tip_jets/fields).
usage: analyse.py [CACHE_DIR]   -> results/*.csv, *.txt, composites.png
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "tip_jets", "fields")
CACHE = sys.argv[1] if len(sys.argv) > 1 else "/mnt/project-files/era5-cache"
RES = os.path.join(HERE, "results")
CF = (59.77, -43.9); RE = 6371.0; NM = 1.852
LAT = 90 - 0.25 * (60 + np.arange(120)); LON = -75 + 0.25 * np.arange(280)   # window, degrees east in -180..180
LA, LO = np.meshgrid(LAT, LON, indexing="ij")
lsm = np.load(os.path.join(CACHE, "lsm_nh_0p25.npy"))[60:180, 1140:1420]
OCEAN = lsm < 0.5


def dist(la0, lo0, la, lo):
    la0, lo0, la, lo = map(np.radians, (la0, lo0, la, lo))
    return RE * np.arccos(np.clip(np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0), -1, 1))


def bearing(la0, lo0, la, lo):
    la0, lo0, la, lo = map(np.radians, (la0, lo0, la, lo))
    y = np.sin(lo - lo0) * np.cos(la); x = np.cos(la0) * np.sin(la) - np.sin(la0) * np.cos(la) * np.cos(lo - lo0)
    return (np.degrees(np.arctan2(y, x)) + 360) % 360


def inb(x, a, b):          # circular interval [a,b)
    return ((x - a) % 360) < ((b - a) % 360)


def season(t):
    return t.year if t.month >= 6 else t.year - 1


box = pd.read_csv(os.path.join(RES, "stage0_box_gust.csv"), parse_dates=["time"])
pt = pd.read_csv(os.path.join(RES, "pull_times.csv"), parse_dates=["time"])
ev = box[box.box_gust_kt >= 54].copy().reset_index(drop=True)
D_CF = dist(CF[0], CF[1], LA, LO)
GHM = (~OCEAN) & (LA >= 64) & (LA <= 78) & (LO >= -55) & (LO <= -30)
UPS = OCEAN & (LA >= 56) & (LA <= 63) & (LO >= -60) & (LO <= -46)
LOWM = OCEAN & (LA >= 50) & (LA <= 70) & (LO >= -60) & (LO <= -10)
KT = 1 / 0.514444


def load(t):
    z = np.load(os.path.join(WORK, t.strftime("%Y%m%d%H") + ".npz"))
    return {k: z[k] for k in z.files}


rows = []
for i, r in ev.iterrows():
    f = load(r.time)
    ri = int(round((90 - r.max_lat) / 0.25)) - 60; ci = int(round((r.max_lon + 75) / 0.25))
    u = f["u10"][ri-1:ri+2, ci-1:ci+2].mean(); v = f["v10"][ri-1:ri+2, ci-1:ci+2].mean()
    wdir = (np.degrees(np.arctan2(-u, -v)) + 360) % 360       # direction wind blows FROM
    brg = bearing(CF[0], CF[1], r.max_lat, r.max_lon)
    cls = "other"
    if inb(brg, 135, 225) and inb(wdir, 225, 315): cls = "S"
    elif inb(brg, 45, 135) and inb(wdir, 315, 45): cls = "E"
    elif inb(brg, 225, 315) and inb(wdir, 45, 135): cls = "W"
    msl = f["msl"] / 100.0
    m = np.where(LOWM, msl, np.inf); k = np.unravel_index(np.argmin(m), m.shape)
    gy, gx = np.gradient(msl, 0.25 * 111.0, 0.25 * 111.0)                            # hPa per km, before the cos(lat) correction
    gm = np.hypot(gy, gx / np.cos(np.radians(LA))) * 100                             # hPa per 100 km
    near = OCEAN & (dist(r.max_lat, r.max_lon, LA, LO) <= 200)
    rows.append(dict(time=r.time, box_gust_kt=r.box_gust_kt, max_lat=r.max_lat, max_lon=r.max_lon, wdir=wdir, wspd_kt=np.hypot(u, v) * KT,
                     bearing=brg, cls=cls, low_mslp=m[k], low_lat=LA[k], low_lon=LO[k], LOW_DEPTH=-m[k],
                     LOW_DIST=dist(CF[0], CF[1], LA[k], LO[k]), GRAD=gm[near].mean(), GH=msl[GHM].mean(),
                     COLD=np.nanmean((f["sst"] - f["t2m"])[UPS])))
E = pd.DataFrame(rows)
E["season"] = E.time.apply(season)
E["hf"] = E.box_gust_kt >= 71.7
# episodes: link events <= 36 h apart
E = E.sort_values("time").reset_index(drop=True)
E["ep"] = (E.time.diff() > pd.Timedelta(hours=36)).cumsum()
pk = E.loc[E.groupby("ep").box_gust_kt.idxmax()].copy()
pk["hf_ep"] = E.groupby("ep").hf.max().reindex(pk.ep).values
pk["n_times"] = E.groupby("ep").size().reindex(pk.ep).values
E.to_csv(os.path.join(RES, "event_times.csv"), index=False)
pk.to_csv(os.path.join(RES, "episodes.csv"), index=False)

out = []
P = lambda *a: out.append(" ".join(str(x) for x in a))
P("ERA5 PROXY. Event times", len(E), "; episodes", len(pk), "; seasons", pk.season.nunique())
tab = pk.groupby("cls").agg(episodes=("ep", "size"), hf=("hf_ep", "sum"), med_gust=("box_gust_kt", "median"), med_wspd=("wspd_kt", "median"))
tab["hf_share"] = tab.hf / tab.episodes
P(tab.round(3).to_string())
P("by half (even seasons = fit):"); pk["half"] = np.where(pk.season % 2 == 0, "even", "odd")
P(pk.groupby(["cls", "half"]).agg(n=("ep", "size"), hf=("hf_ep", "sum")).to_string())

# ---- tests
FACT = ["LOW_DEPTH", "LOW_DIST", "GRAD", "GH", "COLD"]
rng = np.random.default_rng(15)


def fit(X, y):
    X1 = np.column_stack([np.ones(len(y)), X])
    def nll(b):
        z = X1 @ b; return np.sum(np.logaddexp(0, z) - y * z) + 1e-3 * np.sum(b[1:] ** 2)
    return minimize(nll, np.zeros(X1.shape[1]), method="BFGS").x


def est(d, fac, cols):
    X = d[cols].values; return fit(X, d.hf_ep.values.astype(float))[1 + cols.index(fac)]


def boot(d, fac, cols, n=1000):
    ss = d.season.unique(); g = {s: d[d.season == s] for s in ss}; r = []
    for _ in range(n):
        dd = pd.concat([g[s] for s in rng.choice(ss, len(ss))])
        if dd.hf_ep.nunique() < 2: continue
        r.append(est(dd, fac, cols))
    return np.array(r)


def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.empty(len(p)); m = len(p)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def mdor(n, prev):   # approximate, 80% power, two-sided 0.05, episodes treated as independent
    return float(np.exp((1.96 + 0.84) / np.sqrt(n * prev * (1 - prev)))) if 0 < prev < 1 else np.nan


Z = pk[pk.cls.isin(["S", "E", "W"])].copy()
for f_ in FACT:
    Z[f_ + "_z"] = Z.groupby("cls")[f_].transform(lambda x: (x - x.mean()) / x.std())
Z["half"] = np.where(Z.season % 2 == 0, "even", "odd")
res = []
for swap, (fh, th) in enumerate([("even", "odd"), ("odd", "even")]):
    for c in ["S", "E", "W"]:
        d = Z[Z.cls == c]; a = d[d.half == fh]; b = d[d.half == th]
        rare = min(a.hf_ep.sum(), (~a.hf_ep.astype(bool)).sum())
        for f_ in FACT:
            cols = [f_ + "_z"] + ([] if f_ == "LOW_DEPTH" else ["LOW_DEPTH_z"])
            row = dict(swap=swap, cls=c, factor=f_, fit_half=fh, n_fit=len(a), hf_fit=int(a.hf_ep.sum()), n_test=len(b), hf_test=int(b.hf_ep.sum()),
                       tested=bool(rare >= 30 and min(b.hf_ep.sum(), (~b.hf_ep.astype(bool)).sum()) >= 30), mdor=mdor(len(a), a.hf_ep.mean()))
            if a.hf_ep.nunique() > 1 and b.hf_ep.nunique() > 1:
                bf = boot(a, f_ + "_z", cols); bt = boot(b, f_ + "_z", cols)
                row.update(or_fit=np.exp(est(a, f_ + "_z", cols)), lo=np.exp(np.percentile(bf, 2.5)), hi=np.exp(np.percentile(bf, 97.5)),
                           p_fit=2 * min((bf <= 0).mean(), (bf >= 0).mean()), or_test=np.exp(est(b, f_ + "_z", cols)),
                           p_test_one=(bt <= 0).mean() if est(a, f_ + "_z", cols) > 0 else (bt >= 0).mean())
            res.append(row)
R = pd.DataFrame(res)
R["q_fit"] = np.nan
for s in (0, 1):
    m = (R.swap == s) & R.p_fit.notna()
    R.loc[m, "q_fit"] = bh(R.loc[m, "p_fit"].values)
R["supported"] = R.tested & (R.q_fit < 0.05) & (R.p_test_one < 0.10) & ((R.or_fit > 1) == (R.or_test > 1))
R.to_csv(os.path.join(RES, "tests.csv"), index=False)
P("\nDO NOT QUOTE: the pre-registered tests did not meet the 30-per-half rule (tested = False everywhere); see README.")
P("\nTESTS (OR per SD of factor for HF vs storm force, adjusted for low depth except LOW_DEPTH; BH over 15 within each swap direction)")
P(R.round(3).to_string())
# pooled (all seasons) description
P("\nPOOLED (all 22 seasons, descriptive): class x factor median at HF vs S")
for c in ["S", "E", "W"]:
    d = Z[Z.cls == c]
    P(c, d.groupby("hf_ep")[FACT].median().round(2).to_string())
open(os.path.join(RES, "results.txt"), "w").write("\n".join(out))
print("\n".join(out))
