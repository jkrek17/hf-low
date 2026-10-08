"""Station check: do the cases, and the ingredient-favoured times, show stronger observed wind?

ERA5 proxy cases against NOAA ISD hourly reports at Prins Christian Sund (04390099999,
60.0N 43.1W) and Tasiilaq (04360099999, 65.6N 37.6W), held-out seasons only
(PREREGISTRATION.md, family F5). Needs analyse.py to have run (it writes zall.pkl and
bfit.npy in $ERA5_WORK/greenland_jets). ISD files are read from $ERA5_WORK/isd/<station>_<year>.csv
(https://noaa-global-hourly-pds.s3.amazonaws.com/<year>/<station>.csv). Writes results/station.txt.

Tests, 12 in one BH family: for each metric (mean wind in kt, share of reports >= 34 kt)
  A  cases against controls: all cases at both stations, tip cases at Prins Christian Sund,
     barrier cases at Tasiilaq (4 contrasts)
  B  top quintile of the ingredient score against the rest, at both stations (2 contrasts)
The plan's "8 tests" under-counted the tip/barrier variants; the family is 12.
"""
import os, sys, glob
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
ISD = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "isd")
RES = os.path.join(HERE, "results")
sys.path.insert(0, HERE)
import analyse as an
KT = 1 / 0.514444
BAD = set("2367")
NBOOT = int(next((a for a in sys.argv[1:] if a.isdigit()), 2000))
STN = {"PCS": "04390099999", "TAS": "04360099999"}
out = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)


def parse(path):
    d = pd.read_csv(path, usecols=lambda c: c in ("DATE", "WND"), dtype=str)
    w = d.WND.str.split(",", expand=True)
    spd = pd.to_numeric(w[3], errors="coerce")
    ok = (spd != 9999) & ~w[4].isin(BAD)
    d["wind"] = np.where(ok, spd / 10 * KT, np.nan)
    d["t"] = pd.to_datetime(d.DATE)
    return d[["t", "wind"]].dropna()


def series(code):
    fs = sorted(glob.glob(f"{ISD}/{code}_*.csv"))
    s = pd.concat([parse(p) for p in fs if os.path.getsize(p) > 1000]).sort_values("t")
    return s.reset_index(drop=True)


def at(s, times):
    """Wind of the report nearest each time, within +-1 h (NaN when none)."""
    st = s.t.values.astype("datetime64[ns]")
    tv = pd.DatetimeIndex(times).values.astype("datetime64[ns]")
    i = np.searchsorted(st, tv)
    c0, c1 = np.clip(i - 1, 0, len(st) - 1), np.clip(i, 0, len(st) - 1)
    d0, d1 = np.abs(st[c0] - tv), np.abs(st[c1] - tv)
    best = np.where(d0 <= d1, c0, c1)
    dd = np.minimum(d0, d1)
    w = s.wind.values[best].astype(float)
    return np.where(dd <= np.timedelta64(1, "h"), w, np.nan)


def wmean(v, w):
    ok = np.isfinite(v)
    return np.average(v[ok], weights=w[ok]) if ok.sum() else np.nan


def contrast(z, v, a, b):
    """Weighted mean wind and share >= 34 kt, group a minus group b, from rows with a report."""
    r = []
    for f in (lambda x: x, lambda x: np.where(np.isfinite(x), (x >= 34).astype(float), np.nan)):
        x = f(v)
        r.append(wmean(x[a], z.w.values[a]) - wmean(x[b], z.w.values[b]))
    return r


def main():
    z = pd.read_pickle(f"{WORK}/zall.pkl")
    bfit = np.load(f"{WORK}/bfit.npy")
    cols = an.BASE + an.ING + an.MOT
    z = z[z.half == "heldout"].reset_index(drop=True)
    ing = an.ING + an.MOT
    z["score"] = sum(bfit[cols.index(c) + 1] * z[c] for c in ing)
    t = pd.to_datetime(z.time.astype(str), format="%Y%m%d%H")
    o = np.argsort(z.score.values)
    cw = np.cumsum(z.w.values[o]) / z.w.sum()
    thr = z.score.values[o][np.searchsorted(cw, 0.8)]
    z["top"] = z.score >= thr
    say("Station check (ERA5 proxy cases vs ISD reports), held-out seasons (odd start years), cases and controls as sampled")
    say(f"held-out rows {len(z)}, cases {int(z.case.sum())} (tip {int(z.tip.sum())}, barrier {int(z.barrier.sum())}); "
        f"top quintile of the ingredient score (weighted) = score >= {thr:.2f}, {int(z.top.sum())} rows")
    wind = {}
    for k, code in STN.items():
        s = series(code)
        wind[k] = at(s, t)
        z[k] = wind[k]
        cov = np.isfinite(wind[k])
        say(f"{k}: reports {len(s)}, years {s.t.dt.year.min()}-{s.t.dt.year.max()}; matched within +-1 h for "
            f"{cov.sum()} of {len(z)} rows (cases {int((cov & (z.case == 1)).sum())} of {int(z.case.sum())})")
    sidx = an.seasons_idx(z)
    rng = np.random.default_rng(an.SEED + 40)
    spec = [("A all cases", "PCS", z.case == 1, z.case == 0), ("A all cases", "TAS", z.case == 1, z.case == 0),
            ("A tip cases", "PCS", z.tip, z.case == 0), ("A barrier cases", "TAS", z.barrier, z.case == 0),
            ("B top quintile", "PCS", z.top, ~z.top), ("B top quintile", "TAS", z.top, ~z.top)]
    obs = []
    for name, st, a, b in spec:
        obs.append(contrast(z, wind[st], a.values, b.values))
    bs = np.zeros((NBOOT, len(spec), 2))
    zz = z.reset_index(drop=True)
    for i in range(NBOOT):
        ix = an.boot_idx(sidx, rng)
        zb = zz.iloc[ix].reset_index(drop=True)
        for j, (name, st, a, b) in enumerate(spec):
            bs[i, j] = contrast(zb, wind[st][ix], a.values[ix], b.values[ix])
    ps, lines = [], []
    for j, (name, st, a, b) in enumerate(spec):
        for m, lab in enumerate(("mean wind (kt)", "share >= 34 kt")):
            e = obs[j][m]
            bb = bs[:, j, m]
            bb = bb[np.isfinite(bb)]
            p = min(1.0, 2 * min((bb >= 0).mean(), (bb <= 0).mean()) + 1 / max(len(bb), 1))
            if np.nanstd(bb) == 0:
                p = 1.0                                   # degenerate: no report at or above 34 kt in either group
            ps.append(p)
            lines.append((name, st, lab, e, np.percentile(bb, 2.5), np.percentile(bb, 97.5), p,
                          int((a.values & np.isfinite(wind[st])).sum()), int((b.values & np.isfinite(wind[st])).sum())))
    q = an.bh(np.array(ps))
    say("\nF5 (12 tests, BH): group minus comparison; season-block bootstrap 95% CI (held-out seasons)")
    for (name, st, lab, e, lo, hi, p, na, nb), qq in zip(lines, q):
        fmt = (lambda v: f"{v:+.1f}") if lab.startswith("mean") else (lambda v: f"{v:+.3f}")
        say(f"  {name:16s} {st}  {lab:15s} {fmt(e)} [{fmt(lo)}, {fmt(hi)}]  p {p:.3f} q {qq:.3f}  (n with report: {na} vs {nb})")
    # POST HOC (added after the first station results; not in the F5 family): weighted regression of the
    # station wind on the baseline and ingredient terms, season-block bootstrap on the GH coefficient
    say("\nPOST HOC, not pre-registered: station wind (kt) regressed on BASE + ingredients (weighted, held-out rows with a report)")
    cols = an.BASE + an.ING + an.MOT
    for k in STN:
        ok = np.isfinite(wind[k])
        zk = z[ok].reset_index(drop=True)
        yk = wind[k][ok]
        X = an.mat(zk, cols)
        w = zk.w.values
        beta = an.wls(X, yk, w)
        sid = an.seasons_idx(zk)
        r2 = np.random.default_rng(an.SEED + 50)
        bb = []
        for _ in range(NBOOT // 2):
            ix = an.boot_idx(sid, r2)
            bb.append(an.wls(X[ix], yk[ix], w[ix]))
        bb = np.array(bb)
        txt = []
        for c in an.ING:
            j = cols.index(c) + 1
            txt.append(f"{c} {beta[j]:+.2f} kt/SD [{np.percentile(bb[:, j], 2.5):+.2f}, {np.percentile(bb[:, j], 97.5):+.2f}]")
        say(f"  {k} (n {len(zk)}): " + "; ".join(txt))
    say("\ngroup means (weighted, kt):")
    for k in STN:
        for nm, g in (("cases", z.case == 1), ("tip", z.tip), ("barrier", z.barrier), ("controls", z.case == 0),
                      ("top quintile", z.top), ("rest", ~z.top)):
            v = wind[k]
            gg = g.values
            say(f"  {k} {nm:13s} mean {wmean(v[gg], z.w.values[gg]):.1f}, share >= 34 kt "
                f"{wmean(np.where(np.isfinite(v), (v >= 34) * 1.0, np.nan)[gg], z.w.values[gg]):.3f}, n {int((gg & np.isfinite(v)).sum())}")
    open(f"{RES}/station.txt", "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
