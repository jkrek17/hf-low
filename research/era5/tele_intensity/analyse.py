"""Do NAO and PNA change how intense HF lows get, once the index's own pressure
signature and the shift in storm position are taken out? (ERA5 proxy, pipeline A)

Sample: the 4,157 events of the pipeline A catalog (research/era5/hf_history),
HF-equivalent by the 800 km gust index >= 71.7 kt. Everything here is a proxy,
and is intensity GIVEN that a storm reached the HF-equivalent threshold.

Outcomes, per event
  minp       lowest in-domain central pressure (hPa), = catalog minp
  anom_clim  minp - calendar-month climatological MSLP at that position
  depth_ring minp - mean MSLP over the 900-1100 km annulus at the same time
  bg         ring - clim: the large-scale background anomaly around the storm
             (anom_clim = depth_ring + bg exactly, and so do the fitted slopes)
  ndr_max    largest 24 h normalised deepening rate on the track (Bergerons,
             sign: positive = deepening)
  g800max    largest in-domain 800 km gust index (kt)        seasons 2004-05 on
  hf_hours   6 h x in-domain fixes with gust index >= 71.7 kt  seasons 2004-05 on
Pressure outcomes use all 47 seasons (depth did not drift before 2001,
STATUS.md); gust outcomes start 2004-05 (the gust drift gate and RECORD_START).

Predictors: CPC daily NAO and PNA, averaged over days -10..-4 before the
track's first fix (PRIMARY: set before the storm exists, so the storm cannot
feed it), and over the 5 days ending on the day of the deepest fix (CONTRAST
only). Both indices enter together (distinct regions), z-scored on the sample.

Model, per basin and outcome: OLS of the outcome on NAO + PNA + calendar-month
fixed effects (month of the deepest fix) + a linear season term. "position"
variant adds lat, lat^2, lon, lon^2 of the deepest fix, i.e. the effect at a
fixed place: the gap between the two is the part carried by the track shift.

Inference: season-block bootstrap SE (whole seasons resampled), season-block
permutation p (each event takes the index for the same day-of-season in another
season under a random permutation of seasons; the two-sided p counts |b*| >= |b|),
leave-one-season-out range of the slope. Benjamini-Hochberg q over the primary
table (lagged index, both basins, both indices, all outcomes, no position).

usage: analyse.py CATALOG_CSV TRACKS_CSV BACKGROUND_CSV CPC_DIR OUT_DIR [NBOOT]
CPC_DIR holds norm.daily.{nao,pna}.index.b500101.current.ascii (CPC, daily).
"""
import os, sys, json
import numpy as np, pandas as pd

HF = 71.7
SIN60 = np.sin(np.radians(60))
GUST_FROM = 2004
OUTCOMES = ["minp", "anom_clim", "depth_ring", "bg", "ndr_max", "g800max", "hf_hours"]
GUST = {"g800max", "hf_hours"}
UNITS = {"minp": "hPa", "anom_clim": "hPa", "depth_ring": "hPa", "bg": "hPa",
         "ndr_max": "Bergeron", "g800max": "kt", "hf_hours": "h"}


def read_cpc(path):
    x = pd.read_fwf(path, widths=[4, 3, 3, 999], header=None, names=["y", "m", "d", "v"])
    x["v"] = pd.to_numeric(x.v, errors="coerce")
    x.loc[x.v < -90, "v"] = np.nan
    s = pd.Series(x.v.values, index=pd.to_datetime(dict(year=x.y, month=x.m, day=x.d)))
    return s.asfreq("D")


def season_of(d):
    return np.where(d.month >= 6, d.year, d.year - 1)


def by_season_day(series, s0, s1):
    """Array [season, day-of-season 0..365] of a daily series; seasons start 1 June."""
    A = np.full((s1 - s0 + 1, 366), np.nan)
    for s in range(s0, s1 + 1):
        v = series.reindex(pd.date_range(f"{s}-06-01", f"{s + 1}-05-31")).values
        A[s - s0, :len(v)] = v
    return A


def track_metrics(T):
    T = T.sort_values(["track", "time"])
    out = []
    for tid, g in T.groupby("track", sort=False):
        msl, lat = g.msl.values, g.lat.values
        if len(g) > 4:
            ml = np.radians(0.5 * (lat[:-4] + lat[4:]))
            ndr = (msl[:-4] - msl[4:]) / 24 * SIN60 / np.abs(np.sin(ml))
            ndr_max = float(np.nanmax(ndr))
        else:
            ndr_max = np.nan
        d = g[g.basin.notna()]
        out.append(dict(track=tid, ndr_max=ndr_max, g800max=d.g800.max(), hf_hours=6.0 * (d.g800 >= HF).sum()))
    return pd.DataFrame(out)


def design(d, xs, position):
    cols = [np.ones(len(d))]
    cols += [d[x].values for x in xs]
    for m in sorted(d.month.unique())[1:]:
        cols.append((d.month.values == m).astype(float))
    cols.append(d.season.values - d.season.mean())
    if position:
        la = d.lat_min.values - d.lat_min.mean(); lo = d.lon_c.values - d.lon_c.mean()
        cols += [la, la ** 2, lo, lo ** 2]
    return np.column_stack(cols)


def ols(d, y, xs, position):
    X = design(d, xs, position)
    b = np.linalg.lstsq(X, d[y].values, rcond=None)[0]
    return b[1:1 + len(xs)]


def zs(d, xs, ref):
    d = d.copy()
    for x in xs:
        d[x] = (d[x] - ref[x][0]) / ref[x][1]
    return d


def fit_all(d, y, xs, position, I, s0, rng, nboot):
    d = d.dropna(subset=[y] + xs)
    ref = {x: (d[x].mean(), d[x].std()) for x in xs}
    dz = zs(d, xs, ref)
    b = ols(dz, y, xs, position)
    seasons = np.sort(d.season.unique())
    groups = {s: dz[dz.season == s] for s in seasons}
    bb = np.array([ols(pd.concat([groups[s] for s in rng.choice(seasons, len(seasons))]), y, xs, position)
                   for _ in range(nboot)])
    se = bb.std(axis=0, ddof=1)
    loso = np.array([ols(dz[dz.season != s], y, xs, position) for s in seasons])
    # season-block permutation: index from the same day-of-season in a permuted season
    ge_s = d.gen_season.values - s0
    dos = d.gen_dos.values
    allseas = np.arange(I[xs[0]].shape[0])
    cnt = np.zeros(len(xs))
    nperm = 0
    for _ in range(nboot):
        perm = rng.permutation(allseas)
        dp = d.copy()
        for x in xs:
            dp[x] = I[x][perm[ge_s], dos]
        dp = dp.dropna(subset=xs)
        if len(dp) < 0.9 * len(d):
            continue
        dp = zs(dp, xs, {x: (dp[x].mean(), dp[x].std()) for x in xs})
        cnt += np.abs(ols(dp, y, xs, position)) >= np.abs(b)
        nperm += 1
    p = (cnt + 1) / (nperm + 1)
    return dict(n=len(d), seasons=len(seasons), b=b, se=se, p=p, loso_lo=loso.min(0), loso_hi=loso.max(0), nperm=nperm)


def bh(ps):
    ps = np.asarray(ps); n = len(ps); o = np.argsort(ps)
    q = np.empty(n); prev = 1.0
    for rank, i in reversed(list(enumerate(o, 1))):
        prev = min(prev, ps[i] * n / rank); q[i] = prev
    return q


def main():
    cat_csv, tracks_csv, bg_csv, cpc, out = sys.argv[1:6]
    nboot = int(sys.argv[6]) if len(sys.argv) > 6 else 1000
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(20261008)
    C = pd.read_csv(cat_csv)
    E = C[C.role == "event"][["track", "basin", "season", "start", "minp"]].copy()
    T = pd.read_csv(tracks_csv, dtype={"time": str})
    B = pd.read_csv(bg_csv, dtype={"tmin": str})
    E = E.merge(track_metrics(T[T.track.isin(E.track)]), on="track").merge(
        B.drop(columns="minp"), on="track", how="left")
    assert np.allclose(E.minp, C.set_index("track").loc[E.track, "minp"].values)
    E["anom_clim"] = E.minp - E.clim
    E["depth_ring"] = E.minp - E.ring
    E["bg"] = E.ring - E.clim
    tmin = pd.to_datetime(E.tmin, format="%Y%m%d%H")
    E["month"] = tmin.dt.month
    E["lon_c"] = np.where(E.basin == "atl", (E.lon_min + 180) % 360 - 180, E.lon_min)
    gen = pd.to_datetime(E.start.astype(str), format="%Y%m%d%H").dt.normalize()
    E["gen_season"] = season_of(gen.dt)
    E["gen_dos"] = (gen - pd.to_datetime(E.gen_season.astype(str) + "-06-01")).dt.days.values

    s0, s1 = 1977, 2026
    I = {}
    for k in ("NAO", "PNA"):
        s = read_cpc(os.path.join(cpc, f"norm.daily.{k.lower()}.index.b500101.current.ascii"))
        lag = s.rolling(7, min_periods=7).mean().shift(4)          # mean of days -10..-4
        I[k] = by_season_day(lag, s0, s1)
        E[k] = I[k][E.gen_season - s0, E.gen_dos]
        now = s.rolling(5, min_periods=5).mean()                    # days -4..0 of the deepest fix
        E[k + "_now"] = now.reindex(tmin.dt.normalize()).values

    res = []
    for basin in ("atl", "pac"):
        for era in ("all", "2004+"):
            for y in OUTCOMES:
                if y in GUST and era == "all":
                    continue
                lo = GUST_FROM if era == "2004+" else 0
                d = E[(E.basin == basin) & (E.season >= lo)]
                for variant, xs, pos in (("lagged", ["NAO", "PNA"], False),
                                         ("lagged+position", ["NAO", "PNA"], True),
                                         ("at-deepest", ["NAO_now", "PNA_now"], False)):
                    if variant == "at-deepest":
                        # permutation needs a season-day table for the 'now' index; skip p there
                        dd = d.dropna(subset=[y] + xs)
                        ref = {x: (dd[x].mean(), dd[x].std()) for x in xs}
                        dz = zs(dd, xs, ref)
                        b = ols(dz, y, xs, False)
                        seasons = np.sort(dd.season.unique()); groups = {s: dz[dz.season == s] for s in seasons}
                        bb = np.array([ols(pd.concat([groups[s] for s in rng.choice(seasons, len(seasons))]), y, xs, False)
                                       for _ in range(nboot)])
                        r = dict(n=len(dd), seasons=len(seasons), b=b, se=bb.std(0, ddof=1), p=[np.nan] * 2,
                                 loso_lo=[np.nan] * 2, loso_hi=[np.nan] * 2, nperm=0)
                    else:
                        r = fit_all(d, y, xs, pos, I, s0, rng, nboot)
                    for k, x in enumerate(("NAO", "PNA")):
                        res.append(dict(basin=basin, era=era, outcome=y, variant=variant, index=x,
                                        slope=r["b"][k], se=r["se"][k], t=r["b"][k] / r["se"][k], p_perm=r["p"][k],
                                        loso_lo=r["loso_lo"][k], loso_hi=r["loso_hi"][k],
                                        n=r["n"], seasons=r["seasons"], nperm=r["nperm"]))
                print(basin, era, y, flush=True)
    R = pd.DataFrame(res)
    prim = (R.variant == "lagged") & (((R.era == "all") & ~R.outcome.isin(GUST)) | ((R.era == "2004+") & R.outcome.isin(GUST)))
    R["q_bh"] = np.nan
    R.loc[prim, "q_bh"] = bh(R.loc[prim, "p_perm"].values)
    R.to_csv(os.path.join(out, "results.csv"), index=False, float_format="%.4g")
    keep = ["track", "basin", "season", "start", "tmin", "lat_min", "lon_min", "month", "minp", "clim", "ring",
            "anom_clim", "depth_ring", "bg", "ndr_max", "g800max", "hf_hours", "NAO", "PNA", "NAO_now", "PNA_now", "src"]
    E[keep].to_csv(os.path.join(out, "event_table.csv"), index=False, float_format="%.4g")

    # position shift of the deepest fix (the track-shift part), lagged index
    pos = []
    for basin in ("atl", "pac"):
        d = E[E.basin == basin].copy()
        d["lat"] = d.lat_min; d["lon"] = d.lon_c
        for y in ("lat", "lon"):
            r = fit_all(d.rename(columns={}), y, ["NAO", "PNA"], False, I, s0, rng, nboot)
            for k, x in enumerate(("NAO", "PNA")):
                pos.append(dict(basin=basin, outcome=y, index=x, slope=r["b"][k], se=r["se"][k],
                                t=r["b"][k] / r["se"][k], p_perm=r["p"][k], loso_lo=r["loso_lo"][k], loso_hi=r["loso_hi"][k], n=r["n"]))
    P = pd.DataFrame(pos)
    P.to_csv(os.path.join(out, "position_shift.csv"), index=False, float_format="%.4g")

    icc = {}
    for basin in ("atl", "pac"):
        for y in OUTCOMES:
            d = E[(E.basin == basin) & (E.season >= (GUST_FROM if y in GUST else 0))].dropna(subset=[y])
            g = d.groupby("season")[y]; k = g.ngroups; N = len(d)
            n0 = (N - (g.size() ** 2).sum() / N) / (k - 1)
            msb = (g.size() * (g.mean() - d[y].mean()) ** 2).sum() / (k - 1)
            msw = ((d[y] - g.transform("mean")) ** 2).sum() / (N - k)
            icc[f"{basin}_{y}"] = float((msb - msw) / (msb + (n0 - 1) * msw))
    json.dump(dict(icc=icc, nboot=nboot, n_events=int(len(E)), arco_events=int((E.src == "arco").sum()),
                   missing_background=int(E.ring.isna().sum())), open(os.path.join(out, "meta.json"), "w"), indent=1)

    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write("ERA5 PROXY, pipeline A events. Slopes per SD of the index, NAO and PNA fitted together,\n"
                "month fixed effects + linear season term. se: season-block bootstrap; p: season-block permutation;\n"
                "LOSO: leave-one-season-out range; q: BH over the primary rows (lagged index, no position).\n\n")
        for variant in ("lagged", "lagged+position", "at-deepest"):
            f.write(f"== {variant} ==\n")
            sub = R[R.variant == variant]
            for r in sub.itertuples():
                f.write(f"{r.basin} {r.era:5s} {r.outcome:10s} {r.index}  {r.slope:+7.3f} {UNITS[r.outcome]:8s} se {r.se:.3f}  t {r.t:+6.2f}"
                        f"  p {r.p_perm:.4f}  LOSO [{r.loso_lo:+.3f}, {r.loso_hi:+.3f}]  q {r.q_bh:.3f}  n {r.n} / {r.seasons} seasons\n")
            f.write("\n")
        f.write("== position of the deepest fix, lagged index (deg per SD) ==\n")
        for r in P.itertuples():
            f.write(f"{r.basin} {r.outcome:4s} {r.index} {r.slope:+.3f} se {r.se:.3f} t {r.t:+.2f} p {r.p_perm:.4f} LOSO [{r.loso_lo:+.3f}, {r.loso_hi:+.3f}]\n")
        f.write("\nICC by season: " + json.dumps({k: round(v, 3) for k, v in icc.items()}) + "\n")
    print(open(os.path.join(out, "summary.txt")).read())


if __name__ == "__main__":
    main()
