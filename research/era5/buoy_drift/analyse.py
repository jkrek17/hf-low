"""Does ERA5's surface wind drift against moored buoys before 2001?

ERA5 is a PROXY record here as everywhere in this project. This script tests
one candidate cause of the gust drift measured in research/era5/drift_check.py
(pipeline A, research/era5/hf_history: +0.665 kt/decade at fixed depth over
1979-2000, mean 64.38 kt; pipeline B, event_fields/criterion/series: +0.645
kt/decade, mean 58.40 kt; both about +1.0 to +1.1 % per decade): that ERA5's
own 10 m winds and gusts strengthen over 1979-2000 relative to the real wind.
The reference is in-situ: moored-buoy wind reports from NOAA's ISD
(fetch_buoys.py), matched to ERA5 at the nearest grid point (extract_era5.py).

Definitions
  pair      one buoy report at 00/06/12/18 UTC, October to March, matched to
            the ERA5 analysis at the same time and the nearest 0.25 deg point.
  season    cold season named by the year it starts (Oct 1985 - Mar 1986 = 1985),
            as in the June-May seasons elsewhere in the project.
  U_b       buoy wind (8-min mean for NDBC hulls) adjusted from an assumed
            5 m anemometer height to 10 m with a neutral log profile,
            z0 = 2e-4 m (factor 1.064). The per-deployment hull history was not
            reachable, so this is a uniform assumption; it scales each
            station by a constant and cannot create a trend (see README.md).
  U_e, G_e  ERA5 10 m wind speed and instantaneous 10 m gust.
  r_wind    ln(U_e / U_b): relative ERA5 wind bias (x100 ~ percent).
  r_gust    ln(G_e / U_b): ERA5 gust against the in-situ wind. Its trend is the
            gust-index drift as far as it is visible at the buoys.
  r_gf      ln(G_e / U_e): ERA5's own gust factor (r_gust = r_wind + r_gf).
  bins      by the pair mean (U_e + U_b)/2, which conditions on neither side
            alone and so does not manufacture a regression-to-the-mean bias:
            >= 5, >= 10, >= 15 and >= 17.5 m/s (17.5 m/s = 34 kt, gale force).

Trend model, per bin: one value per station-season (mean over its pairs,
needs >= 20 pairs, >= 8 in the two upper bins), regressed on season with a
fixed effect per station (or per station segment, see below), slope per
decade, standard error clustered by season, so the sample size behind every t
is seasons. Leave-one-season-out and leave-one-station-out slope ranges are
reported beside it.

Station-specific breaks (hull, anemometer height or sensor changes) are found
in each station's monthly r_wind after removing the cross-station median of
the same month, so that a drift common to all stations (which is what an ERA5
drift would be) is not mistaken for a break. Binary segmentation, accepting
a step only at |t| >= 5, |step| >= 0.03 and >= 6 months each side.

    python3 research/era5/buoy_drift/analyse.py   # writes results/buoy_drift-result.txt
                                                  # and results/station_seasons.csv
"""
import csv
import gzip
import glob
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from drift_check import load_b, t_cdf, t_crit  # noqa: E402  (numpy-only t distribution)
from fetch_buoys import STATIONS  # noqa: E402

WORK = os.path.join(HERE, "work")
RES = os.path.join(HERE, "results")
Z0 = 2e-4
HEIGHT_FACTOR = math.log(10 / Z0) / math.log(5 / Z0)
MS2KT = 1 / 0.514444
PASS_QC = set("01459")
STORM_KM = 500
UNIT_SCREEN = 0.35            # ln(1.42); the knot error is ln(1.944) = 0.66
BINS = [("all pairs", 0.0), (">= 5 m/s", 5.0), (">= 10 m/s", 10.0), (">= 15 m/s", 15.0), (">= 17.5 m/s (gale)", 17.5)]
MIN_PAIRS = {0.0: 20, 5.0: 20, 10.0: 20, 15.0: 8, 17.5: 8}
PRE = (1979, 2000)
POST = (2001, 2004)          # 2004 = Oct-Dec 2004 only; ISD has no moored buoys after 2004
# drift_check-result.txt section 3 and 5, mean ERA5 gust at fixed depth 955-975 hPa:
PIPE = {"A (hf_history, g800)": (0.665, 64.38), "B (event_fields, gust500)": (0.645, 58.40)}


# --- data ----------------------------------------------------------------------

def load_era5():
    files = sorted(f for f in glob.glob(os.path.join(WORK, "era5", "[12]*.npz")))
    t, ws, gu = [], [], []
    for f in files:
        d = np.load(f)
        t.append(d["time"])
        ws.append(d["ws"])
        gu.append(d["gust"])
        idx = d["idx"]
    lsm = np.load(os.path.join(WORK, "era5", "lsm.npz"))
    assert np.array_equal(lsm["idx"], idx)
    return np.concatenate(t), np.concatenate(ws), np.concatenate(gu), lsm["lsm"], len(files)


def load_buoy(sid):
    out = {}
    with gzip.open(os.path.join(WORK, "isd", sid + ".csv.gz"), "rt") as f:
        for r in csv.DictReader(f):
            d = r["DATE"]
            if d[14:19] != "00:00" or d[11:13] not in ("00", "06", "12", "18"):
                continue
            if int(d[5:7]) not in (10, 11, 12, 1, 2, 3):
                continue
            w = r["WND"].split(",")
            if len(w) != 5 or w[4] not in PASS_QC:
                continue
            if w[2] == "C":
                spd = 0.0
            elif w[3] == "9999" or w[2] == "9":
                continue
            else:
                spd = int(w[3]) / 10.0
            out[np.datetime64(d[:13], "h")] = spd
    return out


def season_of(t):
    y = t.astype("datetime64[Y]").astype(int) + 1970
    m = t.astype("datetime64[M]").astype(int) % 12 + 1
    return np.where(m >= 6, y, y - 1), m


# --- statistics ----------------------------------------------------------------

def fe_slope(y, season, groups, cluster):
    """Slope of y on season (per decade) with a fixed effect per group,
    se clustered by `cluster` (CR1)."""
    ug, gi = np.unique(groups, return_inverse=True)
    x = (season - 1990) / 10.0
    X = np.column_stack([x, np.eye(len(ug))[gi]])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    bread = np.linalg.pinv(X.T @ X)
    meat = np.zeros((X.shape[1],) * 2)
    uc = np.unique(cluster)
    for c in uc:
        m = cluster == c
        s = X[m].T @ r[m]
        meat += np.outer(s, s)
    G, n, k = len(uc), len(y), X.shape[1]
    cov = G / (G - 1) * (n - 1) / max(n - k, 1) * bread @ meat @ bread
    b, se = beta[0], math.sqrt(cov[0, 0])
    t = b / se
    return dict(b=b, se=se, t=t, p=2 * (1 - t_cdf(abs(t), G - 1)), G=G, n=n, crit=t_crit(G - 1))


def station_slopes(c, lo, hi, min_seasons=8):
    """One slope per station (segment intercepts within the station), then the
    mean across stations with its se: inference that lets each station have its
    own buoy-side trend, which a season-clustered se does not."""
    out = []
    for k in np.unique(c["k"]):
        m = (c["k"] == k) & (c["s"] >= lo) & (c["s"] <= hi)
        if len(np.unique(c["s"][m])) < min_seasons:
            continue
        ug, gi = np.unique(c["g"][m], return_inverse=True)
        X = np.column_stack([(c["s"][m] - 1990) / 10.0, np.eye(len(ug))[gi]])
        beta, *_ = np.linalg.lstsq(X, c["y"][m], rcond=None)
        out.append(beta[0])
    out = np.array(out)
    n = len(out)
    se = out.std(ddof=1) / math.sqrt(n)
    t = out.mean() / se
    return dict(mean=out.mean(), se=se, t=t, p=2 * (1 - t_cdf(abs(t), n - 1)), n=n,
                median=float(np.median(out)), lo=out.min(), hi=out.max(), neg=int((out < 0).sum()))


def season_effects(y, season, station):
    """Two-way fixed effects: station + season dummies; season effects centred."""
    us, si = np.unique(station, return_inverse=True)
    uy, yi = np.unique(season, return_inverse=True)
    X = np.column_stack([np.eye(len(us))[si], np.eye(len(uy))[yi][:, 1:]])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    eff = np.concatenate([[0.0], beta[len(us):]])
    return uy, eff - eff.mean(), np.bincount(yi)


def segment(d, t_min=5.0, step_min=0.03, n_min=6, depth=4):
    """Binary segmentation of a station's de-commoned monthly series."""
    cuts = []

    def rec(lo, hi, k):
        if k == 0 or hi - lo < 2 * n_min:
            return
        best = None
        for c in range(lo + n_min, hi - n_min + 1):
            a, b = d[lo:c], d[c:hi]
            se = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
            t = (b.mean() - a.mean()) / se if se > 0 else 0.0
            if best is None or abs(t) > abs(best[1]):
                best = (c, t, b.mean() - a.mean())
        if best and abs(best[1]) >= t_min and abs(best[2]) >= step_min:
            cuts.append(best)
            rec(lo, best[0], k - 1)
            rec(best[0], hi, k - 1)

    rec(0, len(d), depth)
    return sorted(cuts)


def storm_section(P, *_):
    """Pairs with a pipeline B cyclone centre within STORM_KM, at the same time."""
    lat0 = {}
    for k, sid in enumerate(STATIONS):
        with gzip.open(os.path.join(WORK, "isd", sid + ".csv.gz"), "rt") as f:
            r = next(csv.DictReader(f))
        lat0[k] = (float(r["LATITUDE"]), float(r["LONGITUDE"]))
    _, _, tt, st, sea, sids, r_wind, r_gust, r_gf, good = _
    moments = {}
    with gzip.open(os.path.join(os.path.dirname(HERE), "hf_probability.csv.gz"), "rt") as fh:
        for r in csv.DictReader(fh):
            moments.setdefault(np.datetime64(r["time_utc"][:13], "h"), []).append(
                (float(r["lat"]), float(r["lon"]), float(r["centre_hpa"])))
    P("10. NEAR DEEP LOWS: pairs with a pipeline B cyclone centre within %d km at the same time" % STORM_KM)
    P("    (B's scored moments, centre <= 990 hPa, Sep-May; pair-level fit, station FE, se clustered by season)")
    for lab, lo_p, hi_p in (("centre 955-975 hPa", 955, 975), ("centre <= 990 hPa", 0, 990)):
        near = np.zeros(len(st), bool)
        for i, (t, k) in enumerate(zip(tt, st)):
            for la, lo, c in moments.get(t, ()):
                if lo_p <= c <= hi_p:
                    a1, b1 = math.radians(lat0[k][0]), math.radians(lat0[k][1])
                    a2, b2 = math.radians(la), math.radians(lo)
                    d = 6371 * math.acos(min(1, math.sin(a1) * math.sin(a2) + math.cos(a1) * math.cos(a2) * math.cos(b1 - b2)))
                    if d <= STORM_KM:
                        near[i] = True
                        break
        for name, val in (("r_wind", r_wind), ("r_gust", r_gust), ("r_gf", r_gf)):
            for a, b in ((PRE[0], PRE[1]), (1985, PRE[1])):
                m = near & good & (sea >= a) & (sea <= b)
                r = fe_slope(val[m], sea[m], st[m], sea[m])
                P("   %-20s %-7s %d-%d  %+7.3f %%/decade (se %.3f) t=%+5.2f p=%.3f  pairs %d seasons %d stations %d  mean %+.3f" % (
                    lab, name, a, b, 100 * r["b"], 100 * r["se"], r["t"], r["p"], m.sum(), r["G"],
                    len(np.unique(st[m])), val[m].mean()))
        m = near & good
        P("   %-20s pairs by epoch: 1979-1990 %d, 1991-2000 %d, 2001-2004 %d" % (
            lab, (m & (sea <= 1990)).sum(), (m & (sea > 1990) & (sea <= 2000)).sum(), (m & (sea > 2000)).sum()))


# --- main ------------------------------------------------------------------------

def main():
    os.makedirs(RES, exist_ok=True)
    lines = []
    P = lambda s="": lines.append(s)
    T, WS, GU, LSM, nfiles = load_era5()
    sea_all, _ = season_of(T)
    order = np.argsort(T)
    T, WS, GU = T[order], WS[order], GU[order]
    tpos = {t: i for i, t in enumerate(T)}
    sids = list(STATIONS)
    lsm = LSM.reshape(len(sids), 9)

    P("ERA5 surface wind against moored buoys, 1979-2004 (research/era5/buoy_drift/analyse.py)")
    P("=" * 96)
    P("ERA5 is a PROXY record. Buoy reports: NOAA ISD global-hourly (FM-18), 00/06/12/18 UTC,")
    P("October-March. ERA5: 10 m wind speed (ARCO 6-hourly store) and instantaneous 10 m gust (the")
    P("variable pipelines A and B use), nearest 0.25 deg point. Buoy wind x %.3f (5 m -> 10 m, log" % HEIGHT_FACTOR)
    P("profile, z0 = 2e-4 m; assumed for every station). Slopes are per decade of season, with")
    P("station fixed effects and season-clustered se; t is judged on seasons - 1 df.")
    P("ERA5 months read: %d (%d analysis times)." % (nfiles, len(T)))
    P()

    # station screen and pairs
    rows = []                          # one per pair
    P("1. STATIONS (nearest-point land-sea mask, max over the 3 x 3 neighbourhood)")
    kept = []
    for k, sid in enumerate(sids):
        label, basin = STATIONS[sid]
        b = load_buoy(sid)
        ok = lsm[k, 4] == 0 and lsm[k].max() < 0.05
        n = 0
        for t, ub in b.items():
            i = tpos.get(t)
            if i is None:
                continue
            ue, ge = float(WS[i, 9 * k + 4]), float(GU[i, 9 * k + 4])
            if ok:
                rows.append((k, t, ub * HEIGHT_FACTOR, ue, ge))
            n += 1
        seas = sorted({int(season_of(np.array([t]))[0][0]) for t in b})
        P("   %s %-8s %s  lsm %.2f / max %.2f  pairs %6d  seasons %s-%s  %s" % (
            sid, label, basin, lsm[k, 4], lsm[k].max(), n, seas[0] if seas else "-", seas[-1] if seas else "-",
            "kept" if ok else "DROPPED (land in the 3 x 3 neighbourhood)"))
        if ok:
            kept.append(k)
    P("   %d of %d stations kept." % (len(kept), len(sids)))
    P()

    st = np.array([r[0] for r in rows])
    tt = np.array([r[1] for r in rows])
    ub = np.array([r[2] for r in rows])
    ue = np.array([r[3] for r in rows])
    ge = np.array([r[4] for r in rows])
    good = (ub >= 3) & (ue >= 3)
    rw_all = np.where(good, np.log(np.maximum(ue, 1e-3) / np.maximum(ub, 1e-3)), np.nan)
    mkey = tt.astype("datetime64[M]")

    # ISD unit errors. For buoys whose GTS reports were coded in knots, ISD holds
    # stretches where the speed is too low by the knot factor (1.944): ERA5/buoy
    # jumps to about 1.9 and drops back to about 1.0 in one month. A station-month
    # whose median ln(U_e/U_b) exceeds UNIT_SCREEN in size is dropped whole.
    P("1b. ISD UNIT-ERROR SCREEN: station-months with |median ln(U_e/U_b)| > %.2f dropped" % UNIT_SCREEN)
    drop = np.zeros(len(st), bool)
    for k in kept:
        bad = []
        for m in np.unique(mkey[st == k]):
            sel = (st == k) & (mkey == m)
            v = rw_all[sel & good]
            if len(v) >= 10 and abs(np.median(v)) > UNIT_SCREEN:
                drop |= sel
                bad.append((m, float(np.median(v))))
        if bad:
            P("   %s %-8s %3d months dropped, %s to %s, median ERA5/buoy %.2f (knot factor 1.944)" % (
                sids[k], STATIONS[sids[k]][0], len(bad), bad[0][0], bad[-1][0], math.exp(np.median([b for _, b in bad]))))
    P("   %d of %d pairs dropped." % (drop.sum(), len(drop)))
    P()
    keep = ~drop
    st, tt, ub, ue, ge, good, rw_all, mkey = (a[keep] for a in (st, tt, ub, ue, ge, good, rw_all, mkey))
    sea, mon = season_of(tt)
    pm = (ub + ue) / 2
    basin = np.array([STATIONS[sids[k]][1] for k in st])

    # monthly series for break detection, all-wind log ratio (pairs with both >= 3 m/s)
    months = np.unique(mkey)
    M = np.full((len(sids), len(months)), np.nan)
    mi = np.searchsorted(months, mkey)
    for k in kept:
        for j in range(len(months)):
            m = (st == k) & (mi == j) & good
            if m.sum() >= 30:
                M[k, j] = rw_all[m].mean()
    seg_of = {}
    P("2. STATION BREAKS in monthly ln(U_e/U_b), cross-station median removed (|t| >= 5, |step| >= 0.03)")
    for k in kept:
        cols = np.where(~np.isnan(M[k]))[0]
        others = np.nanmedian(np.delete(M, k, axis=0)[:, cols], axis=0)
        d = M[k, cols] - others
        ok = ~np.isnan(d)
        cols, d = cols[ok], d[ok]
        cuts = segment(d)
        bounds = [months[cols[c]] for c, _, _ in cuts]
        seg_of[k] = bounds
        P("   %s %-8s %3d months  %s" % (sids[k], STATIONS[sids[k]][0], len(d),
                                       "; ".join("%s step %+.3f (t %+.1f)" % (months[cols[c]], s, t) for c, t, s in cuts) or "no break"))
    segid = np.array(["%d:%d" % (k, int(np.searchsorted(np.array(seg_of[k], "datetime64[M]"), m, side="right")))
                      for k, m in zip(st, mkey)])
    P()

    # station-season cells
    def cells(val, thr, sel=None):
        m = (pm >= thr) & np.isfinite(val)
        if sel is not None:
            m &= sel
        keys = {}
        for k, s, g, v in zip(st[m], sea[m], segid[m], val[m]):
            keys.setdefault((k, s, g), []).append(v)
        out = [(k, s, g, float(np.mean(v)), len(v)) for (k, s, g), v in keys.items() if len(v) >= MIN_PAIRS[thr]]
        if not out:
            return None
        k, s, g, y, n = map(np.array, zip(*out))
        return dict(k=k, s=s, g=g, y=y, n=n)

    with np.errstate(divide="ignore", invalid="ignore"):
        r_wind = np.log(ue / ub)
        r_gust = np.log(ge / ub)
        r_gf = np.log(ge / ue)
    diff_kt = (ue - ub) * MS2KT

    def report(name, val, scale, unit, sel=None, csv_rows=None):
        P("   " + name)
        res = {}
        for lab, thr in BINS:
            v = val if thr > 0 or name.startswith("U_e - U_b") else val
            c = cells(v, thr, sel)
            if c is None:
                continue
            w = (c["s"] >= PRE[0]) & (c["s"] <= PRE[1])
            out = []
            for tag, grp in (("station FE", c["k"]), ("segment FE", c["g"])):
                r = fe_slope(c["y"][w], c["s"][w], grp[w], c["s"][w])
                loso = [fe_slope(c["y"][w & (c["s"] != s)], c["s"][w & (c["s"] != s)], grp[w & (c["s"] != s)],
                                 c["s"][w & (c["s"] != s)])["b"] for s in np.unique(c["s"][w])]
                lost = [fe_slope(c["y"][w & (c["k"] != k)], c["s"][w & (c["k"] != k)], grp[w & (c["k"] != k)],
                                 c["s"][w & (c["k"] != k)])["b"] for k in np.unique(c["k"][w])]
                flag = "  SIGNIFICANT" if abs(r["t"]) > r["crit"] else ""
                out.append(r)
                P("     %-20s %-10s %+7.3f %s/decade (se %.3f) t=%+5.2f p=%.3f  seasons %d cells %d  "
                  "LOSO %+.3f..%+.3f  LOStation %+.3f..%+.3f%s" % (
                      lab, tag, scale * r["b"], unit, scale * r["se"], r["t"], r["p"], r["G"], r["n"],
                      scale * min(loso), scale * max(loso), scale * min(lost), scale * max(lost), flag))
            ss = station_slopes(c, *PRE)
            P("     %-20s across stations (own slope each, segments within): mean %+.3f %s/decade (se %.3f) t=%+5.2f p=%.3f "
              "n=%d stations; median %+.3f, range %+.3f..%+.3f, %d negative" % (
                  "", scale * ss["mean"], unit, scale * ss["se"], ss["t"], ss["p"], ss["n"], scale * ss["median"],
                  scale * ss["lo"], scale * ss["hi"], ss["neg"]))
            # level by epoch from two-way season effects over all seasons
            uy, eff, cnt = season_effects(c["y"], c["s"], c["g"])
            e = lambda a, b: eff[(uy >= a) & (uy <= b)]
            pre, post = e(*PRE), e(*POST)
            if len(post) >= 2 and len(pre) >= 2:
                dlt = post.mean() - pre.mean()
                se = math.sqrt(pre.var(ddof=1) / len(pre) + post.var(ddof=1) / len(post))
                P("     %-20s level 2001-04 minus 1979-2000 (segment FE, season effects): %+.3f %s (se %.3f, t %+.2f, %d vs %d seasons)" % (
                    "", scale * dlt, unit, scale * se, dlt / se, len(post), len(pre)))
            res[lab] = (out, uy, eff, cnt, c)
            if csv_rows is not None:
                for k, s, g, y, n in zip(c["k"], c["s"], c["g"], c["y"], c["n"]):
                    csv_rows.append([name, lab, sids[k], STATIONS[sids[k]][0], STATIONS[sids[k]][1], int(s), g, "%.5f" % y, int(n)])
        return res

    csv_rows = []
    P("3. FULL-DISTRIBUTION BIAS, pooled over kept stations and all seasons (U_b at 10 m)")
    for lab, thr in BINS:
        m = pm >= thr
        P("   %-20s pairs %7d  mean U_e - U_b %+.2f kt  median ln(U_e/U_b) %+.3f  ln(G_e/U_b) %+.3f  ln(G_e/U_e) %+.3f" % (
            lab, m.sum(), diff_kt[m].mean(), np.nanmedian(r_wind[m & good]), np.nanmedian(r_gust[m & good]), np.nanmedian(r_gf[m & good])))
    for a, b in ((1979, 1990), (1991, 2000), (2001, 2004)):
        m = (sea >= a) & (sea <= b)
        h = m & (pm >= 15)
        P("   seasons %d-%d       pairs %7d  mean U_e - U_b %+.2f kt (all), %+.2f kt (>= 15 m/s, %d pairs); stations %d" % (
            a, b, m.sum(), diff_kt[m].mean(), diff_kt[h].mean(), h.sum(), len(np.unique(st[m]))))
    P()
    P("4. TRENDS 1979-2000 (x100 for log ratios, so units are percent per decade)")
    R = {}
    R["wind"] = report("r_wind = ln(U_e/U_b)  [ERA5 10 m wind vs buoy]", r_wind, 100, "%", good, csv_rows)
    R["gust"] = report("r_gust = ln(G_e/U_b)  [ERA5 gust vs buoy wind]", r_gust, 100, "%", good, csv_rows)
    R["gf"] = report("r_gf = ln(G_e/U_e)    [ERA5 internal gust factor]", r_gf, 100, "%", good, csv_rows)
    R["diff"] = report("U_e - U_b in kt        [absolute wind bias]", diff_kt, 1, "kt", None, csv_rows)
    P()
    P("5. BY BASIN, r_gust >= 15 m/s and r_wind >= 15 m/s, segment FE")
    for bs in ("Atl", "Pac"):
        for name, val in (("r_wind", r_wind), ("r_gust", r_gust)):
            c = cells(val, 15.0, good & (basin == bs))
            if c is None:
                continue
            w = (c["s"] >= PRE[0]) & (c["s"] <= PRE[1])
            r = fe_slope(c["y"][w], c["s"][w], c["g"][w], c["s"][w])
            P("   %s %-7s %+7.3f %%/decade (se %.3f) t=%+5.2f p=%.3f seasons %d cells %d" % (
                bs, name, 100 * r["b"], 100 * r["se"], r["t"], r["p"], r["G"], r["n"]))
    P()
    first = {k: sea[st == k].min() for k in np.unique(st)}
    last = {k: sea[st == k].max() for k in np.unique(st)}
    longrec = np.array([first[k] <= 1981 and last[k] >= 1999 for k in st])
    P("5b. LONG-RECORD STATIONS ONLY (seasons from 1981 or earlier to 1999 or later): %s" % ", ".join(
        sorted({STATIONS[sids[k]][0] for k in st[longrec]})))
    for name, val in (("r_wind", r_wind), ("r_gust", r_gust), ("r_gf", r_gf)):
        for thr in (5.0, 15.0):
            c = cells(val, thr, good & longrec)
            w = (c["s"] >= PRE[0]) & (c["s"] <= PRE[1])
            r = fe_slope(c["y"][w], c["s"][w], c["g"][w], c["s"][w])
            P("   %-7s >= %4.1f m/s  %+7.3f %%/decade (se %.3f) t=%+5.2f p=%.3f seasons %d cells %d" % (
                name, thr, 100 * r["b"], 100 * r["se"], r["t"], r["p"], r["G"], r["n"]))
    P()
    P("6. COMPARISON WITH THE GUST-INDEX DRIFT (drift_check-result.txt, 1979-2000, depth 955-975 hPa)")
    for pname, (kt, mean) in PIPE.items():
        P("   pipeline %-26s %+.3f kt/decade on a mean of %.2f kt = %+.2f %%/decade" % (pname, kt, mean, 100 * kt / mean))
    for key, lab in (("gust", ">= 15 m/s"), ("gust", ">= 17.5 m/s (gale)"), ("wind", ">= 15 m/s"), ("gust", ">= 5 m/s")):
        if lab not in R[key]:
            continue
        for r, tag in zip(R[key][lab][0], ("station FE", "segment FE")):
            pct = 100 * r["b"]
            se = 100 * r["se"]
            P("   buoys r_%s %-18s %-10s %+.2f %%/decade (95%% CI %+.2f..%+.2f); as kt/decade on 58.40 kt: %+.3f; "
              "power to detect +1.07 %%/decade %.2f" % (
                  key, lab, tag, pct, pct - r["crit"] * se, pct + r["crit"] * se, pct / 100 * 58.40,
                  0.5 * (1 + math.erf((1.07 / se - r["crit"]) / math.sqrt(2)))))
    P()
    P("7. SEASON EFFECTS (two-way FE: station segment + season), r_gust and r_wind >= 15 m/s, percent")
    for key in ("wind", "gust"):
        if ">= 15 m/s" in R[key]:
            _, uy, eff, cnt, _ = R[key][">= 15 m/s"]
            P("   r_%s: " % key + " ".join("%d:%+.1f(%d)" % (y, 100 * e, n) for y, e, n in zip(uy, eff, cnt)))
    P()

    # Pipeline B's gust at fixed depth, season by season, from the committed
    # per-moment file (the same selection as drift_check.py section 3).
    bs, _, _, bc, _, bg = load_b()
    m = (bc >= 955) & (bc <= 975)
    bseas = np.unique(bs[m])
    bmean = np.array([bg[m & (bs == s)].mean() for s in bseas])
    P("8. SHAPE OF THE DRIFT: buoy season effects against pipeline B's 500 km gust at 955-975 hPa")
    pre = (bseas >= PRE[0]) & (bseas <= PRE[1])
    ref = bmean[pre].mean()
    P("   B, season mean gust500 minus its 1979-2000 mean (kt): " + " ".join(
        "%d:%+.2f" % (s, v - ref) for s, v in zip(bseas[pre], bmean[pre])))
    for key in ("wind", "gust"):
        for lab in (">= 5 m/s", ">= 15 m/s"):
            _, uy, eff, _, _ = R[key][lab]
            for a, b in ((1979, 2000), (1985, 2000)):
                yy = [y for y in range(a, b + 1) if y in set(uy) and y in set(bseas)]
                e = np.array([eff[list(uy).index(y)] for y in yy])
                g = np.array([bmean[list(bseas).index(y)] for y in yy])
                P("   corr(r_%s %-9s season effect, B gust500), %d-%d: %+.2f (n = %d seasons)" % (
                    key, lab, a, b, np.corrcoef(e, g)[0, 1], len(yy)))
    P()
    P("9. SENSITIVITY, CHOSEN AFTER SEEING SECTION 7: seasons 1985-2000 only, segment FE")
    P("   Section 7 and station_seasons.csv show steps of 10-20 % at single stations between 1979 and")
    P("   1985 at different times at different stations (46005 in 1985, MidNomad and 44005 in 1983,")
    P("   46003 in 1980), which a change in ERA5 would not do; the break search caught only some of")
    P("   them. This window drops those seasons. It was picked after looking, so read it as a check.")
    for name, val in (("r_wind", r_wind), ("r_gust", r_gust), ("r_gf", r_gf)):
        for thr in (5.0, 15.0):
            c = cells(val, thr, good)
            w = (c["s"] >= 1985) & (c["s"] <= PRE[1])
            r = fe_slope(c["y"][w], c["s"][w], c["g"][w], c["s"][w])
            loso = [fe_slope(c["y"][w & (c["s"] != s)], c["s"][w & (c["s"] != s)], c["g"][w & (c["s"] != s)],
                             c["s"][w & (c["s"] != s)])["b"] for s in np.unique(c["s"][w])]
            P("   %-7s >= %4.1f m/s  %+7.3f %%/decade (se %.3f) t=%+5.2f p=%.3f seasons %d cells %d  LOSO %+.3f..%+.3f" % (
                name, thr, 100 * r["b"], 100 * r["se"], r["t"], r["p"], r["G"], r["n"], 100 * min(loso), 100 * max(loso)))
            ss = station_slopes(c, 1985, PRE[1])
            P("   %-7s >= %4.1f m/s  across stations: mean %+.3f %%/decade (se %.3f) t=%+5.2f p=%.3f n=%d; median %+.3f, "
              "range %+.3f..%+.3f, %d negative" % (name, thr, 100 * ss["mean"], 100 * ss["se"], ss["t"], ss["p"], ss["n"],
                                                   100 * ss["median"], 100 * ss["lo"], 100 * ss["hi"], ss["neg"]))
    for a, b in ((1979, 2000), (1985, 2000)):
        w = (bseas >= a) & (bseas <= b)
        x = bseas[w].astype(float)
        y = 100 * np.log(bmean[w])
        xm = x - x.mean()
        sl = (xm * (y - y.mean())).sum() / (xm ** 2).sum()
        res = y - y.mean() - sl * xm
        se = math.sqrt((res ** 2).sum() / (len(x) - 2) / (xm ** 2).sum())
        P("   pipeline B gust500 at 955-975 hPa, %d-%d: %+.3f %%/decade (se %.3f) t=%+5.2f, n = %d seasons" % (
            a, b, 10 * sl, 10 * se, sl / se, len(x)))
    P()
    storm_section(P, bs, bc, tt, st, sea, sids, r_wind, r_gust, r_gf, good)
    with open(os.path.join(RES, "buoy_drift-result.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(RES, "station_seasons.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "bin", "isd_usaf", "label", "basin", "season", "segment", "mean", "pairs"])
        w.writerows(csv_rows)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
