"""Does the ERA5 gust drift before 2001 reach thresholded event counts?

Recomputes, from committed files only, the three findings the session behind
c9dc994 reported without a result file, reproduces its event-level
reconciliation, and then tests both ERA5 pipelines for what that
reconciliation did not test directly: whether a thresholded count of events
drifts before 2001 when storm depth is held fixed.

Both ERA5 records are PROXIES. Nothing before the archive (2001-02, complete
from RECORD_START = 2004-05) can be validated; the tests before 2001 here are
internal-consistency tests of the reanalysis, not validation.

Inputs (all committed; nothing is downloaded):
  data/hf_lows/HF_Data_-_{Atl,Pac}.csv            archive, from the spreadsheet
  docs/data/hf-lows.json                          the built archive payload
  research/era5/hf_probability.csv.gz             pipeline B, every scored 6-hourly
                                                  cyclone moment (centre <= 990 hPa),
                                                  Sep-May, 1979-2025
  research/era5/criterion-result.txt              pipeline B report (for comparison)
  research/era5/hf_history/results/*              pipeline A catalog, fixes, season table

Pipeline A: research/era5/hf_history/ (800 km track-max gust index, 71.7 kt,
June-May seasons). Pipeline B: event_fields.py / criterion.py / series.py
(500 km gust per moment, Sep-May seasons). Season = the year the season starts.

    python3 research/era5/drift_check.py   # writes drift_check-result.txt beside it
"""
import collections
import csv
import gzip
import io
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
A_DIR = os.path.join(HERE, "hf_history", "results")
THR_A = 71.7
PRE = list(range(1979, 2001))      # 22 seasons with no archive
OVL = list(range(2004, 2026))      # 22 seasons from RECORD_START


# --- statistics -----------------------------------------------------------------

def t_cdf(t, df):
    """Student t CDF by Simpson integration of the density (no scipy here)."""
    if t == 0:
        return 0.5
    c = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
    f = lambda x: c * (1 + x * x / df) ** (-(df + 1) / 2)
    a, b, n = 0.0, abs(t), 2000
    h = (b - a) / n
    s = f(a) + f(b) + sum((4 if i % 2 else 2) * f(a + i * h) for i in range(1, n))
    area = s * h / 3
    return 0.5 + area if t > 0 else 0.5 - area


def t_crit(df, q=0.975):
    lo, hi = 0.0, 20.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if t_cdf(mid, df) < q else (lo, mid)
    return (lo + hi) / 2


def ols(x, y):
    """Slope per decade, its se, t, two-sided p, n; plain OLS on season."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    xm = x - x.mean()
    sxx = (xm ** 2).sum()
    b = (xm * (y - y.mean())).sum() / sxx
    r = y - y.mean() - b * xm
    se = math.sqrt((r ** 2).sum() / (n - 2) / sxx)
    t = b / se if se > 0 else float("nan")
    p = 2 * (1 - t_cdf(abs(t), n - 2))
    return dict(b=10 * b, se=10 * se, t=t, p=p, n=n, mean=float(y.mean()), crit=t_crit(n - 2))


def ols_multi(X, y):
    """OLS with an intercept; returns coefficients and their se."""
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    dof = len(y) - X.shape[1]
    cov = (r @ r / dof) * np.linalg.inv(X.T @ X)
    return beta, np.sqrt(np.diag(cov))


def cluster_se_slope(x, y, groups, covars=()):
    """Slope of y on x (with covariates) and a season-clustered se, so that the
    sample size behind the se is seasons, not cyclones."""
    X = np.column_stack([np.ones(len(y)), x] + list(covars))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    bread = np.linalg.inv(X.T @ X)
    meat = np.zeros((X.shape[1], X.shape[1]))
    ug = np.unique(groups)
    for g in ug:
        m = groups == g
        s = X[m].T @ r[m]
        meat += np.outer(s, s)
    G = len(ug)
    k = X.shape[1]
    adj = G / (G - 1) * (len(y) - 1) / (len(y) - k)
    cov = adj * bread @ meat @ bread
    return beta[1], math.sqrt(cov[1, 1]), G


def line(P, label, r, unit="/decade"):
    flag = "  SIGNIFICANT" if abs(r["t"]) > r["crit"] else ""
    P("   %-60s n=%2d  %+8.3f %s (se %.3f)  t=%+5.2f  p=%.3f  mean %.2f%s" % (
        label, r["n"], r["b"], unit, r["se"], r["t"], r["p"], r["mean"], flag))


# --- inputs -----------------------------------------------------------------------

def archive_from_json():
    d = json.load(open(os.path.join(ROOT, "docs", "data", "hf-lows.json")))
    lf = {k: i for i, k in enumerate(d["lowFields"])}
    ff = {k: i for i, k in enumerate(d["fixFields"])}
    per = collections.defaultdict(lambda: [0, 0])
    per_low = collections.defaultdict(lambda: [0, 0])
    for L in d["lows"]:
        hf = sum(1 for f in L[lf["fixes"]] if f[ff["cat"]] == "HF")
        if not hf:
            continue
        s = L[lf["season"]]
        per[s][0] += 1
        per[s][1] += hf
        if L[lf["cls"]] == "low":
            per_low[s][0] += 1
            per_low[s][1] += hf
    return per, per_low


def archive_from_csv():
    per = collections.defaultdict(lambda: [0, 0])
    for b in ("Atl", "Pac"):
        ev = collections.Counter()
        for r in csv.DictReader(open(os.path.join(ROOT, "data", "hf_lows", "HF_Data_-_%s.csv" % b))):
            if r["Category"].strip() == "HF":
                ev[r["ID"].strip()] += 1
        for k, n in ev.items():
            s = int(k[:4])
            per[s][0] += 1
            per[s][1] += n
    return per


def load_b():
    rows = []
    with gzip.open(os.path.join(HERE, "hf_probability.csv.gz"), "rt") as fh:
        for r in csv.DictReader(fh):
            rows.append((int(r["season"]), r["basin"], r["track"], float(r["centre_hpa"]),
                         float(r["pmin150_hpa"]), float(r["gust500_kt"])))
    season = np.array([r[0] for r in rows])
    basin = np.array([r[1] for r in rows])
    track = np.array([r[2] for r in rows])
    centre = np.array([r[3] for r in rows])
    pmin = np.array([r[4] for r in rows])
    g500 = np.array([r[5] for r in rows])
    return season, basin, track, centre, pmin, g500


def b_tracks(season, basin, track, centre, g500):
    """Per linked track over its scored moments: season of first moment, basin,
    deepest smoothed centre, and the largest 500 km sea gust."""
    T = {}
    for s, b, k, c, g in zip(season, basin, track, centre, g500):
        t = T.get(k)
        if t is None:
            T[k] = [s, b, c, g]
        else:
            t[2] = min(t[2], c)
            t[3] = max(t[3], g)
    ks = list(T)
    return (np.array([T[k][0] for k in ks]), np.array([T[k][1] for k in ks]),
            np.array([T[k][2] for k in ks]), np.array([T[k][3] for k in ks]))


def load_a():
    cat = list(csv.DictReader(open(os.path.join(A_DIR, "era5_hf_catalog.csv"))))
    fixes = list(csv.DictReader(open(os.path.join(A_DIR, "era5_hf_catalog_tracks.csv"))))
    counts = list(csv.DictReader(open(os.path.join(A_DIR, "era5_hf_counts_by_season.csv"))))
    return cat, fixes, counts


# --- the report -------------------------------------------------------------------

def run():
    buf = io.StringIO()

    def P(*a):
        print(*a, file=buf)
        print(*a, flush=True)

    P("Gust drift before 2001: index value or thresholded counts? (drift_check.py)")
    P("=" * 78)
    P("Both ERA5 records are PROXIES and cannot be validated before the archive. Every")
    P("test before 2001 below is internal to the reanalysis. OLS on season, slopes per")
    P("decade, t against the two-sided 95% critical value for n-2 df. n is seasons.")
    P("Pipeline A = research/era5/hf_history (800 km track-max gust, 71.7 kt, Jun-May).")
    P("Pipeline B = event_fields/criterion/series (500 km gust per moment, Sep-May).")
    P("")

    # ---------------------------------------------------------------- finding 1
    P("1. ARCHIVE RECORDING PRACTICE: HF fixes per event, seasons 2004-2025")
    P("   Reported: -0.181 /decade, t = -2.54 (session, from docs/data/hf-lows.json);")
    P("   c9dc994 / criterion-result.txt: -0.20 /decade, t = -2.53 (in-domain 'low' events).")
    pj, pjl = archive_from_json()
    pc = archive_from_csv()
    for lab, per in (("docs/data/hf-lows.json, every event with an HF fix", pj),
                     ("docs/data/hf-lows.json, class 'low' only", pjl),
                     ("data/hf_lows CSVs, raw archive IDs keyed by (basin, ID)", pc)):
        P("   -- %s" % lab)
        for lo in (2004, 2006):
            ss = [s for s in sorted(per) if lo <= s <= 2025]
            line(P, "events per season, %d-2025" % lo, ols(ss, [per[s][0] for s in ss]))
            line(P, "HF fixes per season, %d-2025" % lo, ols(ss, [per[s][1] for s in ss]))
            line(P, "HF fixes PER EVENT, %d-2025" % lo, ols(ss, [per[s][1] / per[s][0] for s in ss]))
    ss = [s for s in sorted(pj) if 2004 <= s <= 2025]
    P("   per season (json, all classes): season events HF-fixes per-event")
    for s in ss:
        P("     %d  %4d  %5d  %.2f" % (s, pj[s][0], pj[s][1], pj[s][1] / pj[s][0]))
    P("")

    # ---------------------------------------------------------------- finding 2
    cat, fixes, counts = load_a()
    P("2. RECORD_START = 2004: pipeline A's season table in the archive's first seasons")
    P("   Reported: 2001-02 holds 1 archive event against 85 ERA5 events.")
    P("   season    ERA5 Atl  ERA5 Pac  ERA5 both | archive Atl  archive Pac  both   ERA5/archive")
    for r in counts:
        y = int(r["season"][:4])
        if 2001 <= y <= 2006:
            e = int(r["atl_era5"]) + int(r["pac_era5"])
            a = int(r["atl_archive"]) + int(r["pac_archive"])
            P("   %s   %8s  %8s  %9d | %11s  %11s  %4d   %s" % (
                r["season"], r["atl_era5"], r["pac_era5"], e, r["atl_archive"], r["pac_archive"], a,
                "%.2f" % (e / a) if a else "inf"))
    P("   2001-2006 transfer row in skill.txt: n_obs 205, n_fc 421, bias 2.05 (seasons 2001-02..2005-06).")
    P("   The archive totals above for those seasons sum to %d; the shortfall is in 2001-02 to 2003-04." % sum(
        int(r["atl_archive"]) + int(r["pac_archive"]) for r in counts if 2001 <= int(r["season"][:4]) <= 2005))
    P("")

    # ---------------------------------------------------------------- finding 3 (pipeline B)
    season, basin, track, centre, pmin, g500 = load_b()
    P("3. PIPELINE B: ERA5 500 km sea gust at fixed storm depth, per cyclone moment")
    P("   Reported: +0.65 kt/decade, t = +2.40 over 1979-2000 (session check, band on the smoothed")
    P("   centre, all candidates); criterion-result.txt: +0.64, t = 2.45 (band on pmin, scored rows).")
    P("   hf_probability.csv.gz holds the scored rows only (centre <= 990 hPa): %d rows." % len(season))
    for blab, pv in (("smoothed centre (centre_hpa)", centre), ("pmin within 150 km (pmin150_hpa)", pmin)):
        band = (pv >= 955) & (pv <= 975)
        for lo, hi in ((1979, 2000), (2004, 2025), (1979, 2025)):
            ss = list(range(lo, hi + 1))
            line(P, "mean gust500, %s 955-975, %d-%d" % (blab, lo, hi),
                 ols(ss, [g500[(season == s) & band].mean() for s in ss]), "kt/decade")
    for lo, hi in ((1979, 2000), (2004, 2025)):
        ss = list(range(lo, hi + 1))
        line(P, "moments with centre <= 960 hPa, %d-%d" % (lo, hi), ols(ss, [((season == s) & (centre <= 960)).sum() for s in ss]))
        line(P, "moments with pmin <= 960 hPa, %d-%d" % (lo, hi), ols(ss, [((season == s) & (pmin <= 960)).sum() for s in ss]))
        line(P, "moments with gust500 >= 64 kt (scored rows), %d-%d" % (lo, hi), ols(ss, [((season == s) & (g500 >= 64)).sum() for s in ss]))
    P("   (The session's +22.5 and +56.6 counted every candidate, including centre > 990 hPa, which")
    P("    are not in the committed file; the counts here are over scored rows and differ for that reason.)")
    P("")

    # ---------------------------------------------------------------- reconciliation
    P("4. THE RECONCILIATION AS THE SESSION RAN IT (pipeline A season counts)")
    yr = lambda r: int(r["season"][:4])
    for key, lab in (("atl_era5", "Atlantic"), ("pac_era5", "Pacific")):
        for lo, hi in ((1979, 1996), (1997, 2025), (1979, 2000), (1979, 2025)):
            rr = [r for r in counts if lo <= yr(r) <= hi]
            line(P, "pipeline A %s events per season, %d-%d" % (lab, lo, hi),
                 ols([yr(r) for r in rr], [int(r[key]) for r in rr]))
    for lo in (2001, 2004, 2006):
        rr = [r for r in counts if r["atl_archive"] and yr(r) >= lo]
        line(P, "pipeline A ERA5 minus archive, both basins, %d-2025" % lo,
             ols([yr(r) for r in rr], [int(r["atl_era5"]) + int(r["pac_era5"]) - int(r["atl_archive"]) - int(r["pac_archive"]) for r in rr]))
    P("   What this tests: whether A's counts ramp before 1997 and whether A tracks the archive from")
    P("   2004. What it does not test: whether gust drift at fixed depth moves A's counts. A flat")
    P("   count series is also what a real decline in storms plus a gust drift would give. Section 7")
    P("   asks whether this test could have seen the bias the drift implies.")
    P("")

    # ---------------------------------------------------------------- event level, pipeline B
    ts, tb, tc, tg = b_tracks(season, basin, track, centre, g500)
    P("5. PIPELINE B AT EVENT LEVEL: tracks, not moments")
    P("   %d linked tracks with at least one scored moment. Track index = max gust500 over its" % len(ts))
    P("   scored moments; track depth = its deepest smoothed centre.")
    pre = np.isin(ts, PRE)
    ovl = np.isin(ts, OVL)
    P("   5a. Track index at fixed track depth (track depth 955-975 hPa)")
    band = (tc >= 955) & (tc <= 975)
    for lo, hi in ((1979, 2000), (2004, 2025)):
        ss = list(range(lo, hi + 1))
        line(P, "mean track-max gust500, track depth 955-975, %d-%d" % (lo, hi),
             ols(ss, [tg[(ts == s) & band].mean() for s in ss]), "kt/decade")
    P("   5b. Track index on season with depth held fixed by regression (tracks deeper than 990 hPa")
    P("       excluded by construction; covariates: depth, depth^2, basin). se clustered on season.")
    drift_b = {}
    for lo, hi in ((1979, 2000), (2004, 2025)):
        m = (ts >= lo) & (ts <= hi)
        d0 = tc[m] - 970.0
        b, se, G = cluster_se_slope(ts[m].astype(float), tg[m], ts[m], [d0, d0 ** 2, (tb[m] == "pac").astype(float)])
        drift_b[(lo, hi)] = (10 * b, 10 * se)
        P("   %-60s G=%2d  %+8.3f kt/decade (se %.3f)  t=%+5.2f" % ("track index | depth, %d-%d" % (lo, hi), G, 10 * b, 10 * se, b / se))
    P("       The same regression within track-depth bands. The pooled slope is weighted toward the")
    P("       many shallow tracks; events near the threshold are mostly 955-975 hPa deep.")
    for dlo, dhi in ((900, 955), (955, 965), (965, 975), (975, 990)):
        for lo, hi in ((1979, 2000), (2004, 2025)):
            m = (ts >= lo) & (ts <= hi) & (tc > dlo) & (tc <= dhi)
            d0 = tc[m] - (dlo + dhi) / 2.0
            b, se, G = cluster_se_slope(ts[m].astype(float), tg[m], ts[m], [d0, (tb[m] == "pac").astype(float)])
            drift_b[(lo, hi, dlo, dhi)] = (10 * b, 10 * se)
            P("   %-60s G=%2d  %+8.3f kt/decade (se %.3f)  t=%+5.2f  tracks %d" % (
                "track index | depth, depth %d-%d hPa, %d-%d" % (dlo, dhi, lo, hi), G, 10 * b, 10 * se, b / se, m.sum()))
    P("   5c. Thresholded event counts, depth-standardised. For a gust cut, each season's expected")
    P("       count is the sum over its tracks of the exceedance rate for the track's depth bin (2 hPa")
    P("       bins), rates taken from a reference period. Observed minus expected is the part of the")
    P("       count that depth does not explain. If the gust drift reaches counts, it trends upward")
    P("       before 2001. Depth-only counts (tracks deeper than 960 hPa) are shown alongside.")
    bins = np.floor((tc - 900) / 2).astype(int)
    for cut in (64.0, 69.5, 72.3):
        ex = tg >= cut
        for refl, refm in (("ref 1979-2000", pre), ("ref 2004-2025", ovl)):
            rate = {}
            for k in np.unique(bins):
                mm = refm & (bins == k)
                rate[k] = ex[mm].mean() if mm.sum() >= 5 else np.nan
            # bins too thin in the reference fall back to the pooled rate
            pooled = {k: ex[bins == k].mean() for k in np.unique(bins)}
            exp_t = np.array([rate[k] if not np.isnan(rate[k]) else pooled[k] for k in bins])
            for lo, hi in ((1979, 2000), (2004, 2025)):
                ss = list(range(lo, hi + 1))
                obs = [ex[ts == s].sum() for s in ss]
                exp = [exp_t[ts == s].sum() for s in ss]
                if refl == "ref 1979-2000":
                    if lo == 1979:
                        line(P, "events, track index >= %.1f kt, %d-%d" % (cut, lo, hi), ols(ss, obs))
                    else:
                        line(P, "events, track index >= %.1f kt, %d-%d" % (cut, lo, hi), ols(ss, obs))
                line(P, "  observed - expected | depth (%s), %d-%d" % (refl, lo, hi),
                     ols(ss, [o - e for o, e in zip(obs, exp)]))
    for lo, hi in ((1979, 2000), (2004, 2025)):
        ss = list(range(lo, hi + 1))
        line(P, "depth-only: tracks deeper than 960 hPa, %d-%d" % (lo, hi), ols(ss, [((ts == s) & (tc <= 960)).sum() for s in ss]))
    P("")

    # ---------------------------------------------------------------- pipeline A's own index
    P("6. PIPELINE A'S OWN INDEX AGAINST THE DRIFT")
    P("   The committed catalog holds every track at or above 71.7 kt (events, complete) and one")
    P("   matched sub-threshold track per event (null cases, a sample). Sub-threshold tracks that are")
    P("   not null cases are not committed, so a fixed-depth test of A's index cannot use the whole")
    P("   population; each test below says what selection it carries.")
    ev = [r for r in cat if r["role"] == "event"]
    nu = [r for r in cat if r["role"] == "null_case"]
    role = {r["track"]: r["role"] for r in cat}
    tseason = {r["track"]: int(r["season"]) for r in cat}
    P("   events %d, null cases %d. Null share by track depth (all seasons):" % (len(ev), len(nu)))
    for lo, hi in ((0, 950), (950, 960), (960, 970), (970, 980), (980, 990), (990, 1010)):
        ne = sum(1 for r in ev if lo < float(r["minp"]) <= hi)
        nn = sum(1 for r in nu if lo < float(r["minp"]) <= hi)
        P("     minp %4d-%4d hPa: events %4d  nulls %4d" % (lo, hi, ne, nn))
    # 6a per-fix, fixed depth
    P("   6a. Per fix: A's 800 km gust (g800) at fixed MSLP 955-975 hPa, in-domain fixes of")
    P("       catalog tracks. Selection: tracks with index >= 71.7 are all present; shallower")
    P("       sub-threshold tracks only as null cases.")
    fx = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in fixes:
        if not r["basin"]:
            continue
        msl = float(r["msl"])
        if 955 <= msl <= 975:
            tr = r["track"]
            s = tseason.get(tr)
            if s is None:
                continue
            fx[role[tr]][s].append(float(r["g800"]))
            fx["all"][s].append(float(r["g800"]))
    for rl in ("all", "event"):
        for lo, hi in ((1979, 2000), (2004, 2025)):
            ss = list(range(lo, hi + 1))
            line(P, "mean g800 at 955-975 hPa, %s tracks, %d-%d" % ("catalog" if rl == "all" else "event", lo, hi),
                 ols(ss, [np.mean(fx[rl][s]) for s in ss]), "kt/decade")
    # 6b track level, deep tracks where nulls are rare
    P("   6b. Per track: index on season with depth held fixed, tracks deeper than 960 hPa (where")
    P("       null cases are about a tenth of catalog tracks, so the truncation at 71.7 kt bites")
    P("       least). Covariates depth, depth^2, basin; se clustered on season.")
    drift_a = {}
    for rl, rows in (("events + nulls", ev + nu), ("events only", ev)):
        for lo, hi in ((1979, 2000), (2004, 2025)):
            rr = [r for r in rows if lo <= int(r["season"]) <= hi and float(r["minp"]) <= 960]
            x = np.array([int(r["season"]) for r in rr], float)
            y = np.array([float(r["gust800_kt"]) for r in rr])
            d0 = np.array([float(r["minp"]) for r in rr]) - 950
            bs = np.array([r["basin"] == "pac" for r in rr], float)
            b, se, G = cluster_se_slope(x, y, x, [d0, d0 ** 2, bs])
            drift_a[(rl, lo, hi)] = (10 * b, 10 * se)
            P("   %-60s G=%2d  %+8.3f kt/decade (se %.3f)  t=%+5.2f  tracks %d" % (
                "index | depth, %s, %d-%d" % (rl, lo, hi), G, 10 * b, 10 * se, b / se, len(rr)))
    # 6c shallow vs deep events
    P("   6c. Which events carry the count? Deep events (minp <= 960) need little help from the gust")
    P("       field; shallow events (minp > 975) exist only because their gust clears 71.7 kt. A")
    P("       gust drift should show first in the shallow ones.")
    for lab, f in (("deep events, minp <= 960", lambda r: float(r["minp"]) <= 960),
                   ("middle events, 960 < minp <= 975", lambda r: 960 < float(r["minp"]) <= 975),
                   ("shallow events, minp > 975", lambda r: float(r["minp"]) > 975)):
        for lo, hi in ((1979, 2000), (2004, 2025)):
            ss = list(range(lo, hi + 1))
            line(P, "%s, %d-%d" % (lab, lo, hi), ols(ss, [sum(1 for r in ev if int(r["season"]) == s and f(r)) for s in ss]))
    for lo, hi in ((1979, 2000), (2004, 2025)):
        ss = list(range(lo, hi + 1))
        sh = [sum(1 for r in ev if int(r["season"]) == s and float(r["minp"]) > 975) for s in ss]
        al = [sum(1 for r in ev if int(r["season"]) == s) for s in ss]
        line(P, "shallow share of events, %d-%d" % (lo, hi), ols(ss, [a / b for a, b in zip(sh, al)]))
    P("")

    # ---------------------------------------------------------------- size of the effect, power
    P("7. HOW LARGE A COUNT BIAS DOES THE DRIFT IMPLY, AND COULD SECTION 4 HAVE SEEN IT?")
    P("   A uniform shift of the index by d kt moves the count above a threshold by about f * d, f")
    P("   the number of tracks per kt at the threshold. A's events are complete above 71.7 kt, so f")
    P("   is measured directly from events with index in [71.7, 73.7).")
    ss = PRE
    yrs = len(ss)
    span = (max(ss) - min(ss)) / 10.0
    for bl, key in (("Atlantic", "atl"), ("Pacific", "pac"), ("both", None)):
        sel = [r for r in ev if int(r["season"]) in ss and (key is None or r["basin"] == key)]
        nf = sum(1 for r in sel if 71.7 <= float(r["gust800_kt"]) < 73.7)
        f = nf / 2.0 / yrs
        mean = len(sel) / yrs
        P("   %-9s pre-2001: %.1f events per season; f = %.2f events per kt per season (%d events in 2 kt)" % (bl, mean, f, nf))
    sel = [r for r in ev if int(r["season"]) in ss]
    f_both = sum(1 for r in sel if 71.7 <= float(r["gust800_kt"]) < 73.7) / 2.0 / yrs
    mean_both = len(sel) / yrs
    P("")
    P("   Implied count trend = f * (index drift per decade), both basins, 1979-2000:")
    cases = [("B moment gust, 955-975 (c9dc994 measure)", 0.65, None),
             ("B track index | depth, 1979-2000 (5b)", drift_b[(1979, 2000)][0], drift_b[(1979, 2000)][1]),
             ("A index | depth, events + nulls, 1979-2000 (6b)", drift_a[("events + nulls", 1979, 2000)][0], drift_a[("events + nulls", 1979, 2000)][1]),
             ("A index | depth, events only, 1979-2000 (6b)", drift_a[("events only", 1979, 2000)][0], drift_a[("events only", 1979, 2000)][1])]
    for lab, d, se in cases:
        P("   %-52s drift %+.2f kt/dec%s -> %+5.1f events/decade, %+5.1f%% of the mean over 1979-2000" % (
            lab, d, "" if se is None else " (se %.2f)" % se, f_both * d, 100 * f_both * d * span / mean_both))
    rr = [r for r in counts if 1979 <= yr(r) <= 2000]
    obs = ols([yr(r) for r in rr], [int(r["atl_era5"]) + int(r["pac_era5"]) for r in rr])
    line(P, "observed: pipeline A events, both basins, 1979-2000", obs)
    P("   Power of that observed test against a true slope equal to the c9dc994 drift's implied")
    P("   slope (two-sided 5%%, se from the data): %.0f%%." % (100 * power(f_both * 0.65, obs["se"], obs["n"] - 2)))
    P("   The minimum slope it detects with 80%% power: %.1f events/decade." % ((obs["crit"] + 0.84) * obs["se"]))
    for key, lab in (("atl_era5", "Atlantic"), ("pac_era5", "Pacific")):
        o = ols([yr(r) for r in rr], [int(r[key]) for r in rr])
        f_b = sum(1 for r in sel if r["basin"] == key[:3] and 71.7 <= float(r["gust800_kt"]) < 73.7) / 2.0 / yrs
        P("   %-8s 1979-1996 test as the session ran it: se %.2f/decade over 18 seasons; implied slope %+.1f;"
          " power %.0f%%" % (lab, ols([yr(r) for r in counts if 1979 <= yr(r) <= 1996],
                                       [int(r[key]) for r in counts if 1979 <= yr(r) <= 1996])["se"],
                            f_b * 0.65, 100 * power(f_b * 0.65, ols([yr(r) for r in counts if 1979 <= yr(r) <= 1996],
                                                                    [int(r[key]) for r in counts if 1979 <= yr(r) <= 1996])["se"], 16)))
    tb72 = (tg >= 71.3) & (tg < 73.3) & pre
    P("   Pipeline B at its 72.3 kt track cut: f = %.2f tracks per kt per season over 1979-2000;"
      " c9dc994's +0.65 kt/decade implies %+.1f events/decade." % (tb72.sum() / 2.0 / yrs, tb72.sum() / 2.0 / yrs * 0.65))
    P("")
    P("8. READING")
    for t in READING.strip("\n").split("\n"):
        P(t)
    out = os.path.join(HERE, "drift_check-result.txt")
    with open(out, "w") as fh:
        fh.write(buf.getvalue())
    print("wrote", out)


READING = """
   Written for the numbers above as committed; if a rerun changes them, rewrite this.

   Finding 1 (archive fixes per event) REPRODUCES: -0.181 /decade, t = -2.54 from the
   payload, -0.180, t = -2.51 from the raw CSVs, -0.217, t = -2.96 for class 'low' only.
   Events and total HF fixes are flat. It weakens from 2006 (t = -2.12; -2.03 from the CSVs,
   p = 0.057). A recording-practice drift on the archive side, as reported.

   Finding 2 (RECORD_START = 2004) REPRODUCES: 2001-02 has 1 archive event against 85 in
   pipeline A; 2002-03 has 22 against 113. The 2001-2006 transfer row is the archive starting
   up. With 2001-03 dropped, A minus archive is flat (-2.2 /decade, t = -0.54).

   Finding 3 (gust at fixed depth) REPRODUCES: +0.645 kt/decade, t = +2.40, 1979-2000, flat
   2004-2025 (t = +0.84). It holds on pipeline A's own index too: A's 800 km gust at fixed
   MSLP 955-975 hPa rises +0.665 kt/decade, t = +2.44, 1979-2000, and is flat after 2004
   (t = +0.61). So the index VALUE is not era-comparable in either pipeline.

   The reconciliation REPRODUCES as computed but DOES NOT SHOW what it was used to show.
   Its evidence that counts do not drift was a per-basin trend over 1979-1996 with t = +0.82
   and -0.02. Against the count bias the drift implies (about +3 per basin per decade), those
   two tests had 15% and 12% power, so a null result there says almost nothing. Over the 22
   pre-archive seasons the same counts rise: both basins +5.1 /decade (t = +1.98, p = 0.06),
   Pacific +5.6 (t = +2.19). The c9dc994 drift, applied to A's measured density at its
   threshold, predicts +5.7 /decade. The observed rise is the size the drift predicts.

   That does not make the rise an artefact. Depth-only counts rise in the same years: B tracks
   deeper than 960 hPa +5.8 /decade (t = +1.96), A deep events +2.7 (t = +1.63). Holding depth
   fixed at event level, the drift estimates are small and uncertain: B track index | depth
   +0.12 (se 0.19); A index | depth +0.21 (se 0.42) with nulls, +0.58 (se 0.46) events only;
   B's depth-standardised excess at 72.3 kt +5.9 /decade (t = +1.95), but the same excess
   runs +5.3 (t = +1.79) over 2004-2025, where the per-moment gust is flat. Shallow events,
   which only the gust can make, show no rise (share t = -0.77). By band, the per-track drift
   sits in 955-965 hPa tracks (+0.98, t = +1.99) and not in shallower ones.

   Answer to decision 1, on these files: the index VALUE drifts before 2001 in both pipelines.
   Whether thresholded COUNTS drift is NOT SETTLED, in either direction. The data bound it:
   from about +1 to about +6 events per decade over 1979-2000 in pipeline A (roughly 2% to
   14% of the mean across the 22 seasons), and they cannot separate that from a real rise in
   deep storms of the same size. Neither c9dc994's "cannot be carried back" nor the walk-back's
   "counts are fine" is supported. For a trend that crosses 2001: a ramp of that size through
   2000 and flat after adds +0.5 to +2.5 /decade to a 1979-2025 OLS trend, against pipeline
   A's own whole-record trend of about +2.4 (Atlantic +1.7, Pacific +0.7). The artefact can be
   anything from a fifth of that trend to all of it.

   What would settle it, not done here: the full sub-threshold track population for A (only
   events and null cases are committed), and an independent wind record before 2001
   (scatterometer from 1991, ship and buoy reports) to test the gust field directly.
"""


def power(effect, se, df):
    c = t_crit(df)
    z = effect / se
    return (1 - t_cdf(c - z, df)) + t_cdf(-c - z, df)


if __name__ == "__main__":
    run()
