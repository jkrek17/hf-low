"""Shared pieces for the hemispheric-pattern channel decomposition. ERA5 PROXY, pipeline A.
Plan: PREREGISTRATION.md (committed before any outcome met the index).

Weekly design as hemispheric (PR 41): 30 weeks of 7 days from 1 October, 22 seasons 2004..2025, by the date of the
track's first fix. Poisson models are fitted with a plain IRLS (no statsmodels) so permutation and bootstrap loops are fast.
"""
import datetime as dt
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
ERA = os.path.join(REPO, "research", "era5")
TRACKS = os.path.join(ERA, "hf_history", "results", "all_tracks.csv.gz")
FIXES = os.path.join(ERA, "intensity", "results", "fixes_2004.csv.gz")
OOS = os.path.join(ERA, "hemispheric", "results", "oos_index_{}.csv")
WEEKLY = os.path.join(ERA, "hemispheric", "results", "weekly_table.csv.gz")
THR = 71.7
DEPTH_CUT = {"atl": 966.2, "pac": 965.0}
NWEEK = 30
SEASONS = list(range(2004, 2026))
BASINS = ["atl", "pac"]
SEED = 20261008
MONTHS = [10, 11, 12, 1, 2, 3, 4]


def week_start(s, k):
    return dt.date(s, 10, 1) + dt.timedelta(days=7 * k)


# ------------------------------------------------------------------ data
def load_tracks(first=2004, last=2025):
    T = pd.read_csv(TRACKS, dtype={"start": str, "peak_time": str})
    T = T[T.basin.isin(BASINS) & (T.season >= first) & (T.season <= last)].copy()
    T["date"] = pd.to_datetime(T.start.str[:8], format="%Y%m%d")
    T["hf"] = (T.gust800_kt >= THR).astype(int)
    T["hfd"] = np.where(T.basin == "atl", T.minp <= DEPTH_CUT["atl"], T.minp <= DEPTH_CUT["pac"]).astype(int)
    T["lon_u"] = np.where((T.basin == "atl") & (T.peak_lon <= 10), T.peak_lon + 360, T.peak_lon)
    return T


def add_origin(T):
    """entrant: the low existed >= 12 h before its first in-domain 00/12 UTC fix. local otherwise (and for the tracks
    missing from the fixes table). mature: first in-domain fix below 1000 hPa."""
    f = pd.read_csv(FIXES, usecols=["track", "time", "age", "msl"])
    g = f.sort_values(["track", "time"]).groupby("track").first()
    T = T.copy()
    age = T.track.map(g.age)
    msl = T.track.map(g.msl)
    T["missing_fix"] = age.isna().astype(int)
    T["entrant"] = (age.fillna(0) >= 12).astype(int)
    T["mature"] = ((msl.fillna(1010) < 1000) & (T.entrant == 0)).astype(int)
    return T


def assign_week(T):
    """season label, week index 0..29 and row index (season-2004)*30+k; drops tracks outside Oct 1 + 210 d."""
    d = T.date
    season = np.where(d.dt.month >= 10, d.dt.year, d.dt.year - 1)
    start = pd.Series(pd.to_datetime(dict(year=season, month=10, day=1)).values, index=d.index)
    day = (d - start).dt.days.values
    T = T.assign(sea=season, day=day)
    T = T[(T.sea >= SEASONS[0]) & (T.sea <= SEASONS[-1]) & (T.day >= 0) & (T.day < 7 * NWEEK)].copy()
    T["wk"] = T.day // 7
    T["row"] = (T.sea - SEASONS[0]) * NWEEK + T.wk
    return T


def prev_week_counts(Tall, mask, basin):
    """Counts in the 7 days before each of the 660 weeks, from the full track list (Sep 24-30 for week 0)."""
    t = Tall[(Tall.basin == basin) & mask]
    daily = t.date.value_counts()
    out = np.zeros(len(SEASONS) * NWEEK)
    for i, s in enumerate(SEASONS):
        for k in range(NWEEK):
            S = pd.Timestamp(week_start(s, k))
            out[i * NWEEK + k] = sum(daily.get(S - pd.Timedelta(days=j), 0) for j in range(1, 8))
    return out


def weekly_counts(T, mask):
    return np.bincount(T.row[mask].values, minlength=len(SEASONS) * NWEEK).astype(float)


def index(basin):
    o = pd.read_csv(OOS.format(basin))
    assert len(o) == len(SEASONS) * NWEEK and (o.season.values == np.repeat(SEASONS, NWEEK)).all()
    x = o.idx.values
    return (x - x.mean()) / x.std(), o.y.values


def month_dummies():
    mid = [(week_start(s, k) + dt.timedelta(days=3)).month for s in SEASONS for k in range(NWEEK)]
    return np.array([[1.0 if m == mm else 0.0 for mm in MONTHS] for m in mid])


def month_of_date(dates):
    return np.array([[1.0 if m == mm else 0.0 for mm in MONTHS] for m in dates.dt.month])


# ------------------------------------------------------------------ models
def pois(X, y, ridge=1e-6, iters=60):
    """Poisson IRLS; tiny ridge keeps sparse month cells finite. Returns beta."""
    b = np.zeros(X.shape[1])
    b[:] = 0.0
    mu0 = max(y.mean(), 1e-3)
    # start at the constant-rate solution through the dummies
    eta = np.full(len(y), np.log(mu0))
    mu = np.exp(eta)
    z = eta
    for _ in range(iters):
        W = mu
        A = X.T @ (X * W[:, None]) + ridge * np.eye(X.shape[1])
        rhs = X.T @ (W * z)
        bn = np.linalg.solve(A, rhs)
        eta = np.clip(X @ bn, -25, 12)
        mu = np.exp(eta)
        z = eta + (y - mu) / np.maximum(mu, 1e-12)
        if np.max(np.abs(bn - b)) < 1e-9:
            b = bn
            break
        b = bn
    return b


def design(mdum, lprev, x):
    return np.column_stack([mdum, lprev, x])


def slope(mdum, lprev, x, y, ridge=1e-6):
    return pois(design(mdum, lprev, x), y, ridge)[-1]


def sandwich(X, y, b, cl):
    mu = np.exp(np.clip(X @ b, -25, 12))
    A = X.T @ (X * mu[:, None]) + 1e-6 * np.eye(X.shape[1])
    Ai = np.linalg.inv(A)
    G = np.zeros((cl.max() + 1, X.shape[1]))
    np.add.at(G, cl, X * (y - mu)[:, None])
    return Ai, G


def joint_se(outs, cl, contrast=None):
    """Cluster-robust SE of a slope (contrast=None) or of a difference of slopes of two outcomes jointly sandwiched,
    exactly the stacked model of PR 14 (its parameters are the separate fits, its covariance is joint)."""
    Ais, Gs = [], []
    for X, y in outs:
        b = pois(X, y)
        Ai, G = sandwich(X, y, b, cl)
        Ais.append(Ai[-1]); Gs.append(G)
    # score of slope parameter for each outcome per cluster: u_k = G_k @ Ai_last_row
    u = np.column_stack([G @ a for G, a in zip(Gs, Ais)])
    V = u.T @ u * len(np.unique(cl)) / (len(np.unique(cl)) - 1)
    c = np.array(contrast if contrast is not None else [1.0] + [0.0] * (len(outs) - 1))
    return float(np.sqrt(c @ V @ c))


def season_perm(rng, ns):
    return rng.permutation(ns)


def permute_idx(x, perm):
    ns = len(perm)
    return x.reshape(ns, NWEEK)[perm].reshape(-1)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q


def resample_rows(rng, seasons_idx):
    """Row indices of a season-block bootstrap draw plus the cluster label of each row."""
    ns = len(seasons_idx)
    pick = rng.integers(0, ns, ns)
    rows = np.concatenate([np.arange(seasons_idx[j] * NWEEK, (seasons_idx[j] + 1) * NWEEK) for j in pick])
    cl = np.repeat(np.arange(ns), NWEEK)
    return rows, cl
