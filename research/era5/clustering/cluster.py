"""Do hurricane-force (HF) lows cluster in time and along the same track?

Plan: PREREGISTRATION.md (committed before any statistic here was computed).
Archive tier is observations of OPC warnings; pipeline A tiers are an ERA5 PROXY
(800 km ocean gust index, HF-equivalent at 71.7 kt; research/era5/hf_history).

usage: python3 cluster.py TIER OUT_DIR [--nsim N] [--nperm N] [--cpc DIR]
  TIER  arch  archive, 2004-05..2025-26 (tier 1)
        r1    pipeline A HF events, 2004-05..2025-26 (tier 2)
        r2    pipeline A HF events, 1979-80..2003-04, within-era (decision 1)
        r3    pipeline A deep cyclones (minp <= 1000 hPa, >= 8 fixes), 1979-80..2025-26,
              era term; pressure based, not gust based

Measures: M1 block dispersion (N1 month+trend, N2 season dummies), M2 inter-arrival,
M3 Knox space-time, M4 local grid, M5 conditioning on lagged NAO/PNA (see plan).
"""
import argparse
import datetime as dt
import json
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
THR = 71.7
NDAY = 210
HRS = NDAY * 24
SEED = 20261008
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
HFDIR = os.path.join(REPO, "research/era5/hf_history/results")
CPC_DEFAULT = "/mnt/project-files/teleconnection-test/cpc_indices"
RE = 6371.0
LS = (1, 2, 3, 7, 14, 30)
MON = np.array([(dt.date(2001, 10, 1) + dt.timedelta(days=k)).month for k in range(NDAY)])
MONIDX = np.array([{10: 0, 11: 1, 12: 2, 1: 3, 2: 4, 3: 5, 4: 6}[m] for m in MON])
TIERS = {
    "arch": dict(s0=2004, s1=2025, src="arch"),
    "r1": dict(s0=2004, s1=2025, src="hf"),
    "r2": dict(s0=1979, s1=2003, src="hf"),
    "r3": dict(s0=1979, s1=2025, src="deep"),
    # calibration check only: Poisson events with the archive's rough seasonal cycle, positions drawn
    # independently of time. No real outcome enters. E should be about 0 and p about uniform.
    "syn": dict(s0=2004, s1=2025, src="syn"),
}


# ------------------------------------------------------------------ geometry
def hav(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def centroid(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    x, y, z = np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)
    x, y, z = x.mean(), y.mean(), z.mean()
    return np.degrees(np.arctan2(z, np.hypot(x, y))), np.degrees(np.arctan2(y, x)) % 360


def hours_since_oct1(code):
    """YYYYMMDDHH ints -> (season, hours since 1 Oct 00Z of that season)."""
    t = pd.to_datetime(pd.Series(np.asarray(code)).astype(str), format="%Y%m%d%H")
    season = (t.dt.year - (t.dt.month < 10)).values
    h = (t - pd.to_datetime(pd.DataFrame({"year": season, "month": 10, "day": 1}))).dt.total_seconds().values / 3600
    return season, h


# ------------------------------------------------------------------ events
def make_event(basin, season, t, fx, extra=None):
    la, lo = centroid(fx[:, 1], fx[:, 2])
    d = dict(basin=basin, season=int(season), t=float(t), lat=la, lon=lo, fx=fx)
    d.update(extra or {})
    return d


def load_arch():
    d = json.load(open(os.path.join(REPO, "docs/data/hf-lows.json")))
    f = {k: i for i, k in enumerate(d["lowFields"])}
    coll = pd.read_csv(os.path.join(REPO, "data/hf_lows/collision_pairs.csv"))
    second = set(zip(coll.basin, coll.other_event.astype(str)))
    rows = []
    for r in d["lows"]:
        if r[f["cls"]] != "low":
            continue
        fx = sorted((x for x in r[f["fixes"]] if x[3] == "HF"), key=lambda x: x[0])
        if not fx:
            continue
        code = np.array([x[0] for x in fx])
        season, h = hours_since_oct1(code)
        s0 = season[0]
        h = (pd.to_datetime(pd.Series(code).astype(str), format="%Y%m%d%H")
             - pd.Timestamp(int(s0), 10, 1)).dt.total_seconds().values / 3600
        arr = np.column_stack([h, [x[1] for x in fx], [x[2] % 360 for x in fx]])
        rid = str(r[f["id"]])
        rows.append(make_event(r[f["basin"]], s0, h[0], arr, dict(
            id=rid, idok=bool(r[f["idOk"]]), second=(r[f["basin"]], rid[:10]) in second)))
    return rows


def load_pipea(kind):
    T = pd.read_csv(os.path.join(HFDIR, "all_tracks.csv.gz"), dtype={"start": str, "peak_time": str})
    T = T[T.basin.isin(["atl", "pac"])]
    if kind == "hf":
        T = T[T.gust800_kt >= THR]
        C = pd.read_csv(os.path.join(HFDIR, "era5_hf_catalog_tracks.csv"))
        C = C[C.track.isin(set(T.track)) & (C.g800 >= THR)]
        grp = {k: v for k, v in C.groupby("track")}
    else:
        T = T[(T.minp <= 1000) & (T.n_fix >= 8)]
        grp = {}
    season, h = hours_since_oct1(T.peak_time.astype(np.int64).values)
    rows = []
    for (tr, bs, la, lo), s, hh in zip(T[["track", "basin", "peak_lat", "peak_lon"]].itertuples(index=False), season, h):
        if kind == "hf" and tr in grp:
            g = grp[tr]
            fh = (pd.to_datetime(pd.Series(g.time.values).astype(str), format="%Y%m%d%H")
                  - pd.Timestamp(int(s), 10, 1)).dt.total_seconds().values / 3600
            arr = np.column_stack([fh, g.lat.values, g.lon.values % 360])
        else:
            arr = np.array([[hh, la, lo % 360]])
        rows.append(make_event(bs, s, hh, arr, dict(id=str(tr))))
    return rows


def window(rows, basin, s0, s1):
    """Events of a basin inside Oct 1..Apr 28 of seasons s0..s1. Returns (DataFrame-like list, n dropped)."""
    allb = [e for e in rows if e["basin"] == basin and s0 <= e["season"] <= s1]
    keep = [e for e in allb if 0 <= e["t"] < HRS]
    return keep, len(allb) - len(keep)


def daily(ev, s0, s1):
    S = s1 - s0 + 1
    y = np.zeros((S, NDAY), int)
    for e in ev:
        y[e["season"] - s0, int(e["t"] // 24)] += 1
    return y


# ------------------------------------------------------------------ indices
def read_cpc(path):
    vals = {}
    for ln in open(path):
        if len(ln) < 11 or not ln[:4].isdigit():
            continue
        v = float(ln[10:])
        vals[pd.Timestamp(int(ln[:4]), int(ln[5:7]), int(ln[8:10]))] = np.nan if v <= -98 else v
    return pd.Series(vals).sort_index().asfreq("D")


ONI_SEAS = ["DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ"]


def read_oni(cpc):
    oni = {}
    t = json.load(open(os.path.join(REPO, "docs/data/teleconnections.json")))["oni"]
    y0, m0 = map(int, t["start"].split("-"))
    for i, v in enumerate(t["values"]):
        if v is not None:
            oni[(y0 + (m0 - 1 + i) // 12, (m0 - 1 + i) % 12 + 1)] = v
    for fn in sorted(os.listdir(cpc)):
        if fn.startswith("oni"):
            for ln in open(os.path.join(cpc, fn)):
                p = ln.split()
                if len(p) == 4 and p[0] in ONI_SEAS:
                    oni.setdefault((int(p[1]), ONI_SEAS.index(p[0]) + 1), float(p[3]))
    return oni


def index_tables(cpc, s0, s1):
    """(S, NDAY) tables of the lagged (days d-10..d-4) NAO and PNA, and ONI one month earlier."""
    out = {}
    for k in ("NAO", "PNA"):
        s = read_cpc(os.path.join(cpc, f"norm.daily.{k.lower()}.index.b500101.current.ascii"))
        lag = s.rolling(7, min_periods=7).mean().shift(4)
        A = np.array([[lag.get(pd.Timestamp(s_, 10, 1) + pd.Timedelta(days=d), np.nan) for d in range(NDAY)]
                      for s_ in range(s0, s1 + 1)])
        out[k] = A
    oni = read_oni(cpc)
    O = np.array([[oni.get(((pd.Timestamp(s_, 10, 1) + pd.Timedelta(days=d - 30)).year,
                            (pd.Timestamp(s_, 10, 1) + pd.Timedelta(days=d - 30)).month), np.nan)
                   for d in range(NDAY)] for s_ in range(s0, s1 + 1)])
    out["ONI"] = O
    nmiss = {k: int(np.isnan(v).sum()) for k, v in out.items()}
    for k, v in out.items():
        v = np.where(np.isnan(v), np.nanmean(v), v)
        out[k] = (v - v.mean()) / v.std()
    return out, nmiss


# ------------------------------------------------------------------ Poisson GLM and dispersion
def make_X(kind, S, tr, era=None, covs=None, nd=NDAY, harm=False):
    mon = MONIDX[:nd]
    if harm:
        k = np.arange(nd) / NDAY * 2 * np.pi
        M = np.column_stack([np.ones(nd), np.cos(k), np.sin(k), np.cos(2 * k), np.sin(2 * k)])
    else:
        M = np.eye(7)[mon]
    Mt = np.tile(M, (S, 1))
    trr = np.repeat(tr, nd)[:, None]
    if kind == "N2":
        Mt = Mt[:, 1:] if not harm else Mt[:, 1:]
        sd = np.kron(np.eye(S), np.ones((nd, 1)))
        return np.hstack([Mt, sd])
    cols = [Mt, trr]
    if era is not None:
        cols.append(np.repeat(era, nd)[:, None])
    for c in (covs or []):
        cols.append(np.asarray(c)[:, :nd].reshape(-1, 1))
    return np.hstack(cols)


def fit_pois(y, X, iters=40):
    mu = y + 0.1 if y.sum() > 0 else np.full_like(y, 0.1, float)
    eta = np.log(mu)
    b = np.zeros(X.shape[1])
    dev0 = np.inf
    for _ in range(iters):
        z = eta + (y - mu) / mu
        XtW = X.T * mu
        b = np.linalg.solve(XtW @ X + 1e-9 * np.eye(X.shape[1]), XtW @ z)
        eta = np.clip(X @ b, -30, 30)
        mu = np.exp(eta)
        dev = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0) - (y - mu))
        if abs(dev0 - dev) < 1e-9 * (abs(dev) + 1):
            break
        dev0 = dev
    return b, mu


def phi_blocks(y, mu, L, p, off=0):
    S, nd = y.shape
    n = (nd - off) // L * L
    Y = y[:, off:off + n].reshape(S, -1, L).sum(2)
    M = mu[:, off:off + n].reshape(S, -1, L).sum(2)
    return float(((Y - M) ** 2 / M).sum() / (Y.size - p))


def phis(y, mu, p, Ls, off=0):
    return np.array([phi_blocks(y, mu, L, p, off) for L in Ls])


def dispersion(y, X, Ls, rng, nsim, off=0):
    S, nd = y.shape
    b, mu = fit_pois(y.ravel().astype(float), X)
    mu2 = mu.reshape(S, nd)
    p = X.shape[1]
    obs = phis(y, mu2, p, Ls, off)
    sims = np.empty((nsim, len(Ls)))
    for i in range(nsim):
        ys = rng.poisson(mu2)
        _, ms = fit_pois(ys.ravel().astype(float), X)
        sims[i] = phis(ys, ms.reshape(S, nd), p, Ls, off)
    null_mean = sims.mean(0)
    pv = (1 + (sims >= obs).sum(0)) / (1 + nsim)
    return dict(phi=obs, null=null_mean, E=obs - null_mean, p=pv, mu=mu2, p_par=p,
                null_sd=sims.std(0), b=b)


def boot_E(y, kinds, tr, era, covs_by_kind, Ls, null_by_kind, rng, nboot, off=0):
    """Season-block bootstrap of E for each model kind (percentile CIs); returns dict kind -> (nboot, len(Ls))."""
    S, nd = y.shape
    out = {k: np.empty((nboot, len(Ls))) for k in kinds}
    for b in range(nboot):
        ix = rng.integers(S, size=S)
        yb = y[ix]
        for k in kinds:
            cov = [c[ix] for c in covs_by_kind.get(k, [])]
            X = make_X("N2" if k == "N2" else "N1", S, tr[ix], None if era is None else era[ix], cov, nd)
            _, mu = fit_pois(yb.ravel().astype(float), X)
            out[k][b] = phis(yb, mu.reshape(S, nd), X.shape[1], Ls, off) - null_by_kind[k]
    return out


def ci(a, q=(2.5, 97.5)):
    a = a[~np.isnan(a)]
    return tuple(np.percentile(a, q)) if len(a) else (np.nan, np.nan)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    run = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        run = min(run, p[i] * n / rank)
        q[i] = run
    return q


# ------------------------------------------------------------------ M2 inter-arrival
def gap_stats(season, t, lam_cum):
    """season (n,), t hours (n,), lam_cum: function (season, t) -> cumulative intensity. Returns dict of stats."""
    order = np.lexsort((t, season))
    s, tt = season[order], t[order]
    same = s[1:] == s[:-1]
    gap = (tt[1:] - tt[:-1])[same]
    L = lam_cum(s, tt)
    rg = (L[1:] - L[:-1])[same]
    n = len(tt)
    G48 = int((gap <= 48).sum())
    S72 = float((gap <= 72).sum() / max(n, 1))
    cv = float(rg.std(ddof=1) / rg.mean()) if len(rg) > 2 and rg.mean() > 0 else np.nan
    xs = np.sort(rg)
    F = 1 - np.exp(-xs)
    m = len(xs)
    ks = float(max(np.max(np.arange(1, m + 1) / m - F), np.max(F - np.arange(0, m) / m))) if m else np.nan
    return dict(n=n, G48=G48, S72=S72, CV=cv, KS=ks, gaps=len(gap))


def make_lam(mu, s0):
    cum = np.cumsum(mu, axis=1)
    prev = np.hstack([np.zeros((mu.shape[0], 1)), cum[:, :-1]])

    def lam(season, t):
        si = (season - s0).astype(int)
        d = np.clip((t // 24).astype(int), 0, NDAY - 1)
        return prev[si, d] + mu[si, d] * ((t - 24 * d) / 24)
    return lam


def m2(ev, mu, s0, rng, nsim):
    season = np.array([e["season"] for e in ev])
    t = np.array([e["t"] for e in ev])
    lam = make_lam(mu, s0)
    obs = gap_stats(season, t, lam)
    sims = {k: np.empty(nsim) for k in ("n", "G48", "S72", "CV", "KS")}
    S = mu.shape[0]
    sidx = np.repeat(np.arange(S), NDAY)
    didx = np.tile(np.arange(NDAY), S)
    for i in range(nsim):
        c = rng.poisson(mu).ravel()
        si = np.repeat(sidx, c)
        di = np.repeat(didx, c)
        tt = di * 24 + 6 * rng.integers(0, 4, size=len(di))
        g = gap_stats(si + s0, tt.astype(float), lam)
        for k in sims:
            sims[k][i] = g[k]
    res = dict(obs=obs, null_mean={k: float(np.nanmean(v)) for k, v in sims.items()})
    res["p"] = {k: float((1 + np.nansum(sims[k] >= obs[k])) / (1 + nsim)) for k in ("G48", "S72", "CV", "KS")}
    res["ratio_G48"] = obs["G48"] / res["null_mean"]["G48"] if res["null_mean"]["G48"] else np.nan
    return res


# ------------------------------------------------------------------ M3 Knox
TAUS = (48, 72, 120)
DS = (500, 1000, 1500, 2000)


def m3(ev, rng, nperm):
    n = len(ev)
    season = np.array([e["season"] for e in ev])
    t = np.array([e["t"] for e in ev])
    lat = np.array([e["lat"] for e in ev])
    lon = np.array([e["lon"] for e in ev])
    ii, jj, dd = [], [], []
    for s in np.unique(season):
        ix = np.where(season == s)[0]
        a_, b_ = np.triu_indices(len(ix), 1)
        D = hav(lat[ix[a_]], lon[ix[a_]], lat[ix[b_]], lon[ix[b_]])
        k = D <= max(DS)
        ii.append(ix[a_][k]); jj.append(ix[b_][k]); dd.append(D[k])
    ii, jj, D = np.concatenate(ii), np.concatenate(jj), np.concatenate(dd)
    dbin = np.searchsorted(DS, D, side="left")  # 0: <=500, 1: <=1000, ...
    mon = MONIDX[np.clip((t // 24).astype(int), 0, NDAY - 1)]
    grp = (season - season.min()) * 7 + mon
    order = np.argsort(grp, kind="stable")
    t_sorted = t[order]

    def grid(tt):
        dt_ = np.abs(tt[ii] - tt[jj])
        return np.array([np.bincount(dbin[dt_ <= tau], minlength=len(DS)).cumsum() for tau in TAUS])

    obs = grid(t)
    perm = np.empty((nperm,) + obs.shape)
    for b in range(nperm):
        idx = np.lexsort((rng.random(n), grp))
        tp = np.empty(n)
        tp[idx] = t_sorted
        perm[b] = grid(tp)
    mean = perm.mean(0)
    p = (1 + (perm >= obs).sum(0)) / (1 + nperm)
    return dict(obs=obs, null_mean=mean, ratio=obs / np.where(mean > 0, mean, np.nan), p=p, pairs=int(len(ii)),
                pair_idx=(ii, jj, D))


# ------------------------------------------------------------------ M4 local grid
def cells_for(basin):
    lats = np.arange(30, 70, 5.0)
    lons = np.arange(-90, 20, 10.0) if basin == "atl" else np.arange(140, 240, 10.0)
    return [(a, o % 360) for a in lats for o in lons]


def m4(ev, y_base_tr, s0, S, rng, nsim, nperm, basin, tr, era, radius=700.0, minpass=40):
    cells = cells_for(basin)
    n = len(ev)
    fl = np.concatenate([e["fx"][:, 1] for e in ev])
    fo = np.concatenate([e["fx"][:, 2] for e in ev])
    fh = np.concatenate([e["fx"][:, 0] for e in ev])
    ei = np.concatenate([np.full(len(e["fx"]), k) for k, e in enumerate(ev)])
    PASS = np.zeros((n, len(cells)), bool)
    LOC = np.full((n, len(cells)), np.inf)
    for c, (a, o) in enumerate(cells):
        inside = hav(fl, fo, a, o) <= radius
        for k in np.unique(ei[inside]):
            m = inside & (ei == k)
            PASS[k, c] = True
            LOC[k, c] = fh[m].min()
    t = np.array([e["t"] for e in ev])
    season = np.array([e["season"] for e in ev])
    OFF = np.where(PASS, LOC - t[:, None], 0.0)
    sidx = season - s0
    # passages inside the window only
    ld = np.floor((t[:, None] + OFF) / 24)
    inwin = PASS & (ld >= 0) & (ld < NDAY)
    npass = inwin.sum(0)
    use = np.where(npass >= minpass)[0]
    X = make_X("N1", S, tr, era)
    rows, cols = np.nonzero(inwin)
    Y = np.zeros((len(cells), S, NDAY), int)
    np.add.at(Y, (cols, sidx[rows], ld[rows, cols].astype(int)), 1)
    out = []
    mus = {}
    for c in use:
        yc = Y[c].astype(float)
        b, mu = fit_pois(yc.ravel(), X)
        mu2 = mu.reshape(S, NDAY)
        mus[c] = mu2
        p = X.shape[1]
        phi = phi_blocks(Y[c], mu2, 7, p)
        sims = np.empty(nsim)
        for i in range(nsim):
            sims[i] = phi_blocks(rng.poisson(mu2), mu2, 7, p)
        out.append(dict(cell=int(c), lat=cells[c][0], lon=cells[c][1], passages=int(npass[c]),
                        mean_week=float(mu2.sum() / (S * NDAY / 7)), phi=phi, E=phi - sims.mean(),
                        p_par=float((1 + (sims >= phi).sum()) / (1 + nsim))))
    df = pd.DataFrame(out)
    if df.empty:
        return df, dict(field_obs=0, field_p=np.nan)
    df["q_par"] = bh(df.p_par.values)
    # position-shuffle null
    mon = MONIDX[np.clip((t // 24).astype(int), 0, NDAY - 1)]
    grp = (season - season.min()) * 7 + mon
    p_par_cols = X.shape[1]
    MU = np.stack([mus[c] for c in use])  # (C, S, NDAY)
    MUw = MU.reshape(len(use), S, NDAY // 7, 7).sum(-1)
    PHI = np.empty((nperm, len(use)))
    for b in range(nperm):
        donor_order = np.lexsort((rng.random(n), grp))
        donor = np.empty(n, int)
        donor[np.argsort(grp, kind="stable")] = donor_order
        r2, c2 = np.nonzero(PASS[donor][:, use])
        hrs = t[r2] + OFF[donor[r2], use[c2]]
        dd = np.floor(hrs / 24).astype(int)
        ok = (dd >= 0) & (dd < NDAY)
        flat = (c2[ok] * S + sidx[r2[ok]]) * NDAY + dd[ok]
        Yb = np.bincount(flat, minlength=len(use) * S * NDAY).reshape(len(use), S, NDAY // 7, 7).sum(-1)
        PHI[b] = ((Yb - MUw) ** 2 / MUw).sum((1, 2)) / (S * (NDAY // 7) - p_par_cols)
    phi_obs = df.phi.values
    df["p_shuf"] = (1 + (PHI >= phi_obs).sum(0)) / (1 + nperm)
    df["E_shuf"] = phi_obs - PHI.mean(0)
    df["q_shuf"] = bh(df.p_shuf.values)
    # field significance by rank-p within the permutation set
    rank_p = np.empty_like(PHI)
    for c in range(PHI.shape[1]):
        srt = np.sort(PHI[:, c])
        rank_p[:, c] = (nperm - np.searchsorted(srt, PHI[:, c], side="left")) / nperm
    Nb = (rank_p <= 0.05).sum(1)
    Nobs = int((df.p_shuf.values <= 0.05).sum())
    meta = dict(field_obs=Nobs, field_null_mean=float(Nb.mean()), field_null_p95=float(np.percentile(Nb, 95)),
                field_p=float((1 + (Nb >= Nobs).sum()) / (1 + nperm)), cells_tested=int(len(use)),
                cells_total=len(cells))
    return df, meta


# ------------------------------------------------------------------ tier run
def synthetic_rows(seed=7):
    rng = np.random.default_rng(seed)
    pool = {b: [e for e in load_arch() if e["basin"] == b] for b in ("atl", "pac")}
    w = np.array([.08, .14, .18, .20, .17, .14, .09])
    day_w = w[MONIDX] / np.bincount(MONIDX)[MONIDX]
    rows = []
    for b in ("atl", "pac"):
        for s in range(2004, 2026):
            n = rng.poisson(40)
            days = rng.choice(NDAY, size=n, p=day_w / day_w.sum())
            for d in days:
                t = d * 24 + 6 * rng.integers(0, 4)
                src = pool[b][rng.integers(len(pool[b]))]
                fx = src["fx"].copy()
                fx[:, 0] = fx[:, 0] - fx[0, 0] + t
                rows.append(make_event(b, s, t, fx))
    return rows


def tier_data(name):
    cfg = TIERS[name]
    if cfg["src"] == "syn":
        rows = synthetic_rows()
    elif cfg["src"] == "arch":
        rows = load_arch()
    else:
        rows = load_pipea(cfg["src"])
    return rows, cfg


def run(name, out, nsim, nperm, cpc, nboot):
    rng = np.random.default_rng(SEED + sum(map(ord, name)))
    rows, cfg = tier_data(name)
    s0, s1 = cfg["s0"], cfg["s1"]
    S = s1 - s0 + 1
    tr = np.arange(S) - (S - 1) / 2
    era = (np.arange(s0, s1 + 1) >= 2004).astype(float) if name == "r3" else None
    idx, nmiss = index_tables(cpc, s0, s1)
    os.makedirs(out, exist_ok=True)
    res = {"tier": name, "seasons": [s0, s1], "nmiss_index_days": nmiss, "basins": {}}
    tests = []
    for basin in ("atl", "pac"):
        ev, dropped = window(rows, basin, s0, s1)
        y = daily(ev, s0, s1)
        ib = "NAO" if basin == "atl" else "PNA"
        R = dict(n_events=len(ev), n_outside_window=dropped, per_season_mean=float(y.sum() / S))
        # ---- M1 + M5 dispersion models
        covs = {"N1": [], "N2": [], "IDX": [idx[ib]], "IDX2": [idx["NAO"], idx["PNA"], idx["ONI"]]}
        sm_ = idx[ib].mean(1, keepdims=True) * np.ones((1, NDAY))
        sm_ = (sm_ - sm_.mean()) / sm_.std()
        an = idx[ib] - idx[ib].mean(1, keepdims=True)
        an = an / an.std()
        covs["SPLIT"] = [sm_, an]
        D = {}
        for k in ("N1", "N2", "IDX", "IDX2", "SPLIT"):
            X = make_X("N2" if k == "N2" else "N1", S, tr, era, covs[k])
            D[k] = dispersion(y, X, LS, rng, nsim)
        R["M1"] = {k: dict(phi=v["phi"].tolist(), null=v["null"].tolist(), E=v["E"].tolist(), p=v["p"].tolist())
                   for k, v in D.items()}
        # bootstrap CIs
        nullk = {k: D[k]["null"] for k in ("N1", "N2", "IDX")}
        B = boot_E(y, ("N1", "N2", "IDX"), tr, era, {"IDX": [idx[ib]]}, (7, 30), {k: nullk[k][[3, 5]] for k in nullk}, rng, nboot)
        R["M1_ci"] = {k: {"W7": ci(B[k][:, 0]), "W30": ci(B[k][:, 1])} for k in B}
        # share explained by index
        share = {}
        for j, w in enumerate(("W7", "W30")):
            e0, ei_ = B["N1"][:, j], B["IDX"][:, j]
            ok = e0 > 0.02
            sh = np.where(ok, 1 - ei_ / np.where(ok, e0, 1), np.nan)
            share[w] = dict(point=float(1 - D["IDX"]["E"][[3, 5][j]] / D["N1"]["E"][[3, 5][j]]),
                            ci=ci(sh), valid_frac=float(ok.mean()),
                            dE=float(D["N1"]["E"][[3, 5][j]] - D["IDX"]["E"][[3, 5][j]]),
                            dE_ci=ci(e0 - ei_))
        R["share_explained"] = share
        # index coefficient (M-idx): season-block permutation and cluster Wald
        X = make_X("N1", S, tr, era, [idx[ib]])
        b0, mu0 = fit_pois(y.ravel().astype(float), X)
        res_ = (y.ravel() - mu0).reshape(S, NDAY)
        XS = X.reshape(S, NDAY, -1)
        A = X.T @ (X * mu0[:, None])
        g = np.einsum("snp,sn->sp", XS, res_)
        V = np.linalg.inv(A) @ (g.T @ g) @ np.linalg.inv(A)
        coef = float(b0[-1])
        nperm_i = 2000
        cp = np.empty(nperm_i)
        for i in range(nperm_i):
            Xp = make_X("N1", S, tr, era, [idx[ib][rng.permutation(S)]])
            cp[i] = fit_pois(y.ravel().astype(float), Xp)[0][-1]
        R["index_coef"] = dict(index=ib, logRR_per_SD=coef, RR=float(np.exp(coef)),
                               se_cluster=float(np.sqrt(V[-1, -1])),
                               p_perm=float((1 + (np.abs(cp) >= abs(coef)).sum()) / (1 + nperm_i)))
        # dispersion curve for N1
        # ---- M2
        R["M2"] = m2(ev, D["N1"]["mu"], s0, rng, nsim)
        # ---- M3
        K = m3(ev, rng, nperm)
        R["M3"] = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in K.items() if k != "pair_idx"}
        pi, pj, pd_ = K["pair_idx"]
        # ---- M4
        if name in ("arch", "r1", "syn"):
            cdf, meta = m4(ev, None, s0, S, rng, 500, min(nperm, 1000), basin, tr, era)
            cdf.to_csv(os.path.join(out, f"{name}_{basin}_cells.csv"), index=False)
            R["M4"] = meta
        # ---- tests (tier family)
        tests += [
            (basin, "M1 N1 W7", float(D["N1"]["p"][3])),
            (basin, "M1 N1 W30", float(D["N1"]["p"][5])),
            (basin, "M2 gaps<=48h", float(R["M2"]["p"]["G48"])),
            (basin, "M3 Knox 3d/1000km", float(K["p"][1, 1])),
            (basin, "index coef (perm)", float(R["index_coef"]["p_perm"])),
        ]
        # event/pair table for maps
        pd.DataFrame(dict(season=[e["season"] for e in ev], t=[e["t"] for e in ev],
                          lat=[e["lat"] for e in ev], lon=[e["lon"] for e in ev])).to_csv(
            os.path.join(out, f"{name}_{basin}_events.csv"), index=False)
        res["basins"][basin] = R
        np.save(os.path.join(out, f"{name}_{basin}_daily.npy"), y)
    T = pd.DataFrame(tests, columns=["basin", "test", "p"])
    T["q_bh"] = bh(T.p.values)
    T.to_csv(os.path.join(out, f"{name}_tests.csv"), index=False)
    res["tests"] = T.to_dict("records")
    json.dump(res, open(os.path.join(out, f"{name}_results.json"), "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("tier")
    ap.add_argument("out")
    ap.add_argument("--nsim", type=int, default=2000)
    ap.add_argument("--nperm", type=int, default=4000)
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--cpc", default=CPC_DEFAULT)
    a = ap.parse_args()
    run(a.tier, a.out, a.nsim, a.nperm, a.cpc, a.nboot)
