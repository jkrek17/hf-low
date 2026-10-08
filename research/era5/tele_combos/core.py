"""Shared pieces for the teleconnection-combination tests (ERA5 PROXY, pipeline A).

Plan: PREREGISTRATION.md (committed before any fit). Catalog: CATALOG.md.

Units of analysis
  location test : one pipeline A track (HF-equivalent: gust800_kt >= 71.7, or any cyclone), placed on its
                  genesis day (first fix). Outcome = (peak_lon, peak_lat), where the 800 km gust index peaks.
                  Atlantic longitude is unwrapped (lon < 180 -> +360).
  count test    : daily genesis counts per basin, Poisson.
Predictors are standardised over the analysis DAYS (not events); effects are per SD, interaction per SD x SD.

Inference for the interaction coefficient gamma (product of the two standardised predictors):
  p     : sign-flip score test (restricted wild cluster bootstrap with Rademacher weights, seasons as clusters,
          efficient score with nuisance terms projected out); 2 df for location, 1 df for counts.
  CI    : season-pairs bootstrap, percentile; also CR1 clustered SE for reference.
"""
import datetime as dt
import importlib.util
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CPC = os.environ.get("CPC_DIR", "/mnt/project-files/teleconnection-test/cpc_indices")
TRACKS = os.path.join(REPO, "research/era5/hf_history/results/all_tracks.csv.gz")
SPV_CSV = os.path.join(HERE, "results", "spv_10hpa_60n_daily.csv")
SEED = 20261008

_spec = importlib.util.spec_from_file_location("freq_split_mod", os.path.join(REPO, "research/era5/freq_split/split.py"))
fs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fs)

THR = 71.7
BOXES = {  # name: (lat0, lat1, lon0, lon1); Atlantic lon on the unwrapped 262..370 scale
    "pac": {"P-W": (27, 45, 135, 165), "P-C": (35, 55, 165, 200), "P-E": (40, 60, 200, 240)},
    "atl": {"A-W": (27, 50, 262, 305), "A-C": (40, 62, 305, 335), "A-E": (50, 67, 335, 370)},
}


# ----------------------------------------------------------------- predictors
def mjo_daily():
    """Pentad MJO longitude indices -> daily DataFrame (nearest pentad centre), sign flipped so that
    POSITIVE = ENHANCED convection (CPC files are negative = enhanced)."""
    lines = [ln for ln in open(os.path.join(CPC, "proj_norm_order.ascii")) if ln.strip()]
    cols = lines[1].split()[1:]
    rows = []
    for ln in lines[2:]:
        p = ln.split()
        if any("*" in v for v in p[1:]):
            continue
        rows.append([pd.Timestamp(dt.datetime.strptime(p[0], "%Y%m%d"))] + [-float(v) for v in p[1:]])
    P = pd.DataFrame(rows, columns=["date"] + cols).set_index("date")
    P = P / P.loc["1979-01-01":"2025-12-31"].std()          # unit SD per longitude, 1979-2025
    days = pd.date_range(P.index[0], P.index[-1], freq="D")
    nearest = P.index.get_indexer(days, method="nearest")
    D = pd.DataFrame(P.values[nearest], index=days, columns=cols)
    return D


def highpass(s, n=90):
    """Causal high-pass: value minus the mean of the trailing n days (including today)."""
    return s - s.rolling(n, min_periods=n).mean()


def spv_daily():
    if not os.path.exists(SPV_CSV):
        return None
    s = pd.read_csv(SPV_CSV, parse_dates=["date"]).set_index("date")["u"].asfreq("D")
    return s


def predictors():
    """Daily DataFrame of every predictor used, indexed by genesis date.

    *_lag  : the window the plan fixes (days -10..-4 for NAO/PNA; day -7 for ONI month and MJO pentad;
             days -30..-11 for the polar vortex)
    *_now  : secondary, same-time (days -3..+3 for NAO/PNA; day 0 for ONI month and MJO; days -10..+10 SPV)
    """
    out = {}
    for k in ("NAO", "PNA"):
        s = fs.read_cpc(os.path.join(CPC, f"norm.daily.{k.lower()}.index.b500101.current.ascii"))
        out[k + "_lag"] = s.rolling(7, min_periods=7).mean().shift(4)
        out[k + "_now"] = s.rolling(7, min_periods=7, center=True).mean()
    oni = fs.read_oni(CPC, REPO)
    days = out["NAO_lag"].index
    out["ONI_lag"] = pd.Series([oni.get((d.year, d.month), np.nan) for d in days - pd.Timedelta(days=7)], index=days)
    out["ONI_now"] = pd.Series([oni.get((d.year, d.month), np.nan) for d in days], index=days)
    M = mjo_daily()
    wp = highpass(M[["120E", "140E"]].mean(axis=1))        # west Pacific / Maritime Continent
    dl = highpass(M[["160E", "120W"]].mean(axis=1))        # dateline / central-east Pacific
    for nm, s in (("MJOWP", wp), ("MJODL", dl)):
        s = s.reindex(days)
        out[nm + "_lag"] = s.shift(7)
        out[nm + "_now"] = s
    spv = spv_daily()
    if spv is not None:
        doy = spv.groupby(spv.index.dayofyear)
        clim = doy.transform("mean")
        sd = doy.transform("std")
        a = ((spv - clim) / sd).reindex(days)
        out["SPV_lag"] = a.rolling(20, min_periods=20).mean().shift(11)       # days -30..-11
        out["SPV_now"] = a.rolling(21, min_periods=21, center=True).mean()    # days -10..+10
    return pd.DataFrame(out)


# ----------------------------------------------------------------- data
def load_tracks():
    T = fs.load_tracks(TRACKS)
    T["lon"] = np.where(T.peak_lon < 180, T.peak_lon + 360, T.peak_lon) if True else T.peak_lon
    # unwrap only the Atlantic; the Pacific domain (135-240E) never crosses 0
    T["lon"] = np.where(T.basin == "atl", np.where(T.peak_lon < 180, T.peak_lon + 360, T.peak_lon), T.peak_lon)
    T["lat"] = T.peak_lat
    return T


def build(T, I, basin, A, B, s0, s1, sample, now=False, era=False, window="octapr", zmask=None):
    """Return (D, E): D = analysis days with standardised predictors; E = events with the day's predictors.
    sample: 'hf' (gust >= 71.7) or 'all'."""
    suf = "_now" if now else "_lag"
    D = fs.window_days(window, s0, s1)
    raw = I.reindex(D.date)[[A + suf, B + suf]].reset_index(drop=True)
    ok = raw.notna().all(axis=1).values
    D = D[ok].reset_index(drop=True)
    raw = raw[ok].reset_index(drop=True)
    z = (raw - raw.mean()) / raw.std()
    D["A"], D["B"] = z.iloc[:, 0].values, z.iloc[:, 1].values
    D["AB"] = D.A * D.B
    Tw = fs.assign(T, window, s0, s1)
    Tw = Tw[(Tw.basin == basin)]
    if sample == "hf":
        Tw = Tw[Tw.hf == 1]
    E = Tw.merge(D[["season", "day", "A", "B", "AB"]], on=["season", "day"], how="inner")
    E["month"] = E.gen.dt.month
    D["era"] = (D.season <= 2000).astype(float) if era else 0.0
    E["era"] = (E.season <= 2000).astype(float) if era else 0.0
    return D, E


def covariates(df, era):
    X = pd.DataFrame(index=df.index)
    for mo in sorted(df.month.unique())[1:]:
        X[f"m{mo}"] = (df.month.values == mo).astype(float)
    X["trend"] = (df.season.values - 2014.5) / 10.0
    if era:
        X["era"] = df.era.values
    X["const"] = 1.0
    return X


# ----------------------------------------------------------------- OLS location test
def _ols_resid(X, Y):
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return b, Y - X @ b


def location_test(E, era=False, nflip=50000, nboot=2000, seed=SEED):
    """Interaction on (lon, lat). Returns dict with gamma (2), cluster SE, 95% season-bootstrap CI,
    Wald-type 2 df sign-flip p, and the main effects."""
    rng = np.random.default_rng(seed)
    E = E.reset_index(drop=True)
    cov = covariates(E, era)
    Y = E[["lon", "lat"]].values.astype(float)
    g = E.season.values
    seasons = np.unique(g)
    XF = np.column_stack([E.A, E.B, E.AB, cov.values])
    XR = np.column_stack([E.A, E.B, cov.values])
    b, res = _ols_resid(XF, Y)
    # restricted efficient score
    xt = E.AB.values - XR @ np.linalg.lstsq(XR, E.AB.values, rcond=None)[0]
    _, er = _ols_resid(XR, Y)
    U = np.array([(xt[g == s, None] * er[g == s]).sum(axis=0) for s in seasons])     # G x 2
    V = U.T @ U
    Vi = np.linalg.pinv(V)
    q = U.sum(0) @ Vi @ U.sum(0)
    W = rng.choice([-1.0, 1.0], size=(nflip, len(seasons)))
    Us = W @ U
    qs = np.einsum("bi,ij,bj->b", Us, Vi, Us)
    p2 = (1 + (qs >= q - 1e-12).sum()) / (nflip + 1)
    # per-coordinate sign-flip p
    p1 = []
    for k in range(2):
        t = U[:, k].sum() / np.sqrt((U[:, k] ** 2).sum())
        ts = (W @ U[:, k]) / np.sqrt((U[:, k] ** 2).sum())
        p1.append((1 + (np.abs(ts) >= abs(t) - 1e-12).sum()) / (nflip + 1))
    # cluster SE (CR1) of full model
    XtXi = np.linalg.inv(XF.T @ XF)
    Gn = len(seasons)
    se = []
    for k in range(2):
        S = np.array([(XF[g == s] * res[g == s, k:k + 1]).sum(0) for s in seasons])
        Vc = XtXi @ (S.T @ S) @ XtXi * Gn / (Gn - 1) * (len(E) - 1) / (len(E) - XF.shape[1])
        se.append(np.sqrt(np.diag(Vc)))
    se = np.array(se)           # 2 x p
    # season-pairs bootstrap
    idx_by = {s: np.where(g == s)[0] for s in seasons}
    boots = np.empty((nboot, 2, 3))
    for r in range(nboot):
        pick = rng.choice(seasons, size=len(seasons), replace=True)
        ii = np.concatenate([idx_by[s] for s in pick])
        try:
            bb, _ = _ols_resid(XF[ii], Y[ii])
            boots[r] = bb[:3].T
        except Exception:
            boots[r] = np.nan
    lo, hi = np.nanpercentile(boots, [2.5, 97.5], axis=0)       # 2 x 3
    return dict(n_events=len(E), n_seasons=len(seasons),
                A=b[0], B=b[1], gamma=b[2],
                se_A=se[:, 0], se_B=se[:, 1], se_gamma=se[:, 2],
                ci_gamma=np.column_stack([lo[:, 2], hi[:, 2]]),
                p_joint=p2, p_lon=p1[0], p_lat=p1[1],
                y_sd=Y.std(0), coef=b)


# ----------------------------------------------------------------- Poisson count test
def newton(y, X, iters=40, off=None):
    b = np.zeros(X.shape[1])
    b[-1] = np.log(max(y.mean(), 1e-9))
    for _ in range(iters):
        mu = np.exp(X @ b)
        step = np.linalg.solve(X.T @ (X * mu[:, None]) + 1e-10 * np.eye(len(b)), X.T @ (y - mu))
        b += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b


def daily_counts(D, E):
    key = pd.MultiIndex.from_frame(D[["season", "day"]])
    c = E.groupby(["season", "day"]).size().reindex(key, fill_value=0).values
    return c.astype(float)


def count_test(D, y, era=False, nflip=50000, nboot=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    D = D.reset_index(drop=True)
    D["month"] = D.date.dt.month
    cov = covariates(D, era)
    g = D.season.values
    seasons = np.unique(g)
    XF = np.column_stack([D.A, D.B, D.AB, cov.values])
    XR = np.column_stack([D.A, D.B, cov.values])
    bF = newton(y, XF)
    bR = newton(y, XR)
    muR = np.exp(XR @ bR)
    xg = D.AB.values
    w = muR
    proj = XR @ np.linalg.solve(XR.T @ (XR * w[:, None]), XR.T @ (w * xg))
    xt = xg - proj
    r = y - muR
    U = np.array([(xt[g == s] * r[g == s]).sum() for s in seasons])
    t = U.sum() / np.sqrt((U ** 2).sum())
    W = rng.choice([-1.0, 1.0], size=(nflip, len(seasons)))
    ts = (W @ U) / np.sqrt((U ** 2).sum())
    p = (1 + (np.abs(ts) >= abs(t) - 1e-12).sum()) / (nflip + 1)
    # cluster SE
    muF = np.exp(XF @ bF)
    H = XF.T @ (XF * muF[:, None])
    Hi = np.linalg.inv(H)
    S = np.array([(XF[g == s] * (y - muF)[g == s, None]).sum(0) for s in seasons])
    Gn = len(seasons)
    Vc = Hi @ (S.T @ S) @ Hi * Gn / (Gn - 1)
    se = np.sqrt(np.diag(Vc))
    idx_by = {s: np.where(g == s)[0] for s in seasons}
    boots = np.empty(nboot)
    for k in range(nboot):
        pick = rng.choice(seasons, size=len(seasons), replace=True)
        ii = np.concatenate([idx_by[s] for s in pick])
        try:
            boots[k] = newton(y[ii], XF[ii])[2]
        except Exception:
            boots[k] = np.nan
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return dict(n_events=int(y.sum()), n_seasons=len(seasons), A=bF[0], B=bF[1], gamma=bF[2],
                se_A=se[0], se_B=se[1], se_gamma=se[2], ci_gamma=(lo, hi), p=p, coef=bF)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    run = 1.0
    for rank, i in enumerate(o[::-1]):
        k = n - rank
        run = min(run, p[i] * n / k)
        q[i] = run
    return q
