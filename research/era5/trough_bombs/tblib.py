"""Shared pieces for the RA-24 trough / bomb test (see PREREGISTRATION.md). ERA5 proxy, pipeline A.

Builds the analysis table, the pooled logistic model, the season-clustered (CR1) inference and the season
bootstrap. Nothing here picks a trough definition: D1 = `trough_up`, D2 = `trough_up_eddy`, both always run.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
INT = os.path.join(HERE, "..", "intensity")
sys.path.insert(0, INT)
sys.path.insert(0, os.path.join(HERE, "..", "jet_trough"))
import model as FW                      # Prep, EDGES
from scipy import stats

COV = ["msl", "dp12", "young", "lat", "doy_c", "doy_s"]
DEFS = {"D1": "trough_up", "D2": "trough_up_eddy"}
REF = (2004, 2014)
NDR_MAX = 3.0
YR0 = 1991


def prep_fixes(F):
    """Same derived columns and BOMB definition as `intensity/model.load` (without the environment table)."""
    D = F.copy()
    D["young"] = D.dp12.isna().astype(float)
    D["dp12"] = D.dp12.fillna(0.0)
    D["pac"] = (D.basin == "pac").astype(float)
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    doy = t.dt.dayofyear.values
    D["doy_c"], D["doy_s"] = np.cos(2 * np.pi * doy / 365.25), np.sin(2 * np.pi * doy / 365.25)
    D["month"] = D.time.str[4:6].astype(int)
    ok = D.ndr24.notna() & (D.ndr24.abs() <= NDR_MAX)
    D["elig"] = ok
    D["bomb"] = (ok & (D.ndr24 >= 1.0)).astype(int)
    D["yr"] = (D.season - YR0) / 10.0
    return D.reset_index(drop=True)


def add_lead24(D):
    """BOMB at the same track's fix 24 h later (outcome window (t+24, t+48 h]); NaN where there is none."""
    t = pd.to_datetime(D.time, format="%Y%m%d%H")
    key = D.assign(t=t)[["track", "t", "elig", "bomb"]]
    nxt = key.assign(t=key.t - pd.Timedelta(hours=24)).rename(columns={"elig": "elig_l", "bomb": "bomb_l"})
    m = D.assign(t=t).merge(nxt, on=["track", "t"], how="left")
    D["elig_l24"] = m.elig_l.fillna(False).astype(bool).values
    D["bomb_l24"] = m.bomb_l.fillna(0).astype(int).values
    return D


class Std:
    """Reference-season standardisation: pooled covariates (Prep) and the per-basin trough variable."""
    def __init__(self, R):
        R = R[R.elig]
        self.prep = FW.Prep(R[COV].values.astype(float))
        self.focal = {}
        for col in set(DEFS.values()) | {"jet_up", "trough_loc"}:
            for b in ("atl", "pac"):
                v = R.loc[(R.basin == b), col].dropna().values
                lo, hi = np.percentile(v, [1, 99])
                vc = np.clip(v, lo, hi)
                self.focal[(col, b)] = (lo, hi, vc.mean(), vc.std())

    def z(self, S, col):
        out = np.full(len(S), np.nan)
        for b in ("atl", "pac"):
            m = (S.basin == b).values
            lo, hi, mu, sd = self.focal[(col, b)]
            out[m] = (np.clip(S.loc[m, col].values, lo, hi) - mu) / sd
        return out

    def sd(self, col, b):
        return self.focal[(col, b)][3]


def design(S, std, col, pooled=True, year=True, extra=(), inter=False):
    """Columns: 1, standardised covariates, [pac], [yr], extras (jet), focal z last."""
    C = std.prep(S[COV].values.astype(float))
    X = [np.ones(len(S)), *C.T]
    if pooled:
        X.append(S.pac.values)
        if inter:
            X += [S.pac.values * c for c in C.T]
    if year:
        X.append(S.yr.values)
    for e in extra:
        X.append(std.z(S, e))
    if col is not None:
        X.append(std.z(S, col))
    return np.column_stack(X)


def fit_w(X, y, w=None, iters=50, ridge=1e-8):
    """Weighted unpenalised (ridge 1e-8) logistic regression by Newton steps with step halving."""
    n, p = X.shape
    w = np.ones(n) if w is None else w
    b = np.zeros(p)
    pen = np.full(p, ridge); pen[0] = 0
    def ll(b):
        z = X @ b
        return np.sum(w * (y * z - np.logaddexp(0, z))) - 0.5 * np.sum(pen * b * b)
    cur = ll(b)
    for _ in range(iters):
        mu = 1 / (1 + np.exp(-np.clip(X @ b, -30, 30)))
        g = X.T @ (w * (y - mu)) - pen * b
        H = (X * (w * mu * (1 - mu))[:, None]).T @ X + np.diag(pen)
        step = np.linalg.solve(H, g)
        t = 1.0
        while t > 1e-4:
            new = ll(b + t * step)
            if new >= cur - 1e-10:
                break
            t /= 2
        b = b + t * step
        if abs(new - cur) < 1e-9:
            cur = new
            break
        cur = new
    return b, cur


def cluster_inference(X, y, b, clusters, k=-1):
    """CR1 variance for coefficient k; t reference on G-1 df. Returns se, p (two-sided)."""
    mu = 1 / (1 + np.exp(-np.clip(X @ b, -30, 30)))
    H = (X * (mu * (1 - mu))[:, None]).T @ X
    Hi = np.linalg.inv(H)
    sc = X * (y - mu)[:, None]
    u = np.unique(clusters)
    G = len(u)
    S = np.zeros((G, X.shape[1]))
    idx = np.searchsorted(u, clusters)
    np.add.at(S, idx, sc)
    M = S.T @ S
    n, p = X.shape
    V = Hi @ M @ Hi * (G / (G - 1)) * ((n - 1) / (n - p))
    se = np.sqrt(V[k, k])
    t = b[k] / se
    return se, 2 * stats.t.sf(abs(t), G - 1), V


def summarize(S, y, X, k=-1):
    """Odds ratio per SD of the focal variable (last column), CR1 95% interval, two-sided p."""
    b, ll = fit_w(X, y)
    se, p, V = cluster_inference(X, y, b, S.season.values, k)
    G = S.season.nunique()
    tq = stats.t.ppf(0.975, G - 1)
    return dict(beta=b[k], se=se, OR=np.exp(b[k]), lo=np.exp(b[k] - tq * se), hi=np.exp(b[k] + tq * se), p=p,
                n=len(y), events=int(y.sum()), seasons=G)


def season_boot(S, y, X, B=2000, seed=0, k=-1):
    """Resample whole seasons (weights = draw counts); returns the 2.5/97.5 percentiles of the OR."""
    rng = np.random.default_rng(seed)
    u = np.unique(S.season.values)
    idx = np.searchsorted(u, S.season.values)
    out = []
    b0, _ = fit_w(X, y)
    for _ in range(B):
        c = np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u))
        w = c[idx].astype(float)
        b, _ = fit_w(X, y, w)
        out.append(b[k])
    out = np.array(out)
    return float(np.exp(np.percentile(out, 2.5))), float(np.exp(np.percentile(out, 97.5)))


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty(len(p))
    r = p[o] * len(p) / (np.arange(len(p)) + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    return q


def load_tables(work):
    """New seasons (1979-80..2003-04) and reference seasons (2004-05..2014-15) with trough features.
    `work` holds fixes_1979_2003.csv.gz and features_1979_2003.csv.gz; the reference fixes and PR 78's
    committed features are read from the repository."""
    R0 = os.path.join(HERE, "..")
    N = pd.read_csv(os.path.join(work, "fixes_1979_2003.csv.gz"), dtype={"time": str})
    NF = pd.read_csv(os.path.join(work, "features_1979_2003.csv.gz"), dtype={"time": str})
    N = N.merge(NF, on=["track", "time"], how="left", validate="one_to_one")
    cols = ["track", "time", "basin", "season", "lat", "lon", "msl", "dp12", "ndr24", "age"]
    R = pd.read_csv(os.path.join(R0, "intensity", "results", "fixes_2004.csv.gz"), dtype={"time": str}, usecols=cols)
    R = R[R.season <= REF[1]]
    RF = pd.read_csv(os.path.join(R0, "jet_trough", "results", "features_2004.csv.gz"), dtype={"time": str})
    R = R.merge(RF, on=["track", "time"], how="left", validate="one_to_one")
    return add_lead24(prep_fixes(N)), add_lead24(prep_fixes(R))
