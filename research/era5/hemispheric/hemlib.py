"""Library for the hemispheric pattern search (ERA5 proxy). See PREREGISTRATION.md; nothing here
chooses anything from an outcome except the ridge penalty, by leave-one-season-out inside the
training seasons."""
import os, json, datetime as dt
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
VARS = ("z500", "u250", "sst")
ALLVARS = ("z500", "u250", "mslp", "sst")
LAT = (-87.1875 + 5.625 * np.arange(32))[20:32]
NWEEK = 30
K = 10
LAMBDAS = [1, 3, 10, 30, 100, 300, 1000, 3000, np.inf]       # np.inf = no PCs
SEED = 20261008
BASINS = ("atl", "pac")
DEPTH_CUT = {"atl": 966.2, "pac": 965.0}
HF_KT = 71.7


# ------------------------------------------------------------------ fields and anomalies
def load_fields(workdir, seasons=range(1979, 2026)):
    F = {}
    for s in seasons:
        z = np.load(os.path.join(workdir, f"season_{s}.npz"))
        dates = [dt.date.fromisoformat(x) for x in z["dates"]]
        F[s] = (dates, {v: z[v] for v in ALLVARS})
    return F


def _mdkey(d):
    return (d.month, d.day)


def climatology(F, clim_seasons=range(1991, 2020)):
    """Calendar-day mean over seasons 1991-92 .. 2019-20, 15-day running mean, per variable."""
    ref = [d for d in F[2001][0] if not (d.month == 2 and d.day == 29)]       # non-leap sequence
    keys = [_mdkey(d) for d in ref]
    pos = {k: i for i, k in enumerate(keys)}
    C = {}
    for v in ALLVARS:
        acc = np.zeros((len(keys), 12, 64)); n = np.zeros(len(keys))
        for s in clim_seasons:
            dates, o = F[s]
            for j, d in enumerate(dates):
                k = _mdkey(d)
                if k in pos:
                    acc[pos[k]] += np.nan_to_num(o[v][j], nan=0.0); n[pos[k]] += 1
        c = acc / n[:, None, None]
        # running mean, truncated at the ends
        w = 15
        sm = np.empty_like(c)
        for i in range(len(c)):
            lo, hi = max(0, i - 7), min(len(c), i + 8)
            sm[i] = c[lo:hi].mean(0)
        C[v] = sm
    return keys, pos, C


def anomalies(F, trend_seasons=range(1979, 2022)):
    """Per season: dict var -> [ndays, 12, 64] anomaly (clim and linear trend removed). SST NaN preserved."""
    keys, pos, C = climatology(F)
    def clim_at(v, d):
        if d.month == 2 and d.day == 29:
            return 0.5 * (C[v][pos[(2, 28)]] + C[v][pos[(3, 1)]])
        return C[v][pos[_mdkey(d)]]
    A = {}
    for s, (dates, o) in F.items():
        A[s] = {v: np.stack([o[v][j] - clim_at(v, d) for j, d in enumerate(dates)]) for v in ALLVARS}
    # trend of the Oct-Apr seasonal mean
    smean = {v: [] for v in ALLVARS}
    ts = list(trend_seasons)
    for s in ts:
        dates = F[s][0]
        m = [j for j, d in enumerate(dates) if dt.date(s, 10, 1) <= d <= dt.date(s + 1, 4, 30)]
        for v in ALLVARS:
            smean[v].append(np.mean(A[s][v][m], axis=0))
    x = np.array(ts, float); x0 = x - x.mean()
    for v in ALLVARS:
        Y = np.stack(smean[v])                                        # season, 12, 64
        b = (x0[:, None, None] * (Y - Y.mean(0))).sum(0) / (x0 ** 2).sum()
        a = Y.mean(0)
        for s in A:
            A[s][v] = A[s][v] - (a + b * (s - x.mean()))
    return A


def sst_mask(A):
    bad = np.zeros((12, 64), bool)
    for s in A:
        bad |= np.isnan(A[s]["sst"]).any(0)
    return ~bad


def week_start(s, k):
    return dt.date(s, 10, 1) + dt.timedelta(days=7 * k)


def window_map(F, A, s, k, v, a, b):
    """Mean anomaly of var v over days S+a .. S+b relative to the week start S."""
    dates = F[s][0]
    S = week_start(s, k)
    i0 = dates.index(S + dt.timedelta(days=a)); i1 = dates.index(S + dt.timedelta(days=b))
    return np.mean(A[s][v][i0:i1 + 1], axis=0)


def weekly_maps(F, A, seasons, a, b):
    """dict var -> [nseason*NWEEK, 12, 64], rows ordered season then week."""
    return {v: np.stack([window_map(F, A, s, k, v, a, b) for s in seasons for k in range(NWEEK)]) for v in VARS}


# ------------------------------------------------------------------ EOFs
class EOFs:
    """Leading K EOFs per field of training weekly maps; project any weekly maps to unit-variance scores."""

    def __init__(self, maps, mask, k=K):
        self.mask = mask
        w = np.sqrt(np.cos(np.radians(LAT)))[:, None] * np.ones((1, 64))
        self.w = w
        self.mean, self.vec, self.sd = {}, {}, {}
        for v in VARS:
            X = self._flat(maps[v], v)
            self.mean[v] = X.mean(0)
            U, S, Vt = np.linalg.svd(X - self.mean[v], full_matrices=False)
            self.vec[v] = Vt[:k]
            sc = (X - self.mean[v]) @ Vt[:k].T
            self.sd[v] = sc.std(0, ddof=1)
            self.var_explained = getattr(self, "var_explained", {})
            self.var_explained[v] = float((S[:k] ** 2).sum() / (S ** 2).sum())

    def _flat(self, M, v):
        X = M * self.w
        if v == "sst":
            return X[:, self.mask]
        return X.reshape(len(M), -1)

    def scores(self, maps):
        return np.hstack([((self._flat(maps[v], v) - self.mean[v]) @ self.vec[v].T) / self.sd[v] for v in VARS])

    def to_map(self, beta, v):
        """Map (12 x 64) of sum_j beta_j * EOF_j (in anomaly units; weights undone) for the PC block of var v."""
        i = VARS.index(v)
        b = beta[i * K:(i + 1) * K] / self.sd[v]
        flat = b @ self.vec[v]
        if v == "sst":
            m = np.full((12, 64), np.nan); m[self.mask] = flat
        else:
            m = flat.reshape(12, 64)
        return m / self.w


# ------------------------------------------------------------------ Poisson with ridge
def deviance(y, mu):
    y = np.asarray(y, float)
    t = np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0)
    return 2 * float(np.sum(t - (y - mu)))


def fit_pois(X, y, pen, beta0=None, iters=60):
    """Poisson log link, penalty 0.5 * sum(pen * beta^2). Newton with step halving."""
    n, p = X.shape
    pen = np.asarray(pen, float)
    b = np.zeros(p) if beta0 is None else beta0.copy()
    if beta0 is None:
        # month dummies come first and sum to one: start every free coefficient at 0 and the dummies at log(mean y)
        b[:] = 0.0
    def nll(bb):
        eta = np.clip(X @ bb, -30, 30)
        return float(np.sum(np.exp(eta) - y * eta) + 0.5 * np.sum(pen * bb * bb))
    if beta0 is None:
        m0 = max(y.mean(), 1e-3)
        b[:7] = np.log(m0)
    cur = nll(b)
    for _ in range(iters):
        eta = np.clip(X @ b, -30, 30); mu = np.exp(eta)
        g = X.T @ (y - mu) - pen * b
        H = (X * mu[:, None]).T @ X + np.diag(pen) + 1e-8 * np.eye(p)
        step = np.linalg.solve(H, g)
        t = 1.0
        while t > 1e-4:
            nb = b + t * step
            new = nll(nb)
            if new <= cur + 1e-10:
                break
            t /= 2
        if t <= 1e-4:
            break
        done = abs(cur - new) < 1e-9 * (1 + abs(cur))
        b, cur = nb, new
        if done:
            break
    return b


def design(df, extra_cols, pcs=None):
    """Month dummies (Oct..Apr), log1p(prev), extra unpenalised columns, then PCs. Returns X, number unpenalised."""
    m = pd.get_dummies(df["month"]).reindex(columns=[10, 11, 12, 1, 2, 3, 4], fill_value=0).values.astype(float)
    cols = [m, df[["lprev"]].values]
    for c in extra_cols:
        cols.append(df[[c]].values)
    nu = sum(c.shape[1] for c in cols)
    if pcs is not None:
        cols.append(pcs)
    return np.hstack(cols), nu


def fit_model(df, extra_cols, pcs, lam, y_col="y"):
    X, nu = design(df, extra_cols, pcs if np.isfinite(lam) else None)
    pen = np.zeros(X.shape[1]); pen[nu:] = lam
    b = fit_pois(X, df[y_col].values.astype(float), pen)
    return dict(beta=b, nu=nu, extra=list(extra_cols), lam=float(lam))


def predict(model, df, pcs, parts=False):
    X, nu = design(df, model["extra"], pcs if np.isfinite(model["lam"]) else None)
    b = model["beta"]
    if parts:
        return X[:, :nu] @ b[:nu], (X[:, nu:] @ b[nu:] if X.shape[1] > nu else np.zeros(len(df)))
    return np.exp(np.clip(X @ b, -30, 30))


def loso_lambda(df, extra_cols, pcs, y_col="y", lambdas=LAMBDAS):
    """Leave-one-season-out deviance per lambda inside the training frame."""
    seasons = sorted(df["season"].unique())
    tot = {}
    for lam in lambdas:
        d = 0.0
        for s in seasons:
            tr, te = (df["season"] != s).values, (df["season"] == s).values
            mdl = fit_model(df[tr], extra_cols, pcs[tr], lam, y_col)
            d += deviance(df.loc[te, y_col].values, predict(mdl, df[te], pcs[te]))
        tot[lam] = d if np.isfinite(d) else float('inf')
    best = min(tot, key=tot.get)
    return best, tot


# ------------------------------------------------------------------ SOM
def som_fit(Xtr, shape=(3, 3), epochs=200, r0=1.5, r1=0.5, seed=SEED):
    rng = np.random.default_rng(seed)
    nr, nc = shape
    gx, gy = np.meshgrid(np.arange(nc), np.arange(nr))
    G = np.stack([gy.ravel(), gx.ravel()], 1).astype(float)
    U, S, Vt = np.linalg.svd(Xtr - Xtr.mean(0), full_matrices=False)
    s1, s2 = S[0] / np.sqrt(len(Xtr)), S[1] / np.sqrt(len(Xtr))
    W = np.zeros((nr * nc, Xtr.shape[1]))
    for i, (r, c) in enumerate(G):
        W[i] = Xtr.mean(0) + 2 * ((c / (nc - 1)) - 0.5) * s1 * Vt[0] + 2 * ((r / (nr - 1)) - 0.5) * s2 * Vt[1]
    for ep in range(epochs):
        rad = r0 + (r1 - r0) * ep / (epochs - 1)
        bmu = np.argmin(((Xtr[:, None, :] - W[None]) ** 2).sum(2), axis=1)
        d2 = ((G[:, None, :] - G[None, :, :]) ** 2).sum(2)                    # node-node
        hh = np.exp(-d2 / (2 * rad ** 2))                                    # [node, bmu node]
        num = np.zeros_like(W); den = np.zeros(len(W))
        for j in range(len(W)):
            m = bmu == j
            if m.any():
                num += hh[:, j][:, None] * Xtr[m].sum(0)[None]
                den += hh[:, j] * m.sum()
        W = np.where(den[:, None] > 0, num / np.maximum(den[:, None], 1e-12), W)
    return W


def som_assign(W, X):
    return np.argmin(((X[:, None, :] - W[None]) ** 2).sum(2), axis=1)


# ------------------------------------------------------------------ counts
def load_archive_counts(seasons, extra_before=1):
    d = json.load(open(os.path.join(REPO, "docs/data/hf-lows.json")))
    f = {k: i for i, k in enumerate(d["lowFields"])}
    days = {b: {} for b in BASINS}
    for r in d["lows"]:
        if r[f["cls"]] != "low" or r[f["basin"]] not in BASINS:
            continue
        s = str(r[f["start"]])
        day = dt.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
        days[r[f["basin"]]][day] = days[r[f["basin"]]].get(day, 0) + 1
    return days


def load_proxy_days(kind):
    """Pipeline A per-day counts by `start`: kind 'hf' (gust800 >= 71.7 kt) or 'depth' (minp cut)."""
    T = pd.read_csv(os.path.join(REPO, "research/era5/hf_history/results/all_tracks.csv.gz"), dtype={"start": str})
    T = T[T.basin.isin(BASINS)]
    days = {b: {} for b in BASINS}
    for b in BASINS:
        t = T[T.basin == b]
        sel = t[t.gust800_kt >= HF_KT] if kind == "hf" else t[t.minp <= DEPTH_CUT[b]]
        for st in sel.start:
            d = dt.date(int(st[:4]), int(st[4:6]), int(st[6:8]))
            days[b][d] = days[b].get(d, 0) + 1
    return days


def weekly_table(days, seasons, basin):
    rows = []
    for s in seasons:
        for k in range(NWEEK):
            S = week_start(s, k)
            y = sum(days[basin].get(S + dt.timedelta(days=i), 0) for i in range(7))
            p = sum(days[basin].get(S - dt.timedelta(days=i), 0) for i in range(1, 8))
            mid = S + dt.timedelta(days=3)
            rows.append(dict(season=s, week=k, month=mid.month, y=y, prev=p, lprev=np.log1p(p)))
    return pd.DataFrame(rows)
