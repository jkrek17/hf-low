"""Shared pieces for the jet / trough threshold analysis (see PREREGISTRATION.md).

ERA5 proxy, pipeline A tracks, the intensity framework's fix table. Nothing here reads the held-out
seasons except `evaluate`, which the analysis calls once per declared model.
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
INT = os.path.join(HERE, "..", "intensity")
sys.path.insert(0, INT)
import model as FW                                    # the framework: load(), Prep, SETS

DEV = bool(os.environ.get("JT_DEV"))                  # mechanics check: fit and "test" both inside the fit seasons
FIT_LAST = 2009 if DEV else 2014                      # fit seasons 2004..2014, held out 2015..2025
COV = ["msl", "dp12", "young", "lat", "doy_c", "doy_s"]
COV_HF = COV + ["g800"]
KNOT_PCT = np.arange(10, 90.01, 2.5)
SPL_PCT = (1, 25, 50, 75, 99)


def load_all():
    D = FW.load(os.path.join(INT, "results", "fixes_2004.csv.gz"), os.path.join(INT, "results", "env_2004.csv.gz"))
    T = pd.read_csv(os.path.join(HERE, "results", "features_2004.csv.gz"), dtype={"time": str})
    D = D.merge(T, on=["track", "time"], how="left", validate="one_to_one")
    D["hf_now"] = D.hf_now.astype(bool)
    D["hf24"] = D.hf24.astype(bool)
    D["hf48"] = D.hf48.astype(bool)
    D["bomb"] = (D.cls == 4).astype(int)
    D["elig"] = D.cls >= 0
    D["fit"] = D.season <= FIT_LAST
    if DEV:
        D = D[D.season <= 2014].reset_index(drop=True)
    return D


def subset(D, outcome, basin):
    """Rows and 0/1 outcome for one basin. BOMB: class-eligible fixes. HFON: fixes not HF at t."""
    S = D[D.basin == basin]
    if outcome == "BOMB":
        S = S[S.elig]
        y = S.bomb.values
    elif outcome == "HFON":
        S = S[~S.hf_now]
        y = S.hf24.astype(int).values
    else:
        raise ValueError(outcome)
    return S, y


def cov_cols(outcome):
    return COV_HF if outcome == "HFON" else COV


# ---------------------------------------------------------------- logistic regression
def fit_logit(X, y, iters=40, ridge=1e-6, beta0=None):
    n, p = X.shape
    b = np.zeros(p) if beta0 is None else beta0.copy()
    pen = np.full(p, ridge); pen[0] = 0.0
    prev = -np.inf
    for _ in range(iters):
        z = X @ b
        mu = 1 / (1 + np.exp(-np.clip(z, -30, 30)))
        ll = np.sum(y * z - np.logaddexp(0, z)) - 0.5 * np.sum(pen * b * b)
        g = X.T @ (y - mu) - pen * b
        H = (X * (mu * (1 - mu))[:, None]).T @ X + np.diag(pen)
        step = np.linalg.solve(H, g)
        t = 1.0
        while t > 1e-4:
            b2 = b + t * step
            z2 = X @ b2
            ll2 = np.sum(y * z2 - np.logaddexp(0, z2)) - 0.5 * np.sum(pen * b2 * b2)
            if ll2 >= ll - 1e-9:
                break
            t /= 2
        b = b2
        if abs(ll2 - prev) < 1e-7:
            break
        prev = ll2
    z = X @ b
    return b, float(np.sum(y * z - np.logaddexp(0, z)))


def prob(X, b):
    return 1 / (1 + np.exp(-np.clip(X @ b, -30, 30)))


def logloss(p, y):
    p = np.clip(p, 1e-9, 1 - 1e-9)
    return -(y * np.log(p) + (1 - y) * np.log(1 - p))


# ---------------------------------------------------------------- design
class Design:
    """Covariates standardised on the fit rows, focal variable clipped to the fit 1st-99th percentile
    and standardised. Builds the five forms: C, lin, hinge(k), spl, and the additive two-variable form."""
    def __init__(self, S, fit_mask, cols, xcols):
        self.cols = cols
        self.prep = FW.Prep(S.loc[fit_mask, cols].values.astype(float)) if cols else None
        self.xs = {}
        for c in xcols:
            v = S.loc[fit_mask, c].values
            lo, hi = np.percentile(v, [1, 99])
            vc = np.clip(v, lo, hi)
            self.xs[c] = (lo, hi, vc.mean(), vc.std())
        self.spl = {c: np.percentile(self.z(c, S.loc[fit_mask, c].values), SPL_PCT) for c in xcols}

    def cov(self, S):
        if self.prep is None:
            return np.ones((len(S), 1))
        return np.column_stack([np.ones(len(S)), self.prep(S[self.cols].values.astype(float))])

    def z(self, c, v):
        lo, hi, mu, sd = self.xs[c]
        return (np.clip(v, lo, hi) - mu) / sd

    def x(self, S, c):
        return self.z(c, S[c].values)

    def rcs(self, c, z):
        t = self.spl[c]
        k = len(t)
        sc = (t[-1] - t[0]) ** 2
        out = [z]
        for j in range(k - 2):
            b = (np.maximum(z - t[j], 0) ** 3
                 - np.maximum(z - t[k - 2], 0) ** 3 * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
                 + np.maximum(z - t[k - 1], 0) ** 3 * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2]))
            out.append(b / sc)
        return np.column_stack(out)

    def form(self, kind, S, c, k=None, C=None):
        C = self.cov(S) if C is None else C
        z = self.x(S, c)
        if kind == "C":
            return C
        if kind == "lin":
            return np.column_stack([C, z])
        if kind == "hinge":
            return np.column_stack([C, z, np.maximum(z - k, 0)])
        if kind == "spl":
            return np.column_stack([C, self.rcs(c, z)])
        raise ValueError(kind)

    def knots(self, S_fit, c):
        zf = self.x(S_fit, c)
        return np.percentile(zf, KNOT_PCT)


def best_hinge(Dg, S, y, c, rows=None):
    """Knot by maximum likelihood on the rows given (default all of S). Returns k, beta, loglik, lls."""
    idx = np.arange(len(S)) if rows is None else rows
    Sf = S.iloc[idx]
    C = Dg.cov(Sf)
    z = Dg.x(Sf, c)
    ks = np.percentile(z, KNOT_PCT)
    best = None
    lls = []
    for k in ks:
        X = np.column_stack([C, z, np.maximum(z - k, 0)])
        b, ll = fit_logit(X, y[idx])
        lls.append(ll)
        if best is None or ll > best[2]:
            best = (k, b, ll)
    return best[0], best[1], best[2], np.array(lls), ks


# ---------------------------------------------------------------- season bootstrap
def season_gain(ll_a, ll_b, seasons, B=2000, seed=0):
    """Total log-loss gain of model b over model a (positive favours b), with a season-bootstrap
    interval and p = share of resamples with gain <= 0."""
    rng = np.random.default_rng(seed)
    us = np.unique(seasons)
    g = np.array([np.sum(ll_a[seasons == s] - ll_b[seasons == s]) for s in us])
    idx = rng.integers(0, len(us), (B, len(us)))
    bs = g[idx].sum(1)
    return dict(gain=float(g.sum()), lo=float(np.percentile(bs, 5)), hi=float(np.percentile(bs, 95)),
                p=float((np.sum(bs <= 0) + 1) / (B + 1)), n_seasons=len(us))


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty(len(p))
    ranked = p[o] * len(p) / (np.arange(len(p)) + 1)
    q[o] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    return q


# ---------------------------------------------------------------- the three declared tests
def three_tests(S, y, c, cols, B=2000, seed=0, with_curves=False):
    """Fit C, lin, hinge (knot by ML on the fit seasons), spl on the fit seasons; score the held-out
    seasons once each. Returns the gains T1 (lin over C), T2 (hinge over lin), T3 (hinge over spl)."""
    fm = (S.season <= FIT_LAST).values
    Dg = Design(S, fm, cols, [c])
    Sf, Se = S[fm], S[~fm]
    yf, ye = y[fm], y[~fm]
    Cf, Ce = Dg.cov(Sf), Dg.cov(Se)
    k, bh_, llh, lls, ks = best_hinge(Dg, Sf, yf, c)
    bC, _ = fit_logit(Cf, yf)
    bL, _ = fit_logit(Dg.form("lin", Sf, c, C=Cf), yf)
    bS, _ = fit_logit(Dg.form("spl", Sf, c, C=Cf), yf)
    pe = {"C": prob(Ce, bC),
          "lin": prob(Dg.form("lin", Se, c, C=Ce), bL),
          "hinge": prob(Dg.form("hinge", Se, c, k=k, C=Ce), bh_),
          "spl": prob(Dg.form("spl", Se, c, C=Ce), bS)}
    ll = {m: logloss(p, ye) for m, p in pe.items()}
    sea = Se.season.values
    out = dict(T1=season_gain(ll["C"], ll["lin"], sea, B, seed),
               T2=season_gain(ll["lin"], ll["hinge"], sea, B, seed),
               T3=season_gain(ll["spl"], ll["hinge"], sea, B, seed),
               T3rev=season_gain(ll["hinge"], ll["spl"], sea, B, seed),
               knot_z=float(k), knot_pct=float(KNOT_PCT[int(np.argmax(lls))]),
               slope_pre=float(bh_[-2]), slope_post=float(bh_[-2] + bh_[-1]),
               b_lin=float(bL[-1]), n_fit=int(len(yf)), n_test=int(len(ye)),
               ev_fit=int(yf.sum()), ev_test=int(ye.sum()),
               ll={m: float(v.sum()) for m, v in ll.items()})
    if with_curves:
        out["_obj"] = dict(Dg=Dg, Sf=Sf, yf=yf, Se=Se, ye=ye, k=k, bh=bh_, bL=bL, bS=bS, bC=bC, lls=lls, ks=ks, pe=pe)
    return out
