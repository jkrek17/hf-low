"""RA-9 first pass: does the gust in HF-strength fixes follow gradient-wind scaling? ERA5 PROXY, pipeline A.

Plan: PREREGISTRATION.md (committed first, 3acbd3e). Inputs are committed tables only:
  hf_structure/results/fixes.csv        5,983 HF-strength fixes (g800 >= 71.7 kt by construction)
  tele_intensity/results/event_table.csv  per-storm ring and clim MSLP at the deepest fix
  docs/data/teleconnections.json         monthly ONI, daily PNA
usage: scaling.py [NBOOT=2000]
"""
import os, sys, json, numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import log_ndtr, ndtr
from scipy.stats import norm
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__)); ERA = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(ERA))
OUT = os.path.join(HERE, "results")
CUT = np.log(71.7)
NBOOT = int(sys.argv[1]) if len(sys.argv) > 1 else 2000


def load():
    f = pd.read_csv(os.path.join(ERA, "hf_structure/results/fixes.csv"))
    e = pd.read_csv(os.path.join(ERA, "tele_intensity/results/event_table.csv"))[["track", "ring", "clim"]]
    f = f.merge(e, on="track", how="left")
    t = pd.to_datetime(f.time.astype(str), format="%Y%m%d%H")
    f["date"] = t
    T = json.load(open(os.path.join(REPO, "docs/data/teleconnections.json")))
    y0, m0 = map(int, T["oni"]["start"].split("-"))
    oni = {}
    for i, v in enumerate(T["oni"]["values"]):
        if v is not None:
            oni[(y0 + (m0 - 1 + i) // 12, (m0 - 1 + i) % 12 + 1)] = v
    pm = pd.Timestamp(T["pna"]["start"]); pv = np.array([np.nan if v is None else v for v in T["pna"]["values"]], float)

    def pna(d, a, b):  # mean over days d-a .. d-b (a > b), UTC calendar days
        k = (d.normalize() - pm).days
        return np.nanmean(pv[k - a:k - b + 1]) if k - a >= 0 and k - b + 1 <= len(pv) else np.nan
    prev = t - pd.DateOffset(months=1)
    f["oni_raw"] = [oni.get((d.year, d.month), np.nan) for d in prev]
    f["pna_raw"] = [pna(d, 10, 4) for d in t]
    f["pna0_raw"] = [pna(d, 0, 0) for d in t]
    f["Dp"] = f.ring - f.msl
    f["Dp_clim"] = f.clim - f.msl
    f["Dp_raw"] = 1013.0 - f.msl
    f["R"] = np.maximum(f.gmax_r, 25.0)
    f["Req"] = np.sqrt(f.a_g48 / (np.pi * f.own_ocean_frac))
    f["y"] = np.log(f.g800)
    f["month"] = f.date.dt.month
    return f


def design(d, dep="Dp", size="R", extra=(), idx=("oni", "pna"), with_size=True):
    cols = {"const": np.ones(len(d)), "lnDp": np.log(d[dep].values)}
    if with_size:
        cols["lnR"] = np.log(d[size].values)
    cols["lnsinlat"] = np.log(np.sin(np.radians(d.lat.values)))
    for k in idx:
        cols[k] = d[k].values
    for k in extra:
        cols[k] = d[k].values
    X = pd.DataFrame(cols)
    for c in X.columns[1:]:
        X[c] = X[c] - X[c].mean()  # centre (slopes unchanged; conditions the optimiser)
    return X


def _inner(X, y, sg, beta0):
    """Truncated-normal location problem at fixed sigma: strictly concave in beta, so damped Newton converges."""
    def ll_g(beta):
        mu = X @ beta; z = (y - mu) / sg; a = (CUT - mu) / sg; lq = log_ndtr(-a)
        ll = (-0.5 * z * z - np.log(sg) - 0.5 * np.log(2 * np.pi) - lq).sum()
        lam = np.exp(-0.5 * a * a - 0.5 * np.log(2 * np.pi) - lq)  # phi(a)/(1-Phi(a))
        return ll, X.T @ ((z - lam) / sg), (1.0 - lam * (lam - a)) / (sg * sg)
    beta = beta0.copy(); ll, g, w = ll_g(beta)
    for _ in range(60):
        step = np.linalg.solve(X.T @ (X * w[:, None]), g)
        t = 1.0
        while t > 1e-8:
            ll2, g2, w2 = ll_g(beta + t * step)
            if ll2 >= ll - 1e-10:
                break
            t *= 0.5
        beta = beta + t * step; done = abs(ll2 - ll) < 1e-9; ll, g, w = ll2, g2, w2
        if done:
            break
    return beta, ll


def trunc_fit(X, y, x0=None):
    """Truncated-normal MLE of ln g800 (truncation point CUT known): the pre-registered primary estimator.
    Profile likelihood: for each sigma the location problem is log-concave (_inner); sigma is found by bounded 1-D search.
    (A joint BFGS on uncentred covariates stopped early on this flat ridge and gave start-dependent answers; see README.)"""
    from scipy.optimize import minimize_scalar
    X = np.asarray(X, float)
    state = {"beta": np.linalg.lstsq(X, y, rcond=None)[0] if x0 is None else np.asarray(x0[:-1], float)}

    def f(ls):
        beta, ll = _inner(X, y, np.exp(ls), state["beta"]); state["beta"] = beta
        return -ll
    r = minimize_scalar(f, bounds=(np.log(0.03), np.log(1.0)), method="bounded", options={"xatol": 1e-5})
    beta, _ = _inner(X, y, np.exp(r.x), state["beta"])
    return np.append(beta, r.x)


def gpd_fit(X, e, x0=None):
    """Generalised Pareto regression for the exceedance e = g800 - 71.7 kt (>= 0): log sigma = X b, common shape xi.
    Returns [b..., xi]. X must have its non-constant columns centred (done in design()) so the optimiser is well scaled."""
    X = np.asarray(X, float); n, p = X.shape
    if x0 is None:
        x0 = np.zeros(p + 1); x0[0] = np.log(max(e.mean(), 1e-3)); x0[-1] = 0.1

    def nll(th):
        b, xi = th[:-1], th[-1]; ls = X @ b; sg = np.exp(ls); z = e / sg
        if abs(xi) < 1e-6:
            return float((ls + z).sum())
        w = 1 + xi * z
        if np.any(w <= 0):
            return 1e12
        return float((ls + (1 + 1 / xi) * np.log(w)).sum())
    bounds = [(None, None)] * p + [(-0.45, 1.0)]
    r = minimize(nll, x0, method="L-BFGS-B", bounds=bounds, options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-7})
    return r.x


def trunc_loglik_fit(X, y, x0):
    """Truncated-normal MLE (the pre-registered estimator); kept only to document that it is not identified here."""
    X = np.asarray(X, float)

    def nll(th):
        b, ls = th[:-1], th[-1]; sg = np.exp(ls); mu = X @ b; z = (y - mu) / sg; a = (CUT - mu) / sg
        return float(-(-0.5 * z * z - ls - 0.5 * np.log(2 * np.pi) - log_ndtr(-a)).sum())
    r = minimize(nll, x0, method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-8})
    return r.x, -r.fun


def ols_fit(X, y):
    return np.append(np.linalg.lstsq(np.asarray(X, float), y, rcond=None)[0], np.nan)


_G = {}


def _boot_one(args):
    seed, = args
    rng = np.random.default_rng(seed); G = _G
    ss = rng.integers(0, len(G["seas"]), len(G["seas"]))
    idx = np.concatenate([G["rows"][i] for i in ss])
    X = G["X"][idx]; y = G["y"][idx]
    try:
        return G["fit"](X, y, G["x0"])[:G["np"]]
    except Exception:
        return np.full(G["np"], np.nan)


def _fit_wrap(fit):
    def f(X, y, x0):
        return fit(X, y, x0) if fit in (gpd_fit, trunc_fit) else fit(X, y)
    return f


def run(name, d, X, y, fit, nboot, pool):
    """Point estimate plus season-block bootstrap. Returns DataFrame term, est, se, lo, hi, p."""
    est = fit(X.values, y) if fit is ols_fit else fit(X.values, y)
    p = X.shape[1]
    _G.clear()
    seas = np.sort(d.season.unique()); rows = [np.where(d.season.values == s)[0] for s in seas]
    _G.update(seas=seas, rows=rows, X=X.values, y=y, fit=_fit_wrap(fit), x0=est if fit in (gpd_fit, trunc_fit) else None, np=p)
    seeds = [(1000 + i,) for i in range(nboot)]
    with Pool(os.cpu_count()) as pl:  # fork after _G is set
        B = np.array(pl.map(_boot_one, seeds, chunksize=10))
    se = np.nanstd(B, axis=0, ddof=1)
    lo, hi = np.nanpercentile(B, [2.5, 97.5], axis=0)
    z = est[:p] / se; pv = 2 * norm.sf(np.abs(z))
    df = pd.DataFrame(dict(model=name, term=X.columns, est=est[:p], se=se, lo=lo, hi=hi, p=pv, n_fix=len(y),
                           n_boot_ok=int(np.isfinite(B[:, 0]).sum()), xi=(est[p] if len(est) > p else np.nan)))
    return df, B


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    r = p[o] * n / (np.arange(n) + 1); r = np.minimum.accumulate(r[::-1])[::-1]; q[o] = np.minimum(r, 1)
    return q


def trunc_diag(f):
    """S0: the pre-registered truncated-normal MLE from several starts. Documents that it is not identified on these data."""
    rows = []
    for b in ("atl", "pac"):
        d = f[f.basin == b]
        X = design(d).values; y = d.y.values
        b0 = np.linalg.lstsq(X, y, rcond=None)[0]
        for s in range(6):
            rng = np.random.default_rng(s)
            x0 = np.append(b0 + rng.normal(size=len(b0)) * 0.02 * (s > 0), np.log(0.06) + rng.normal() * 0.3 * (s > 0))
            th, ll = trunc_loglik_fit(X, y, x0)
            rows.append(dict(basin=b, start=s, loglik=ll, sigma=float(np.exp(th[-1])), lnDp=th[1], lnR=th[2], lnsinlat=th[3],
                             oni=th[4], pna=th[5]))
    D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, "trunc_normal_multistart.csv"), index=False)
    print(D.round(4).to_string())


def prep(f, basin, **flt):
    d = f[f.basin == basin].copy() if basin in ("atl", "pac") else f.copy()
    return d


def main():
    os.makedirs(OUT, exist_ok=True)
    f = load()
    n0 = len(f)
    bad = {"Dp<5": int((f.Dp < 5).sum()), "oni_nan": int(f.oni_raw.isna().sum()), "pna_nan": int(f.pna_raw.isna().sum()),
           "ring_nan": int(f.ring.isna().sum())}
    f = f[(f.Dp >= 5) & f.oni_raw.notna() & f.pna_raw.notna()].copy()
    print("dropped:", bad, "kept", len(f), "of", n0)
    for b in ("atl", "pac"):  # standardise within basin, once
        m = f.basin == b
        for k in ("oni", "pna", "pna0"):
            f.loc[m, k] = (f.loc[m, k + "_raw"] - f.loc[m, k + "_raw"].mean()) / f.loc[m, k + "_raw"].std()
    meta = dict(dropped=bad, n_kept=len(f), sd={b: {k: float(f[f.basin == b][k + "_raw"].std()) for k in ("oni", "pna")}
                                                for b in ("atl", "pac")})
    desc = []
    for b in ("atl", "pac"):
        d = f[f.basin == b]
        desc.append(dict(basin=b, n_fix=len(d), n_storm=d.track.nunique(), n_season=d.season.nunique(),
                         **{f"med_{k}": float(d[k].median()) for k in ("g800", "Dp", "R", "lat")},
                         corr_lnDp_lnR=float(np.corrcoef(np.log(d.Dp), np.log(d.R))[0, 1]),
                         corr_lnDp_g=float(np.corrcoef(np.log(d.Dp), d.y)[0, 1]),
                         corr_lnR_g=float(np.corrcoef(np.log(d.R), d.y)[0, 1]),
                         corr_lnsinlat_g=float(np.corrcoef(np.log(np.sin(np.radians(d.lat))), d.y)[0, 1]),
                         corr_lnDp_lnsinlat=float(np.corrcoef(np.log(d.Dp), np.log(np.sin(np.radians(d.lat))))[0, 1])))
    pd.DataFrame(desc).to_csv(os.path.join(OUT, "descriptive.csv"), index=False)
    json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=1)
    pass  # trunc_diag(f): superseded by the profile-likelihood fit

    jobs = []  # (label, basin, subset-mask fn, design kwargs, fit, kind)
    pool = None
    if True:
        allr = []

        def go(label, d, kind, **kw):
            X = design(d, **{k: v for k, v in kw.items() if k not in ("fit", "yv", "ylog")})
            fit = kw.get("fit", trunc_fit)
            if kind == "size":
                X = X.drop(columns=["lnR"], errors="ignore")
                fit = ols_fit
                yv = np.log(d[kw.get("size", "R")].values); fac = 1.0; scale = "ols_lnR"
            elif fit is ols_fit or fit is trunc_fit:  # primary (truncated normal) and OLS contrast, both on ln g800
                yv = d.y.values; fac = 1.0; scale = "ols_lng" if fit is ols_fit else "trunc_lng"
            else:  # S10: GPD regression of the exceedance over 71.7 kt
                yv = np.maximum(d.g800.values - 71.7, 0.0)
                # d E[g]/d x = mean exceedance * beta, so the elasticity of the mean gust is beta * mean(e) / mean(g)
                fac = float(yv.mean() / d.g800.mean()); scale = "gpd_logscale"
            df, B = run(label, d, X, yv, fit, NBOOT, pool)
            df["kind"] = kind; df["scale"] = scale; df["fac"] = fac
            allr.append(df)
            return df, B
        for b in ("atl", "pac"):
            d = f[f.basin == b].copy()
            tag = b
            # --- primary
            go(f"M1_{tag}", d, "gust", idx=())
            dfM2, B2 = go(f"M2_{tag}", d, "gust")
            go(f"M3_{tag}", d, "size", size="R")
            # ONI shrink when lnR enters: M2 without size
            dfN, BN = go(f"M2nosize_{tag}", d, "gust", with_size=False)
            i = list(dfM2.term).index("oni"); j = list(dfN.term).index("oni")
            diff = B2[:, i] - BN[:, j]
            allr.append(pd.DataFrame(dict(model=f"H5shrink_{tag}", term="oni_withR_minus_noR", est=dfM2.est[i] - dfN.est[j],
                                          se=np.nanstd(diff, ddof=1), lo=np.nanpercentile(diff, 2.5), hi=np.nanpercentile(diff, 97.5),
                                          p=np.nan, n_fix=len(d), n_boot_ok=int(np.isfinite(diff).sum()), kind="gust", scale="gpd_logscale", fac=dfM2.fac.iloc[0] if "fac" in dfM2 else np.nan), index=[0]))
            # --- secondary
            go(f"S1_ols_{tag}", d, "gust", fit=ols_fit)
            go(f"S10_gpd_{tag}", d, "gust", fit=gpd_fit)
            deep = d.sort_values("msl").groupby("track").head(1)
            go(f"S2_onefix_{tag}", deep, "gust")
            m3 = ~(((d.gmax_coast_km < 100) | d.terrain.astype(bool)) | ((d.lat >= 60) & (b == "atl")))
            go(f"S3_noterrain_{tag}", d[m3], "gust")
            go(f"S4a_Dpclim_{tag}", d[d.Dp_clim >= 5], "gust", dep="Dp_clim")  # same >= 5 hPa rule
            go(f"S4b_Dpraw_{tag}", d[d.Dp_raw >= 5], "gust", dep="Dp_raw")
            d5 = d[d.Req > 0]
            go(f"S5_Req_{tag}", d5, "gust", size="Req")
            go(f"S5size_Req_{tag}", d5, "size", size="Req")
            mm = pd.get_dummies(d.month, prefix="m", drop_first=True).astype(float)
            d7 = pd.concat([d, mm.set_index(d.index)], axis=1)
            go(f"S7_month_{tag}", d7, "gust", extra=tuple(mm.columns))
            for st in ("deepening", "mature", "filling"):
                go(f"S8_{st}_{tag}", d[d.stage == st], "gust")
            d9 = d[d.pna0.notna()].copy(); d9["pna"] = d9["pna0"]
            go(f"S9_pnasame_{tag}", d9, "gust")
        d6 = f.copy(); d6["atl"] = (d6.basin == "atl").astype(float)
        go("S6_pooled", d6, "gust", extra=("atl",))
        go("S6size_pooled", d6, "size", extra=("atl",))
    R = pd.concat(allr, ignore_index=True)
    prim = R[R.model.isin([f"M2_{b}" for b in ("atl", "pac")] + [f"M3_{b}" for b in ("atl", "pac")]) & (R.term != "const")
             & (R.term != "lnDp") | (R.model.str.startswith("M2_") & R.term.isin(["lnDp"]))]
    prim = prim[~((prim.model.str.startswith("M3_")) & prim.term.isin(["lnDp", "lnsinlat"]))]
    R["primary"] = False; R.loc[prim.index, "primary"] = True
    R["q_primary"] = np.nan
    R.loc[R.primary, "q_primary"] = bh(R.loc[R.primary, "p"].values)
    ok = R.p.notna() & (R.term != "const")
    R["q_all"] = np.nan; R.loc[ok, "q_all"] = bh(R.loc[ok, "p"].values)
    R["mde80"] = 2.8 * R.se
    for c in ("est", "se", "lo", "hi", "mde80"):
        R["eps_" + c] = R[c] * R.fac.fillna(1.0)
    R.to_csv(os.path.join(OUT, "results.csv"), index=False)
    print(R[R.primary][["model", "term", "est", "se", "lo", "hi", "p", "q_primary"]].round(4).to_string())
    print("primary:", int(R.primary.sum()), "tests;", int((R.q_primary < 0.05).sum()), "pass q<0.05")


if __name__ == "__main__":
    main()
