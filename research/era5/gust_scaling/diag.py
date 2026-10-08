"""Tail-model check (post hoc, logged in PREREGISTRATION.md): which tail describes the HF-fix gusts, truncated normal on
ln g800 or generalised Pareto on g800 - 71.7? Probability-integral-transform KS distance, and observed vs fitted exceedance
quantiles. Uses the M2 design of scaling.py, point estimates only."""
import numpy as np, pandas as pd, json, os
from scipy.special import ndtr
import scaling as S

f = S.load(); f = f[f.Dp >= 5].copy(); rows = []
for _b in ("atl", "pac"):
    _m = f.basin == _b
    for _k in ("oni", "pna"):
        f.loc[_m, _k] = (f.loc[_m, _k + "_raw"] - f.loc[_m, _k + "_raw"].mean()) / f.loc[_m, _k + "_raw"].std()
S.trunc_diag(f)  # six Nelder-Mead starts: does the pre-registered MLE have one maximum?
# profile log-likelihood over sigma (location problem solved at each sigma)
prof = []
for _b in ("atl", "pac"):
    _d = f[f.basin == _b]; _X = S.design(_d).values; _y = _d.y.values; _b0 = np.linalg.lstsq(_X, _y, rcond=None)[0]
    for _sg in (0.07, 0.09, 0.11, 0.13, 0.16, 0.20, 0.30, 0.50):
        _beta, _ll = S._inner(_X, _y, _sg, _b0)
        prof.append(dict(basin=_b, sigma=_sg, loglik=_ll, lnDp=_beta[1], lnR=_beta[2], lnsinlat=_beta[3]))
pd.DataFrame(prof).to_csv(os.path.join(S.OUT, "profile_sigma.csv"), index=False)
for b in ("atl", "pac"):
    d = f[f.basin == b].copy()
    for k in ("oni", "pna"):
        d[k] = (d[k + "_raw"] - d[k + "_raw"].mean()) / d[k + "_raw"].std()
    X = S.design(d).values; y = d.y.values; e = np.maximum(d.g800.values - 71.7, 0)
    th = S.trunc_fit(X, y); mu = X @ th[:-1]; sg = np.exp(th[-1]); a = (S.CUT - mu) / sg
    u_n = (ndtr((y - mu) / sg) - ndtr(a)) / (1 - ndtr(a))
    tg = S.gpd_fit(X, e); sgm = np.exp(X @ tg[:-1]); xi = tg[-1]
    u_g = 1 - (1 + xi * e / sgm) ** (-1 / xi)

    def ks(u):
        u = np.sort(u); n = len(u); return float(max(np.max(np.arange(1, n + 1) / n - u), np.max(u - np.arange(0, n) / n)))
    # simulated fitted quantiles of g800 for the first covariate row mix: use inverse CDF per fix, mean over fixes
    qs = [0.5, 0.9, 0.99]
    obs = np.quantile(d.g800, qs)
    rng = np.random.default_rng(0)
    un = rng.uniform(size=(200, len(d)))
    gn = np.exp(mu + sg * __import__("scipy").stats.norm.ppf(ndtr(a) + un * (1 - ndtr(a))))
    gg = 71.7 + sgm * ((1 - un) ** (-xi) - 1) / xi
    rows.append(dict(basin=b, n=len(d), ks_trunc_normal=ks(u_n), ks_gpd=ks(u_g), xi=float(xi), sigma_ln=float(sg),
                     obs_q50=obs[0], obs_q90=obs[1], obs_q99=obs[2],
                     tn_q50=float(np.quantile(gn, .5)), tn_q90=float(np.quantile(gn, .9)), tn_q99=float(np.quantile(gn, .99)),
                     gpd_q50=float(np.quantile(gg, .5)), gpd_q90=float(np.quantile(gg, .9)), gpd_q99=float(np.quantile(gg, .99)),
                     max_obs=float(d.g800.max()), ks_5pct_crit=float(1.36 / np.sqrt(len(d)))))
R = pd.DataFrame(rows); R.to_csv(os.path.join(S.OUT, "tail_diagnostic.csv"), index=False)
print(R.round(4).T.to_string())
