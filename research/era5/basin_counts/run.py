"""RA-25 stage 1: are weekly HF counts anticorrelated between basins? ERA5 PROXY (pipeline A) plus the archive.
Plan: PREREGISTRATION.md (commit 3b7a49e, before any count correlation). Writes results/*.csv and results/summary.json.
Run: python3 research/era5/basin_counts/run.py
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
HS = os.path.join(ERA, "hemispheric", "results")
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
NW = 30
SEED = 20261008 + 25
NPERM, NBOOT = 20000, 5000
WT = pd.read_csv(os.path.join(HS, "weekly_table.csv.gz")).sort_values(["season", "week"]).reset_index(drop=True)
MONTHS = [10, 11, 12, 1, 2, 3, 4]

SAMPLES = {
    "ARCH": ("archive", "y_{}", range(2004, 2026)),
    "PRX": ("pipeline A gust proxy", "pAhf_{}", range(2004, 2026)),
    "DEP": ("pipeline A depth proxy", "pAdepth_{}", range(1979, 2001)),
    "DEP_fit": ("pipeline A depth proxy", "pAdepth_{}", range(2004, 2026)),
}


def block(sample, basin):
    _, pat, ys = SAMPLES[sample]
    w = WT[WT.season.isin(list(ys))]
    assert len(w) == 22 * NW and w[pat.format(basin)].notna().all(), (sample, basin)
    return w[pat.format(basin)].values.astype(float).reshape(22, NW), w


def month_dummies(seasons):
    import datetime as dt
    rows = []
    for s in seasons:
        for k in range(NW):
            m = (dt.date(s, 10, 1) + dt.timedelta(days=7 * k + 3)).month
            rows.append([1.0 if m == mm else 0.0 for mm in MONTHS])
    return np.array(rows)


def deseason_week(x):
    return x - x.mean(axis=0, keepdims=True)


def corr(a, b):
    a = a.ravel() - a.mean(); b = b.ravel() - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def perm_stats(a, b, rng, nperm=NPERM):
    """r of a against b with whole seasons of b shuffled among seasons; a, b shape (ns, 30)."""
    ns = a.shape[0]
    av = a.ravel() - a.mean()
    na = np.sqrt((av * av).sum())
    out = np.empty(nperm)
    for i in range(nperm):
        bp = b[rng.permutation(ns)].ravel()
        bp = bp - bp.mean()
        out[i] = (av * bp).sum() / (na * np.sqrt((bp * bp).sum()))
    return out


def boot_ci(a, b, rng, nboot=NBOOT, fn=corr):
    ns = a.shape[0]
    v = np.empty(nboot)
    for i in range(nboot):
        ix = rng.integers(0, ns, ns)
        v[i] = fn(a[ix], b[ix])
    return np.percentile(v, [2.5, 97.5])


def test(a, b, rng, name, label, extra=None):
    r = corr(a, b)
    null = perm_stats(a, b, rng)
    p_neg = (1 + (null <= r).sum()) / (1 + len(null))
    p_two = (1 + (np.abs(null) >= abs(r)).sum()) / (1 + len(null))
    lo, hi = boot_ci(a, b, rng)
    env = np.percentile(null, [2.5, 97.5])
    d = dict(test=name, sample=label, r=r, ci_lo=lo, ci_hi=hi, p_neg=p_neg, p_two=p_two, null_lo=env[0], null_hi=env[1])
    if extra:
        d.update(extra)
    return d


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); m = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        m = min(m, p[i] * n / rank); q[i] = m
    return q


def resid_on(y, X):
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    summ = {}
    D = {}
    for s in SAMPLES:
        A, w = block(s, "atl"); P, _ = block(s, "pac")
        D[s] = dict(A=A, P=P, w=w, seasons=sorted(set(w.season)))

    # --- primary
    for s in ("ARCH", "PRX", "DEP"):
        a, b = deseason_week(D[s]["A"]), deseason_week(D[s]["P"])
        D[s]["a"], D[s]["b"] = a, b
        rows.append(dict(family="primary", **test(a, b, rng, "primary", s)))
    # --- S1 month dummies
    for s in ("ARCH", "PRX", "DEP"):
        M = month_dummies(D[s]["seasons"])
        a = resid_on(D[s]["A"].ravel(), M).reshape(22, NW); b = resid_on(D[s]["P"].ravel(), M).reshape(22, NW)
        rows.append(dict(family="S1", **test(a, b, rng, "S1_month", s)))
    # --- S2a interannual: season totals; S2b subseasonal
    for s in ("ARCH", "PRX", "DEP"):
        ta, tb = D[s]["A"].sum(1), D[s]["P"].sum(1)
        r = float(np.corrcoef(ta, tb)[0, 1])
        null = np.array([np.corrcoef(ta, tb[rng.permutation(22)])[0, 1] for _ in range(NPERM)])
        bs = []
        for _ in range(NBOOT):
            ix = rng.integers(0, 22, 22)
            if ta[ix].std() > 0 and tb[ix].std() > 0:
                bs.append(np.corrcoef(ta[ix], tb[ix])[0, 1])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        rows.append(dict(family="S2a", test="S2a_interannual", sample=s, r=r, ci_lo=lo, ci_hi=hi,
                         p_neg=(1 + (null <= r).sum()) / (1 + NPERM), p_two=(1 + (np.abs(null) >= abs(r)).sum()) / (1 + NPERM),
                         null_lo=np.percentile(null, 2.5), null_hi=np.percentile(null, 97.5)))
        a = D[s]["a"] - D[s]["a"].mean(1, keepdims=True); b = D[s]["b"] - D[s]["b"].mean(1, keepdims=True)
        rows.append(dict(family="S2b", **test(a, b, rng, "S2b_subseasonal", s)))
    # --- S3 named indices removed
    for s in ("ARCH", "PRX", "DEP"):
        w = D[s]["w"]
        X = np.column_stack([w.nao.values, w.pna.values])
        a = resid_on(D[s]["a"].ravel(), X).reshape(22, NW); b = resid_on(D[s]["b"].ravel(), X).reshape(22, NW)
        rows.append(dict(family="S3", **test(a, b, rng, "S3_minus_NAO_PNA", s)))
    # --- S4 depth in the fit era
    a, b = deseason_week(D["DEP_fit"]["A"]), deseason_week(D["DEP_fit"]["P"])
    rows.append(dict(family="S4", **test(a, b, rng, "S4_depth_2004_2025", "DEP_fit")))

    T = pd.DataFrame(rows)
    T["q_fam1"] = np.nan
    m = T.family == "primary"
    T.loc[m, "q_fam1"] = bh(T.loc[m, "p_neg"])
    T["q_all"] = bh(T.p_neg.values)
    T.to_csv(os.path.join(OUT, "tests.csv"), index=False, float_format="%.5g")

    # --- S5 lags (descriptive)
    lag_rows = []
    for s in ("ARCH", "PRX", "DEP"):
        a, b = D[s]["a"], D[s]["b"]
        for k in range(-2, 3):
            if k >= 0:
                x, y = a[:, : NW - k], b[:, k:]
            else:
                x, y = a[:, -k:], b[:, : NW + k]
            r = corr(x, y)
            null = []
            for _ in range(2000):
                null.append(corr(x, y[rng.permutation(22)]))
            lag_rows.append(dict(sample=s, k=k, r=r, env_lo=np.percentile(null, 2.5), env_hi=np.percentile(null, 97.5)))
    pd.DataFrame(lag_rows).to_csv(os.path.join(OUT, "lag_profile.csv"), index=False, float_format="%.4g")

    # --- S6 halves
    hr = []
    for s in ("ARCH", "PRX"):
        for lab, sl in (("2004-2014", slice(0, 11)), ("2015-2025", slice(11, 22))):
            hr.append(dict(sample=s, half=lab, r=corr(D[s]["a"][sl], D[s]["b"][sl])))
    pd.DataFrame(hr).to_csv(os.path.join(OUT, "halves.csv"), index=False, float_format="%.4g")

    # --- S7 coverage
    cr = []
    for s in ("ARCH", "PRX", "DEP"):
        A, P = D[s]["A"], D[s]["P"]
        both0 = (A == 0) & (P == 0)
        keep = ~both0
        a, b = D[s]["a"], D[s]["b"]
        r2 = float(np.corrcoef(a[keep], b[keep])[0, 1])
        cr.append(dict(sample=s, share_both_zero=both0.mean(), share_atl_zero=(A == 0).mean(), share_pac_zero=(P == 0).mean(),
                       mean_atl=A.mean(), mean_pac=P.mean(), var_over_mean_atl=A.var() / A.mean(), var_over_mean_pac=P.var() / P.mean(),
                       r_drop_both_zero=r2))
    pd.DataFrame(cr).to_csv(os.path.join(OUT, "coverage.csv"), index=False, float_format="%.4g")

    # --- S8 index seesaw over the same 660 weeks (ARCH sample)
    ia = pd.read_csv(os.path.join(HS, "oos_index_atl.csv")).idx.values.reshape(22, NW)
    ip = pd.read_csv(os.path.join(HS, "oos_index_pac.csv")).idx.values.reshape(22, NW)
    r_idx = corr(ia, ip)
    ratio = []
    for _ in range(NBOOT):
        ix = rng.integers(0, 22, 22)
        ri = corr(ia[ix], ip[ix])
        ratio.append(corr(D["ARCH"]["a"][ix], D["ARCH"]["b"][ix]) / ri)
    # index r after the same week-of-season deseasonalisation (the count r is deseasonalised)
    r_idx_ds = corr(deseason_week(ia), deseason_week(ip))
    s8 = dict(r_index=r_idx, r_index_deseason=r_idx_ds, ratio_counts_over_index=corr(D["ARCH"]["a"], D["ARCH"]["b"]) / r_idx,
              ratio_ci=list(np.percentile(ratio, [2.5, 97.5])))

    # --- power (planted effect) per sample
    prng = np.random.default_rng(SEED + 1)
    grid = [-0.05, -0.10, -0.15, -0.20, -0.25]
    prow = []
    for s in ("ARCH", "PRX", "DEP"):
        a, b = D[s]["a"], D[s]["b"]
        Za = (a - a.mean()) / a.std(); Zb = (b - b.mean()) / b.std()
        for rho in grid:
            hit1 = hit2 = 0
            nsim = 500
            for _ in range(nsim):
                Zp = Zb[prng.permutation(22)]
                Ps = rho * Za + np.sqrt(1 - rho ** 2) * Zp
                r = corr(Za, Ps)
                null = perm_stats(Za, Ps, prng, nperm=1000)
                p = (1 + (null <= r).sum()) / 1001
                hit1 += p < 0.05; hit2 += p < 0.05 / 3
            prow.append(dict(sample=s, rho=rho, power_p05=hit1 / nsim, power_bonf=hit2 / nsim))
    PW = pd.DataFrame(prow)
    PW.to_csv(os.path.join(OUT, "power.csv"), index=False, float_format="%.3f")
    mde = {}
    for s in ("ARCH", "PRX", "DEP"):
        g = PW[PW["sample"] == s].sort_values("rho", ascending=False)
        pw = g.power_p05.values; rh = g.rho.values  # rho from -0.05 down to -0.25, power rising
        mde[s] = float(np.interp(0.8, pw, rh)) if pw.max() >= 0.8 else None
    summ.update(s8=s8, mde_r_80pct_p05=mde)
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)
    print(T.round(4).to_string())
    print(pd.DataFrame(lag_rows).round(3).to_string())
    print(pd.DataFrame(hr).round(3)); print(pd.DataFrame(cr).round(3).to_string())
    print(json.dumps(summ, indent=1, default=float)); print(PW.round(3).to_string())


if __name__ == "__main__":
    main()
