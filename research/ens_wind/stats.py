"""Reliability, calibration slope, Brier skill, leave-one-season-out isotonic calibration and season-block bootstrap, from count tables.
counts: n[k], y[k] for k = members reaching the threshold (k = 0..N). Registered in PREREGISTRATION.md."""
import numpy as np

BINS = [(0, 0), (1e-9, 0.10), (0.10, 0.30), (0.30, 0.60), (0.60, 1.0)]    # fraction bins: 0, (0,10%], (10,30%], (30,60%], (60,100%]
BINNAMES = ["none", "up to 1 in 10", "1 to 3 in 10", "3 to 6 in 10", "over 6 in 10"]
MIN_EVENTS = 30
SEED = 20261028
NBOOT = 2000


def bin_of(N):
    """bin index of each member count k"""
    f = np.arange(N + 1) / N
    b = np.zeros(N + 1, int)
    for i, (lo, hi) in enumerate(BINS):
        if i == 0: b[f == 0] = 0
        else: b[(f > lo) & (f <= hi)] = i
    return b


def reliability(n, y, N):
    b = bin_of(N); nb = np.bincount(b, weights=n, minlength=5); yb = np.bincount(b, weights=y, minlength=5)
    return nb, yb


def brier(n, y, p):
    return float((y * (1 - p) ** 2 + (n - y) * p ** 2).sum())


def pav(p_obs, w):
    """weighted pool-adjacent-violators, non-decreasing fit of p_obs (indexed by k) with weights w"""
    v = list(p_obs.astype(float)); ww = list(w.astype(float)); blocks = [[i] for i in range(len(v))]
    i = 0
    while i < len(v) - 1:
        if v[i] > v[i + 1] + 1e-15:
            tot = ww[i] + ww[i + 1]; v[i] = (v[i] * ww[i] + v[i + 1] * ww[i + 1]) / tot if tot > 0 else v[i]; ww[i] = tot
            blocks[i] += blocks[i + 1]; del v[i + 1], ww[i + 1], blocks[i + 1]; i = max(i - 1, 0)
        else: i += 1
    out = np.zeros(len(p_obs))
    for val, bl in zip(v, blocks): out[bl] = val
    return out


def iso_fit(n, y):
    """calibrated chance for each k from pooled training counts; k with no cases take the neighbour below (monotone fill)"""
    k = np.arange(len(n)); have = n > 0
    p = np.where(have, y / np.maximum(n, 1), 0.0)
    fit = pav(p[have], n[have]); out = np.interp(k, k[have], fit)
    return np.clip(out, 0.0, 1.0)


def slope(n, y, N, it=50):
    """logistic y ~ a + b*logit(clip(k/N, 0.5/N, 1-0.5/N)) over k >= 1 by damped weighted Newton; returns (a, b) or (nan, nan)"""
    k = np.arange(N + 1); f = np.clip(k / N, 0.5 / N, 1 - 0.5 / N); x = np.log(f / (1 - f))
    m = (n > 0) & (k >= 1)       # only cases where at least one member reaches it; the none bin is reported on its own
    if y[m].sum() == 0 or y[m].sum() == n[m].sum(): return np.nan, np.nan
    nn = n[m].astype(float); yy = y[m].astype(float)
    sep = (yy == 0) | (yy == nn); nn = nn + sep; yy = yy + 0.5 * sep          # groups with no events or all events get half an extra case each way, so the fit stays finite
    X = np.column_stack([np.ones(m.sum()), x[m]])
    def dev(be):
        e = np.clip(X @ be, -30, 30); mu = 1 / (1 + np.exp(-e))
        return -2 * float((yy * np.log(mu + 1e-12) + (nn - yy) * np.log(1 - mu + 1e-12)).sum())
    beta = np.array([0.0, 1.0]); d0 = dev(beta)               # start at "the fraction is the chance"; damped Newton
    for _ in range(it):
        eta = X @ beta; mu = 1 / (1 + np.exp(-np.clip(eta, -30, 30))); W = nn * mu * (1 - mu) + 1e-9
        g = X.T @ (yy - nn * mu); H = X.T @ (W[:, None] * X) + 1e-9 * np.eye(2)
        step = np.linalg.solve(H, g); t = 1.0
        while t > 1e-6 and dev(beta + t * step) > d0 + 1e-9: t /= 2
        nb = beta + t * step
        if np.max(np.abs(nb - beta)) < 1e-9: beta = nb; break
        beta = nb; d0 = dev(beta)
    return float(beta[0]), float(beta[1])


def season_tables(n_by_s, y_by_s, N):
    """per held-out season: Brier sums of raw fraction, LOSO isotonic, and training climatology; all in one table."""
    S = sorted(n_by_s)
    raw = np.arange(N + 1) / N; out = {}
    for s in S:
        ntr = sum(n_by_s[t] for t in S if t != s); ytr = sum(y_by_s[t] for t in S if t != s)
        cal = iso_fit(ntr, ytr); clim = ytr.sum() / max(ntr.sum(), 1)
        n, y = n_by_s[s], y_by_s[s]
        out[s] = dict(N=n.sum(), y=y.sum(), b_raw=brier(n, y, raw), b_cal=brier(n, y, cal), b_clim=brier(n, y, np.full(N + 1, clim)))
    return out


def summarize(n_by_s, y_by_s, N, rng=None, nboot=NBOOT):
    """n_by_s, y_by_s: {season: counts array}. Returns dict with pooled values and season-block bootstrap intervals."""
    S = sorted(n_by_s); rng = rng or np.random.default_rng(SEED)
    n = sum(n_by_s[s] for s in S); y = sum(y_by_s[s] for s in S)
    res = dict(N=N, cases=int(n.sum()), events=int(y.sum()), seasons=len(S))
    if res["events"] < MIN_EVENTS:
        res["status"] = "too few events (%d) to calibrate" % res["events"]; return res
    res["status"] = "ok"
    nb, yb = reliability(n, y, N); res["bins_n"] = nb.tolist(); res["bins_y"] = yb.tolist()
    a, b = slope(n, y, N); res["slope"] = b; res["intercept"] = a
    st = season_tables(n_by_s, y_by_s, N); res["by_season"] = {str(s): {k: float(v) for k, v in d.items()} for s, d in st.items()}
    tot = lambda key, idx: sum(st[S[i]][key] for i in idx)
    full = range(len(S))
    res["bss_raw"] = 1 - tot("b_raw", full) / tot("b_clim", full); res["bss_cal"] = 1 - tot("b_cal", full) / tot("b_clim", full)
    res["cal_gain"] = res["bss_cal"] - res["bss_raw"]
    res["cal_table"] = iso_fit(n, y).tolist()
    # bootstrap by season
    slopes = []; gains = []; braw = []; rates = []
    for _ in range(nboot):
        pick = rng.integers(0, len(S), len(S))
        nn = sum(n_by_s[S[i]] for i in pick); yy = sum(y_by_s[S[i]] for i in pick)
        if yy.sum() == 0: continue
        slopes.append(slope(nn, yy, N)[1])
        den = tot("b_clim", pick)
        if den > 0:
            braw.append(1 - tot("b_raw", pick) / den); gains.append((tot("b_raw", pick) - tot("b_cal", pick)) / den)
        bn, by = reliability(nn, yy, N); rates.append(np.where(bn > 0, by / np.maximum(bn, 1), np.nan))
    sl = np.array(slopes); sl = sl[np.isfinite(sl)]
    res["slope_ci"] = np.percentile(sl, [2.5, 97.5]).tolist()
    res["p_overconf"] = float((np.sum(sl >= 1.0) + 1) / (len(sl) + 1))     # one-sided: H1 slope < 1
    res["bss_raw_ci"] = np.percentile(braw, [2.5, 97.5]).tolist(); res["cal_gain_ci"] = np.percentile(gains, [2.5, 97.5]).tolist()
    r = np.array(rates); res["rate_ci"] = [np.nanpercentile(r, [2.5, 97.5], axis=0)[:, i].tolist() for i in range(5)]
    return res


def bh(p, q=0.05):
    """Benjamini-Hochberg: returns list of booleans (rejected)"""
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); thr = q * (np.arange(1, m + 1) / m)
    ok = p[o] <= thr; k = np.max(np.nonzero(ok)[0]) + 1 if ok.any() else 0
    rej = np.zeros(m, bool); rej[o[:k]] = True; return rej.tolist()
