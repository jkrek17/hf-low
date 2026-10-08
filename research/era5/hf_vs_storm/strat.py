"""Stratified HF - SF difference with a season-block bootstrap (hf_vs_storm).

cells(vals, grp, season, stratum) accumulates per (season, stratum, group) sums and non-missing counts.
estimate / boot give diff = sum_m w_m (mean_HF,m - mean_SF,m) / sum_m w_m with w_m the HF count in stratum m (months or depth bins
where both groups are present), and resample the seasons with replacement. grp: 0 = HF, 1 = SF.
"""
import numpy as np


def cells(vals, grp, season, stratum, seasons, nstr):
    S = len(seasons)
    si = np.searchsorted(seasons, season)
    shape = (S, nstr, 2) + vals.shape[1:]
    sums = np.zeros(shape, np.float64)
    cnt = np.zeros(shape, np.float64)
    ok = ~np.isnan(vals)
    v0 = np.where(ok, vals, 0.0)
    for g in (0, 1):
        sel = grp == g
        np.add.at(sums[:, :, g], (si[sel], stratum[sel]), v0[sel])
        np.add.at(cnt[:, :, g], (si[sel], stratum[sel]), ok[sel].astype(float))
    return sums, cnt


def estimate(sums, cnt, c):
    """c: season multiplicities [S]. Returns diff, HF mean, SF mean (weighted to HF strata)."""
    ss = np.tensordot(c, sums, axes=(0, 0))          # [M, 2, ...]
    sc = np.tensordot(c, cnt, axes=(0, 0))
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = ss / sc
    valid = (sc[:, 0] > 0) & (sc[:, 1] > 0)
    w = np.where(valid, sc[:, 0], 0.0)
    wsum = w.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        hf = (np.where(valid, mu[:, 0], 0) * w).sum(0) / wsum
        sf = (np.where(valid, mu[:, 1], 0) * w).sum(0) / wsum
    out = hf - sf
    bad = wsum <= 0
    return np.where(bad, np.nan, out), np.where(bad, np.nan, hf), np.where(bad, np.nan, sf)


def boot(sums, cnt, B, rng, keep=True):
    """Returns observed (diff, hf, sf), bootstrap diffs [B, ...] (if keep) and sign counts."""
    S = sums.shape[0]
    obs = estimate(sums, cnt, np.ones(S))
    pos = np.zeros(obs[0].shape); neg = np.zeros(obs[0].shape); nn = np.zeros(obs[0].shape)
    reps = np.empty((B,) + obs[0].shape, np.float32) if keep else None
    for b in range(B):
        c = np.bincount(rng.integers(0, S, S), minlength=S).astype(float)
        d = estimate(sums, cnt, c)[0]
        if keep:
            reps[b] = d
        pos += d > 0; neg += d < 0; nn += ~np.isnan(d)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = 2 * np.minimum((pos + 1) / (nn + 1), (neg + 1) / (nn + 1))
    p = np.minimum(p, 1.0)
    p = np.where(nn < B * 0.9, np.nan, p)           # unusable where too many resamples lack a stratum
    return obs, reps, p


def bh(p, q=0.05):
    """Benjamini-Hochberg on a flat array with NaN ignored. Returns (qvals, reject)."""
    p = np.asarray(p, float)
    flat = p.ravel()
    ok = ~np.isnan(flat)
    qv = np.full(flat.shape, np.nan)
    pv = flat[ok]
    n = len(pv)
    if n == 0:
        return qv.reshape(p.shape), np.zeros(p.shape, bool)
    o = np.argsort(pv)
    adj = pv[o] * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    r = np.empty(n); r[o] = np.minimum(adj, 1.0)
    qv[ok] = r
    return qv.reshape(p.shape), (qv <= q).reshape(p.shape)
