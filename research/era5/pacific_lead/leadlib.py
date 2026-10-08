"""Shared pieces for RA-5 (Pacific lead on the Atlantic). ERA5 PROXY; pipeline A only in S3/S4. Plan: PREREGISTRATION.md.

Weekly design of PR 41/64: 30 weeks of 7 days from 1 October, seasons labelled by starting year. A Dataset holds, for an
outcome basin O and a predictor basin Q, arrays of shape (seasons, 30). Tests use weeks 2..29 so a lag of two weeks stays in
season. The Poisson fit is the plain IRLS of hem_channels/chanlib.
"""
import datetime as dt
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ERA, "hem_channels"))
import chanlib as C  # noqa: E402

NW = 30
MONTHS = C.MONTHS
SEED = C.SEED + 5
HS = os.path.join(ERA, "hemispheric", "results")
WT = pd.read_csv(C.WEEKLY)
NAMED = ["nao", "pna", "oni", "mjo1", "mjo2"]
BH = C.bh


def zs(x):
    """Standardise; NaN (two MJO weeks in 2004-2025, logged as a deviation) become 0, the mean."""
    x = np.asarray(x, float)
    return np.nan_to_num((x - np.nanmean(x)) / np.nanstd(x), nan=0.0)


def month_dummies(seasons):
    mid = [(C.week_start(s, k) + dt.timedelta(days=3)).month for s in seasons for k in range(NW)]
    m = np.array([[1.0 if v == mm else 0.0 for mm in MONTHS] for v in mid])
    return m.reshape(len(seasons), NW, len(MONTHS))


class Dataset:
    """Arrays (ns, 30): idx[b], cnt[b], prev[b] for b in 'atl','pac'; M (ns,30,7); named (ns,30,5)."""

    def __init__(self, seasons, idx, cnt, prev, named, label):
        self.seasons, self.idx, self.cnt, self.prev, self.named, self.label = seasons, idx, cnt, prev, named, label
        self.M = month_dummies(seasons)
        self.ns = len(seasons)


def _rs(v, ns):
    return np.asarray(v, float).reshape(ns, NW)


def load_2004(counts="archive"):
    """22 seasons, out-of-sample (leave-one-season-out) indices of PR 41; archive counts (primary) or pipeline A HF-equivalent."""
    S = list(range(2004, 2026))
    ns = len(S)
    idx, cnt, prev = {}, {}, {}
    W = WT[WT.season >= 2004].sort_values(["season", "week"]).reset_index(drop=True)
    for b in ("atl", "pac"):
        o = pd.read_csv(os.path.join(HS, f"oos_index_{b}.csv"))
        assert (o.season.values == np.repeat(S, NW)).all()
        assert (W[f"y_{b}"].values == o.y.values).all()
        idx[b] = _rs(zs(o.idx.values), ns)
        if counts == "archive":
            cnt[b] = _rs(o.y.values, ns); prev[b] = _rs(W[f"prev_{b}"].values, ns)
        elif counts == "proxy":
            cnt[b] = _rs(W[f"pAhf_{b}"].values, ns); prev[b] = _rs(W[f"pAhfprev_{b}"].values, ns)
        else:
            raise ValueError(counts)
    named = np.stack([_rs(zs(W[n].values), ns) for n in NAMED], axis=-1)
    return Dataset(S, idx, cnt, prev, named, f"2004-2025 {counts}")


def load_1979_2000():
    """Frozen primary-split patterns (fitted 2004-14) applied to 1979-2000 PC scores; pipeline A depth counts (within-era)."""
    S = list(range(1979, 2001))
    ns = len(S)
    fz = json.load(open(os.path.join(HS, "frozen_primary.json")))
    W = WT[(WT.season >= 1979) & (WT.season <= 2000)].sort_values(["season", "week"]).reset_index(drop=True)
    cols = [f"{v}{i}" for v in "zus" for i in range(1, 11)]
    idx, cnt, prev = {}, {}, {}
    for b in ("atl", "pac"):
        beta = np.array(fz["basins"][b]["models"]["P"]["beta"])
        idx[b] = _rs(zs(W[cols].values @ beta[8:]), ns)
        cnt[b] = _rs(W[f"pAdepth_{b}"].values, ns)
        prev[b] = _rs(W[f"pAdepthprev_{b}"].values, ns)
    named = np.stack([_rs(zs(W[n].values), ns) for n in NAMED], axis=-1)
    return Dataset(S, idx, cnt, prev, named, "1979-2000 pipeline A depth (frozen patterns)")


# ------------------------------------------------------------------ one test
class Spec:
    def __init__(self, name, kind, k, pred, extras=False, O="atl", Q="pac"):
        self.name, self.kind, self.k, self.pred, self.extras, self.O, self.Q = name, kind, k, pred, extras, O, Q


def design(ds, sp, seas=None, perm=None):
    """y, X (tested term last), n_obs. perm: season permutation applied to the predictor basin's series."""
    seas = np.arange(ds.ns) if seas is None else np.asarray(seas)
    O, Q, k = sp.O, sp.Q, sp.k
    q_idx = ds.idx[Q] if perm is None else ds.idx[Q][perm]
    q_cnt = ds.cnt[Q] if perm is None else ds.cnt[Q][perm]
    if sp.pred == "idx":
        xq = q_idx
    else:  # standardised log(1+count) over all seasons of this dataset
        xq = zs(np.log1p(ds.cnt[Q])) if perm is None else zs(np.log1p(ds.cnt[Q]))[perm]
    sl, sk = slice(2, NW), slice(2 - k, NW - k)
    M = ds.M[seas][:, sl].reshape(-1, ds.M.shape[2])
    cols = [M]
    if sp.kind == "pois":
        y = ds.cnt[O][seas][:, sl].ravel()
        cols += [np.log1p(ds.prev[O][seas][:, sl]).reshape(-1, 1), ds.idx[O][seas][:, sl].reshape(-1, 1)]
    else:
        y = ds.idx[O][seas][:, sl].ravel()
        cols += [ds.idx[O][seas][:, 1:NW - 1].reshape(-1, 1)]
    if sp.extras:
        cols += [ds.named[seas][:, sl].reshape(-1, len(NAMED)), ds.named[seas][:, sk].reshape(-1, len(NAMED))]
    cols.append(xq[seas][:, sk].reshape(-1, 1))
    return y, np.column_stack(cols)


def fit(y, X, kind):
    """(slope, partial r or nan)."""
    if kind == "pois":
        return C.pois(X, y)[-1], np.nan
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    df = len(y) - X.shape[1]
    s2 = res @ res / df
    xtx_i = np.linalg.pinv(X.T @ X)
    t = b[-1] / np.sqrt(s2 * xtx_i[-1, -1])
    return b[-1], t / np.sqrt(t * t + df)


def stat(ds, sp, seas=None, perm=None):
    y, X = design(ds, sp, seas, perm)
    return fit(y, X, sp.kind)


def perm_test(ds, sp, perms):
    obs, r = stat(ds, sp)
    nul = np.array([stat(ds, sp, perm=p)[0] for p in perms])
    return obs, r, (1 + np.sum(nul >= obs)) / (1 + len(perms)), (1 + np.sum(nul <= obs)) / (1 + len(perms)), nul


def boot(ds, sp, rng, nb, seas_pool=None):
    pool = np.arange(ds.ns) if seas_pool is None else np.asarray(seas_pool)
    out = np.empty((nb, 2))
    for i in range(nb):
        s = rng.choice(pool, len(pool))
        out[i] = stat(ds, sp, seas=s)
    return out
