"""POST HOC follow-ups on the pre-registered results (written after seeing them; not in PREREGISTRATION.md).

Every number from this script is labelled post hoc in README.md and counts as exploratory.
  a. lower tails: W7, W30 and the 48 h gap count against the same N1 simulations, tail toward regularity
  b. Knox excess by time separation (pairs within 1000 km), season-month permutation null as pre-registered
  c. Knox conditioned on the circulation: times permuted only among events in the same season-month AND
     the same tercile of the basin index (NAO for the Atlantic, PNA for the Pacific), at the lagged value
     (days -10..-4, as the main conditioning) and at the same-time value (days -3..+3, circular: upper bound)
  d. per-season Knox ratio for pairs within 72 h and 1000 km, with a season-block bootstrap interval

usage: posthoc.py TIER OUT_DIR [--nperm N] [--nsim N]
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

import cluster as C

UP = np.array([24, 48, 72, 96, 168, 336, 720])  # right-inclusive upper edges, hours
NB = len(UP)


def pair_table(ev, dmax=1000.0):
    season = np.array([e["season"] for e in ev])
    lat = np.array([e["lat"] for e in ev])
    lon = np.array([e["lon"] for e in ev])
    ii, jj = [], []
    for s in np.unique(season):
        ix = np.where(season == s)[0]
        a, b = np.triu_indices(len(ix), 1)
        D = C.hav(lat[ix[a]], lon[ix[a]], lat[ix[b]], lon[ix[b]])
        k = D <= dmax
        ii.append(ix[a][k])
        jj.append(ix[b][k])
    return np.concatenate(ii), np.concatenate(jj), season


def knox_bins(ev, group, rng, nperm):
    """Per time-separation bin and season: observed pair counts and permutation mean; pooled perm totals."""
    n = len(ev)
    t = np.array([e["t"] for e in ev])
    ii, jj, season = pair_table(ev)
    sp = season[ii] - season.min()
    ns = season.max() - season.min() + 1
    order = np.argsort(group, kind="stable")
    t_sorted = t[order]

    def counts(tt):
        dt_ = np.abs(tt[ii] - tt[jj])
        b = np.searchsorted(UP, dt_, side="left")
        ok = b < NB
        return np.bincount(b[ok] * ns + sp[ok], minlength=NB * ns).reshape(NB, ns)

    obs = counts(t)
    acc = np.zeros_like(obs, float)
    tot = np.empty((nperm, NB))
    for k in range(nperm):
        idx = np.lexsort((rng.random(n), group))
        tp = np.empty(n)
        tp[idx] = t_sorted
        c = counts(tp)
        acc += c
        tot[k] = c.sum(1)
    return obs, acc / nperm, tot


def summarise(obs, mean, tot, rng, nboot=2000):
    nperm = len(tot)
    o = obs.sum(1)
    m = mean.sum(1)
    p = (1 + (tot >= o).sum(0)) / (1 + nperm)
    out = dict(ratio_by_bin=(o / m).tolist(), p_by_bin=p.tolist(), obs_by_bin=o.tolist(), null_by_bin=m.tolist())
    # pooled <= 72 h (bins 0-2), season-level
    so, sm_ = obs[:3].sum(0), mean[:3].sum(0)
    ok = sm_ > 0
    out["ratio_72h"] = float(so.sum() / sm_.sum())
    out["p_72h"] = float((1 + (tot[:, :3].sum(1) >= so.sum()).sum()) / (1 + nperm))
    out["seasons_above_1"] = int((so[ok] > sm_[ok]).sum())
    out["seasons"] = int(ok.sum())
    ns = len(so)
    bs = []
    for _ in range(nboot):
        ix = rng.integers(ns, size=ns)
        bs.append(so[ix].sum() / sm_[ix].sum())
    out["ratio_72h_ci"] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
    return out


def same_time_index(cpc, s0, s1):
    out = {}
    for k in ("NAO", "PNA"):
        s = C.read_cpc(os.path.join(cpc, f"norm.daily.{k.lower()}.index.b500101.current.ascii"))
        c = s.rolling(7, min_periods=7, center=True).mean()
        out[k] = np.array([[c.get(pd.Timestamp(s_, 10, 1) + pd.Timedelta(days=d), np.nan) for d in range(C.NDAY)]
                           for s_ in range(s0, s1 + 1)])
    return out


def tercile_group(ev, table, s0):
    v = np.array([table[e["season"] - s0, int(e["t"] // 24)] for e in ev])
    v = np.where(np.isnan(v), np.nanmedian(v), v)
    q = np.quantile(v, [1 / 3, 2 / 3])
    terc = np.digitize(v, q)
    season = np.array([e["season"] for e in ev])
    t = np.array([e["t"] for e in ev])
    mon = C.MONIDX[np.clip((t // 24).astype(int), 0, C.NDAY - 1)]
    return ((season - season.min()) * 7 + mon) * 3 + terc


def tails(ev, s0, s1, rng, nsim, era):
    S = s1 - s0 + 1
    tr = np.arange(S) - (S - 1) / 2
    y = C.daily(ev, s0, s1)
    X = C.make_X("N1", S, tr, era)
    D = C.dispersion(y, X, (7, 30), rng, nsim)
    # lower-tail p needs the simulated phi again: redo cheaply with the same recipe
    mu = D["mu"]
    sims = np.empty((nsim, 2))
    for i in range(nsim):
        ys = rng.poisson(mu)
        _, ms = C.fit_pois(ys.ravel().astype(float), X)
        sims[i] = C.phis(ys, ms.reshape(S, C.NDAY), X.shape[1], (7, 30))
    obs = D["phi"]
    low = ((1 + (sims <= obs).sum(0)) / (1 + nsim)).tolist()
    r = C.m2(ev, mu, s0, rng, nsim)
    # lower tail for G48 needs the simulated counts; recompute
    lam = C.make_lam(mu, s0)
    g = np.empty(nsim)
    sidx = np.repeat(np.arange(S), C.NDAY)
    didx = np.tile(np.arange(C.NDAY), S)
    for i in range(nsim):
        c = rng.poisson(mu).ravel()
        tt = (np.repeat(didx, c) * 24 + 6 * rng.integers(0, 4, size=c.sum())).astype(float)
        g[i] = C.gap_stats(np.repeat(sidx, c) + s0, tt, lam)["G48"]
    gl = float((1 + (g <= r["obs"]["G48"]).sum()) / (1 + nsim))
    return dict(E_W7=float(D["E"][0]), E_W30=float(D["E"][1]), p_lower_W7=low[0], p_lower_W30=low[1],
                G48=r["obs"]["G48"], G48_null=float(g.mean()), p_lower_G48=gl)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tier")
    ap.add_argument("out")
    ap.add_argument("--nperm", type=int, default=4000)
    ap.add_argument("--nsim", type=int, default=2000)
    ap.add_argument("--cpc", default=C.CPC_DEFAULT)
    a = ap.parse_args()
    rows, cfg = C.tier_data(a.tier)
    s0, s1 = cfg["s0"], cfg["s1"]
    rng = np.random.default_rng(C.SEED + 31 + sum(map(ord, a.tier)))
    lag, _ = C.index_tables(a.cpc, s0, s1)
    now = same_time_index(a.cpc, s0, s1)
    era = (np.arange(s0, s1 + 1) >= 2004).astype(float) if a.tier == "r3" else None
    res = {}
    for basin in ("atl", "pac"):
        ev, _ = C.window(rows, basin, s0, s1)
        ib = "NAO" if basin == "atl" else "PNA"
        t = np.array([e["t"] for e in ev])
        season = np.array([e["season"] for e in ev])
        mon = C.MONIDX[np.clip((t // 24).astype(int), 0, C.NDAY - 1)]
        base = (season - season.min()) * 7 + mon
        R = {}
        R["tails"] = tails(ev, s0, s1, rng, a.nsim, era)
        for name, grp in (("base", base), ("lag_tercile", tercile_group(ev, lag[ib], s0)),
                          ("now_tercile", tercile_group(ev, now[ib], s0))):
            o, m, tot = knox_bins(ev, grp, rng, a.nperm)
            R["knox_" + name] = summarise(o, m, tot, rng)
        res[basin] = R
        print(basin, json.dumps({k: (v if k == "tails" else dict(ratio_72h=v["ratio_72h"], p=v["p_72h"], ci=v["ratio_72h_ci"],
              ssn=f'{v["seasons_above_1"]}/{v["seasons"]}')) for k, v in R.items()}), flush=True)
    os.makedirs(a.out, exist_ok=True)
    json.dump(dict(tier=a.tier, nperm=a.nperm, nsim=a.nsim, bins_upper_edges_h=UP.tolist(), basins=res),
              open(os.path.join(a.out, f"{a.tier}_posthoc.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
