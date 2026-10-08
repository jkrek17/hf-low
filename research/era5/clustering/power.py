"""Power of the clustering tests against a known burst process (see PREREGISTRATION.md, "Power").

A fraction e of events are secondary lows: each parent low (positions resampled from the real
events of the basin and month) has a child with probability q = e / (1 - e) arriving after an
exponential delay (mean 1.5 days, 6-hourly) and displaced by a normal 500 km. Parent intensity is
mu / (1 + q), so the mean count stays on the fitted N1 intensity. The only input taken from the
data is that fitted seasonal-cycle intensity and the set of event positions, never a test outcome.
e = 0 is the size of each test (a calibration check). Un-adjusted alpha = 0.05, one-sided.

usage: power.py TIER OUT_DIR [--nrep N]
"""
import argparse
import json
import os

import numpy as np

import cluster as C

ES = (0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40)


def draw(mu, pool, s0, e, rng):
    q = e / (1 - e)
    S = mu.shape[0]
    lam = mu / (1 + q)
    rows = []
    c = rng.poisson(lam)
    for si in range(S):
        for d in np.nonzero(c[si])[0]:
            for _ in range(c[si, d]):
                t = d * 24 + 6 * rng.integers(0, 4)
                src = pool[rng.integers(len(pool))]
                rows.append((s0 + si, t, src["lat"], src["lon"]))
                if rng.random() < q:
                    dly = 6 * np.round(rng.exponential(36) / 6)
                    tc = t + dly
                    if tc < C.HRS:
                        la = src["lat"] + rng.normal(0, 500 / 111.2)
                        lo = src["lon"] + rng.normal(0, 500 / (111.2 * max(np.cos(np.radians(src["lat"])), 0.2)))
                        rows.append((s0 + si, tc, float(np.clip(la, -89, 89)), lo % 360))
    return [dict(season=r[0], t=r[1], lat=r[2], lon=r[3]) for r in rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tier")
    ap.add_argument("out")
    ap.add_argument("--nrep", type=int, default=150)
    ap.add_argument("--nnull", type=int, default=1000)
    ap.add_argument("--nperm", type=int, default=199)
    a = ap.parse_args()
    rows, cfg = C.tier_data(a.tier)
    s0, s1 = cfg["s0"], cfg["s1"]
    S = s1 - s0 + 1
    tr = np.arange(S) - (S - 1) / 2
    era = (np.arange(s0, s1 + 1) >= 2004).astype(float) if a.tier == "r3" else None
    rng = np.random.default_rng(C.SEED + 99)
    out = {}
    for basin in ("atl", "pac"):
        ev, _ = C.window(rows, basin, s0, s1)
        y = C.daily(ev, s0, s1)
        X = C.make_X("N1", S, tr, era)
        _, mu = C.fit_pois(y.ravel().astype(float), X)
        mu = mu.reshape(S, C.NDAY)
        p = X.shape[1]
        # null reference from the fitted model
        n7 = np.empty(a.nnull)
        n30 = np.empty(a.nnull)
        g48 = np.empty(a.nnull)
        lam = C.make_lam(mu, s0)
        for i in range(a.nnull):
            ys = rng.poisson(mu)
            _, ms = C.fit_pois(ys.ravel().astype(float), X)
            ms = ms.reshape(S, C.NDAY)
            n7[i] = C.phi_blocks(ys, ms, 7, p)
            n30[i] = C.phi_blocks(ys, ms, 30, p)
            si = np.repeat(np.repeat(np.arange(S), C.NDAY), ys.ravel())
            di = np.repeat(np.tile(np.arange(C.NDAY), S), ys.ravel())
            tt = (di * 24 + 6 * rng.integers(0, 4, size=len(di))).astype(float)
            g48[i] = C.gap_stats(si + s0, tt, lam)["G48"]
        crit = dict(W7=np.percentile(n7, 95), W30=np.percentile(n30, 95), G48=np.percentile(g48, 95))
        res = {}
        for e in ES:
            hit = dict(W7=0, W30=0, G48=0, KNOX=0)
            E7 = []
            for _ in range(a.nrep):
                evs = draw(mu, ev, s0, e, rng)
                yy = C.daily([dict(t=x["t"], season=x["season"]) for x in evs], s0, s1)
                _, m_ = C.fit_pois(yy.ravel().astype(float), X)
                m_ = m_.reshape(S, C.NDAY)
                f7, f30 = C.phi_blocks(yy, m_, 7, p), C.phi_blocks(yy, m_, 30, p)
                hit["W7"] += f7 > crit["W7"]
                hit["W30"] += f30 > crit["W30"]
                E7.append(f7 - n7.mean())
                se = np.array([x["season"] for x in evs])
                tt = np.array([x["t"] for x in evs])
                g = C.gap_stats(se, tt, lam)
                hit["G48"] += g["G48"] > crit["G48"]
                K = C.m3(evs, rng, a.nperm)
                hit["KNOX"] += K["p"][1, 1] <= 0.05
            res[str(e)] = {k: v / a.nrep for k, v in hit.items()}
            res[str(e)]["mean_E_W7"] = float(np.mean(E7))
            print(basin, e, res[str(e)], flush=True)
        out[basin] = res
    os.makedirs(a.out, exist_ok=True)
    json.dump(dict(tier=a.tier, nrep=a.nrep, es=ES, power=out), open(os.path.join(a.out, f"{a.tier}_power.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
