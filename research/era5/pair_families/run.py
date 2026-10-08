"""RA-12 test per PREREGISTRATION.md (commit 8fe1f08). usage: python3 run.py [NPERM]"""
import sys, json
import numpy as np
from common import load
from engine import *

NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
T, F = load()


def stats(pr, P):
    s = summarise(pr, P)
    s["f"] = s["rear"] / s["all"] if s["all"] else np.nan
    s["f_obs"] = s["rear_obs"] / s["obs"] if s["obs"] else np.nan
    return s


def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p)
    q = np.empty(m); run = 1.0
    for r, k in enumerate(o[::-1]):
        run = min(run, p[k] * m / (m - r)); q[k] = run
    return q


res, prim = {}, []
for b, G in T.groupby("basin"):
    G = G.reset_index(drop=True)
    P = prep(G, F); st = strata(G); rng = np.random.default_rng(20261008)
    ob = stats(pairs(P), P)
    nl = [stats(pairs(P, perm_shift(P, st, rng)), P) for _ in range(NPERM)]
    nf = np.array([x["f"] for x in nl]); nfo = np.array([x["f_obs"] for x in nl])
    nr = np.array([x["rear"] for x in nl]); na = np.array([x["all"] for x in nl])
    nrobs = np.array([x["rear_obs"] for x in nl]); nobs = np.array([x["obs"] for x in nl])
    ok = ~np.isnan(nf)
    p1 = (1 + np.sum(nf[ok] >= ob["f"] - 1e-12)) / (1 + ok.sum())
    p2 = (1 + np.sum(np.abs(nf[ok] - nf[ok].mean()) >= abs(ob["f"] - nf[ok].mean()) - 1e-12)) / (1 + ok.sum())
    ko = ~np.isnan(nfo)
    p1o = (1 + np.sum(nfo[ko] >= ob["f_obs"] - 1e-12)) / (1 + ko.sum())
    p2o = (1 + np.sum(np.abs(nfo[ko] - nfo[ko].mean()) >= abs(ob["f_obs"] - nfo[ko].mean()) - 1e-12)) / (1 + ko.sum())
    # excess rear pairs at the observed number of pairs, with a binomial-plus-null-noise interval
    mu_f, sd_f = nf[ok].mean(), nf[ok].std(ddof=1)
    ex = ob["rear"] - mu_f * ob["all"]
    sd_ex = np.sqrt(ob["all"] * mu_f * (1 - mu_f) + (sd_f * ob["all"]) ** 2)
    # planted-daughter power: rear = Binomial(all, mu_f) + k vs. permutation distribution of f
    sim = np.random.default_rng(7); crit = np.quantile(nf[ok], 0.975)
    pw = {}
    for k in range(0, 80, 2):
        r = sim.binomial(ob["all"], mu_f, 4000) + k
        pw[k] = float(np.mean(r / (ob["all"] + k) >= crit))
    mde = next((k for k, v in pw.items() if v >= 0.8), None)
    res[b] = dict(tracks=len(G), observed=ob, null_f_mean=float(mu_f), null_f_sd=float(sd_f),
                  null_count_mean=dict(all=float(na.mean()), rear=float(nr.mean()), obs=float(nobs.mean()), rear_obs=float(nrobs.mean())),
                  p_one_sided=float(p1), p_two_sided=float(p2), p_obs_one_sided=float(p1o), p_obs_two_sided=float(p2o),
                  null_f_obs_mean=float(np.nanmean(nfo)), excess_rear=float(ex), excess_rear_se=float(sd_ex),
                  excess_rear_upper95=float(ex + 1.645 * sd_ex), mde_k_80=mde, power_curve=pw, nperm=NPERM)
    prim.append(p1)
    print(b, json.dumps({k: v for k, v in res[b].items() if k != "power_curve"}), flush=True)
q = bh(prim)
for (b, _), qq in zip(sorted(res.items()), q):
    res[b]["q_primary"] = float(qq)
print("BH q", dict(zip(sorted(res), q.round(4))))
json.dump(res, open("results/pair_families.json", "w"), indent=1)
