"""Design power for the pre-registered tests. Reads outcome MARGINALS only: the index series are moved between
seasons (season-block permutation, A and B moved together, day-of-season kept), which destroys any association
of outcome with index. The spread of the interaction coefficient over permutations is the null SE; the
detectable effect at 80% power (two-sided 5%) is 2.8 x SE; a well-powered null needs 1.96 x SE <= SESOI.

usage: power.py OUT_TXT [NPERM]
"""
import sys
import numpy as np
import pandas as pd
import core

NPERM = int(sys.argv[2]) if len(sys.argv) > 2 else 400
TESTS = {
    "T1": ("pac", "ONI", "MJOWP"), "T2": ("pac", "PNA", "MJODL"),
    "T3": ("atl", "NAO", "PNA"), "T4": ("atl", "NAO", "ONI"), "T5": ("atl", "NAO", "SPV"),
}
S0, S1 = 2004, 2025
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def neff_per_season(D, col):
    x = D[col].values
    x = x - np.nanmean(x)
    n = len(x)
    ac = [np.nansum(x[:-k] * x[k:]) / np.nansum(x * x) for k in range(1, 60)]
    tau = 0.0
    for r in ac:
        if r <= 0:
            break
        tau += r
    return 210.0 / (1 + 2 * tau), tau


T = core.load_tracks()
I = core.predictors()
rng = np.random.default_rng(core.SEED)
say("Design power (season-block permutation, %d permutations). Pipeline A, ERA5 proxy. Seasons %d-%d, Oct 1 + 210 d." % (NPERM, S0, S1 + 1))
say("Effective n per season = 210 / (1 + 2 tau), tau = integrated autocorrelation time (days) of the lagged predictor.")
say()
for name, (basin, A, B) in TESTS.items():
    if (B + "_lag") not in I.columns:
        say("%s %s x %s %s: not run (series absent)" % (name, A, B, basin))
        say()
        continue
    D, E = core.build(T, I, basin, A, B, S0, S1, "hf")
    seasons = sorted(D.season.unique())
    ns = len(seasons)
    nA, tA = neff_per_season(D, "A")
    nB, tB = neff_per_season(D, "B")
    # day arrays for permuting
    Aa = np.zeros((ns, 210)); Ba = np.zeros((ns, 210))
    for i, s in enumerate(seasons):
        d = D[D.season == s]
        Aa[i, d.day.values] = d.A.values
        Ba[i, d.day.values] = d.B.values
    sidx = {s: i for i, s in enumerate(seasons)}
    Es = E.season.map(sidx).values
    Ed = E.day.values
    y = core.daily_counts(D, E)
    Dm = D.copy(); Dm["month"] = Dm.date.dt.month
    sdy = E[["lon", "lat"]].std().values
    g_loc, g_cnt = [], []
    cov_e = core.covariates(E.assign(month=E.gen.dt.month), False).values
    cov_d = core.covariates(Dm, False).values
    for r in range(NPERM):
        perm = rng.permutation(ns)
        Ae = Aa[perm[Es], Ed]; Be = Ba[perm[Es], Ed]
        X = np.column_stack([Ae, Be, Ae * Be, cov_e])
        b, *_ = np.linalg.lstsq(X, E[["lon", "lat"]].values, rcond=None)
        g_loc.append(b[2])
        Ad = np.array([Aa[perm[seasons.index(s)], dd] for s, dd in zip(D.season.values, D.day.values)])
        Bd = np.array([Ba[perm[seasons.index(s)], dd] for s, dd in zip(D.season.values, D.day.values)])
        Xc = np.column_stack([Ad, Bd, Ad * Bd, cov_d])
        g_cnt.append(core.newton(y, Xc)[2])
    g_loc = np.array(g_loc); g_cnt = np.array(g_cnt)
    sd_loc = g_loc.std(0); sd_cnt = g_cnt.std()
    sesoi_loc = 0.10 * sdy
    say("%s  %s: %s x %s   events %d, seasons %d, outcome SD (lon, lat) = (%.1f, %.1f) deg" % (name, basin, A, B, len(E), ns, sdy[0], sdy[1]))
    say("   predictor effective n / season: %s %.1f, %s %.1f  (-> %.0f, %.0f over %d seasons); corr(A,B) on analysis days = %.2f"
        % (A, nA, B, nB, nA * ns, nB * ns, ns, np.corrcoef(D.A, D.B)[0, 1]))
    say("   L  null SE of gamma (lon, lat) = (%.2f, %.2f) deg; detectable at 80%% power = (%.2f, %.2f); SESOI = (%.2f, %.2f); well-powered null possible: lon %s, lat %s"
        % (sd_loc[0], sd_loc[1], 2.8 * sd_loc[0], 2.8 * sd_loc[1], sesoi_loc[0], sesoi_loc[1],
           "yes" if 1.96 * sd_loc[0] <= sesoi_loc[0] else "no", "yes" if 1.96 * sd_loc[1] <= sesoi_loc[1] else "no"))
    say("   C  null SE of gamma = %.3f (log RR); detectable RR = %.3f; SESOI RR 1.05; well-powered null possible: %s"
        % (sd_cnt, np.exp(2.8 * sd_cnt), "yes" if 1.96 * sd_cnt <= np.log(1.05) else "no"))
    say()
open(sys.argv[1], "w").write("\n".join(OUT) + "\n")
