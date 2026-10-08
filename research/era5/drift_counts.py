"""Do pipeline A's thresholded counts drift before 2001 at fixed depth?

Runs the test fixed in advance in drift_counts_plan.md (committed fe0d0cc,
before any statistic was computed), on the full pipeline A track population
research/era5/hf_history/results/all_tracks.csv.gz (8b109b5). Pipeline A is an
ERA5 PROXY: an event is a track whose 800 km ocean gust index reaches 71.7 kt.
Nothing before 2001 is validated.

    python3 research/era5/drift_counts.py   # writes drift_counts-result.txt beside it
"""
import csv
import gzip
import io
import math
import os

import numpy as np

from drift_check import ols, t_crit, t_cdf, power

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "hf_history", "results", "all_tracks.csv.gz")
THR = 71.7
PRE = list(range(1979, 2001))
REF = list(range(2004, 2026))


def load():
    with gzip.open(SRC, "rt") as fh:
        rows = list(csv.DictReader(fh))
    season = np.array([int(r["season"]) for r in rows])
    basin = np.array([r["basin"] for r in rows])
    gust = np.array([float(r["gust800_kt"]) for r in rows])
    minp = np.array([float(r["minp"]) for r in rows])
    return season, basin, gust, minp


def bins_of(minp, width):
    lo, hi = 945.0, 1005.0
    b = np.floor((np.clip(minp, lo - 1e-6, hi + 1e-6) - lo) / width).astype(int)
    b[minp < lo] = -1
    b[minp >= hi] = 999
    return b


def strata(basin, minp, width):
    b = bins_of(minp, width)
    return np.array(["%s:%d" % (x, y) for x, y in zip(basin, b)])


def excess(season, stratum, ex, ref_seasons, seasons, mask=None):
    """Observed minus depth-expected events per season; rates from ref_seasons."""
    if mask is None:
        mask = np.ones(len(season), bool)
    ref = np.isin(season, ref_seasons) & mask
    rate = {}
    pooled = {}
    season, stratum, ex = season[mask], stratum[mask], ex[mask]
    ref = ref[mask]
    mask = np.ones(len(season), bool)
    for k in np.unique(stratum):
        m = ref & (stratum == k)
        rate[k] = ex[m].mean() if m.sum() > 0 else np.nan
        pooled[k] = ex[mask & (stratum == k)].mean()
    r = np.array([rate[k] if not np.isnan(rate.get(k, np.nan)) else pooled[k] for k in stratum])
    nofill = sum(1 for k in np.unique(stratum[mask & np.isin(season, seasons)]) if np.isnan(rate.get(k, np.nan)))
    O = np.array([ex[mask & (season == s)].sum() for s in seasons], float)
    E = np.array([r[mask & (season == s)].sum() for s in seasons], float)
    return O, E, nofill


def logit_fe(y, strat, x, groups):
    """Logistic regression y ~ stratum dummies + x, by IRLS; season-clustered
    sandwich se for the x coefficient. Strata with no variation are dropped
    (they carry no information on x)."""
    keep = np.ones(len(y), bool)
    for k in np.unique(strat):
        m = strat == k
        if y[m].min() == y[m].max():
            keep &= ~m
    y, strat, x, groups = y[keep], strat[keep], x[keep], groups[keep]
    ks = np.unique(strat)
    X = np.zeros((len(y), len(ks) + 1))
    idx = {k: i for i, k in enumerate(ks)}
    for i, k in enumerate(strat):
        X[i, idx[k]] = 1.0
    X[:, -1] = x
    beta = np.zeros(X.shape[1])
    for _ in range(50):
        eta = X @ beta
        p = 1 / (1 + np.exp(-eta))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None])
        g = X.T @ (y - p)
        step = np.linalg.solve(H, g)
        beta += step
        if np.abs(step).max() < 1e-10:
            break
    p = 1 / (1 + np.exp(-(X @ beta)))
    W = p * (1 - p)
    Hinv = np.linalg.inv(X.T @ (X * W[:, None]))
    meat = np.zeros_like(Hinv)
    G = np.unique(groups)
    for gg in G:
        m = groups == gg
        sc = X[m].T @ (y[m] - p[m])
        meat += np.outer(sc, sc)
    cov = len(G) / (len(G) - 1) * Hinv @ meat @ Hinv
    b, se = beta[-1], math.sqrt(cov[-1, -1])
    ame = float(np.mean(W)) * b       # average marginal effect of x on P(exceed)
    return b, se, ame, int(keep.sum()), len(ks)


def run():
    buf = io.StringIO()

    def P(*a):
        print(*a, file=buf)
        print(*a, flush=True)

    def L(label, r, unit="events/decade"):
        lo, hi = r["b"] - r["crit"] * r["se"], r["b"] + r["crit"] * r["se"]
        P("   %-58s n=%2d  %+7.2f %s  95%% CI [%+6.2f, %+6.2f]  t=%+5.2f  p=%.3f" % (
            label, r["n"], r["b"], unit, lo, hi, r["t"], r["p"]))
        return lo, hi

    season, basin, gust, minp = load()
    ex = (gust >= THR).astype(float)
    P("Pipeline A thresholded counts at fixed depth, 1979-2000 (drift_counts.py)")
    P("=" * 78)
    P("Plan committed before running: drift_counts_plan.md (fe0d0cc). Pipeline A is an ERA5")
    P("PROXY; nothing before 2001 is validated. OLS slopes per decade on season, n = seasons.")
    P("")
    P("0. DATA")
    P("   %d tracks, seasons %d-%d; events (index >= %.1f kt): %d (atl %d, pac %d)." % (
        len(season), season.min(), season.max(), THR, int(ex.sum()),
        int(ex[basin == "atl"].sum()), int(ex[basin == "pac"].sum())))
    for lab, ss in (("1979-2000", PRE), ("2004-2025", REF)):
        m = np.isin(season, ss)
        P("   %s: %d tracks, %d events, %.1f events per season" % (lab, m.sum(), int(ex[m].sum()), ex[m].sum() / len(ss)))
    P("")

    st5 = strata(basin, minp, 5)
    P("1. PRIMARY TEST T1: depth-standardised excess, basin x 5 hPa strata, rates from 2004-2025")
    O, E, nf = excess(season, st5, ex, REF, PRE)
    P("   season  observed  expected  excess")
    for s, o, e in zip(PRE, O, E):
        P("   %d   %7.0f  %8.1f  %+6.1f" % (s, o, e, o - e))
    P("   strata in 1979-2000 with no 2004-2025 tracks (pooled rate used): %d" % nf)
    P("   mean excess 1979-2000: %+.2f events per season (observed %.1f, expected %.1f)" % (
        (O - E).mean(), O.mean(), E.mean()))
    t1 = ols(PRE, O - E)
    t1_lo, t1_hi = L("T1: excess (obs - exp | depth), 1979-2000", t1)
    L("    observed events alone, 1979-2000", ols(PRE, O))
    L("    expected from depth alone, 1979-2000", ols(PRE, E))
    O2, E2, _ = excess(season, st5, ex, PRE, PRE)
    L("T1 with rates from 1979-2000 itself", ols(PRE, O2 - E2))
    O3, E3, _ = excess(season, st5, ex, REF, REF)
    L("T1 control: excess 2004-2025 (rates 2004-2025)", ols(REF, O3 - E3))
    P("")

    P("2. SECONDARY TEST T2: logistic, exceed ~ stratum FE + season, 1979-2000, season-clustered se")
    m = np.isin(season, PRE)
    b, se, ame, nk, ns = logit_fe(ex[m], st5[m], (season[m] - 1989.5) / 10.0, season[m])
    tpd = m.sum() / len(PRE)
    crit = t_crit(len(PRE) - 2)
    P("   season coefficient %+.4f per decade (log-odds), se %.4f, z=%+.2f; %d tracks in %d strata with variation" % (
        b, se, b / se, nk, ns))
    P("   implied: %+.2f events/decade (95%% CI %+.2f to %+.2f), from AME x %.0f tracks per season" % (
        ame * tpd, (b - crit * se) / b * ame * tpd if b else float("nan"), (b + crit * se) / b * ame * tpd if b else float("nan"), tpd))
    m2 = np.isin(season, REF)
    b2, se2, ame2, _, _ = logit_fe(ex[m2], st5[m2], (season[m2] - 2014.5) / 10.0, season[m2])
    P("   control 2004-2025: coefficient %+.4f, se %.4f, z=%+.2f; implied %+.2f events/decade" % (
        b2, se2, b2 / se2, ame2 * m2.sum() / len(REF)))
    P("   Descriptive, per 5 hPa bin, both basins pooled: trend in share reaching 71.7 kt, 1979-2000")
    P("   (no multiplicity correction)")
    b5 = bins_of(minp, 5)
    for k in sorted(np.unique(b5)):
        lab = "<945" if k == -1 else (">=1005" if k == 999 else "%d-%d" % (945 + 5 * k, 950 + 5 * k))
        sh = []
        nn = 0
        for s in PRE:
            mm = (season == s) & (b5 == k)
            nn += mm.sum()
            sh.append(ex[mm].mean() if mm.sum() else np.nan)
        sh = np.array(sh)
        ok = ~np.isnan(sh)
        if ok.sum() >= 10 and nn >= 100 and 0 < np.nanmean(sh) < 1:
            r = ols(np.array(PRE)[ok], sh[ok])
            P("     minp %-7s tracks %5d  mean share %.3f  trend %+.4f /decade  t=%+5.2f" % (lab, nn, np.nanmean(sh), r["b"], r["t"]))
    P("")

    P("3. POWER (as planned): shift every 1979-2000 index by d * (season - 1989.5)/10 kt")
    base = t1["b"]
    ds = np.linspace(0, 1.5, 151)
    induced = []
    for d in ds:
        g2 = gust.copy()
        mm = np.isin(season, PRE)
        g2[mm] = gust[mm] + d * (season[mm] - 1989.5) / 10.0
        e2 = (g2 >= THR).astype(float)
        Ox, Ex, _ = excess(season, st5, e2, REF, PRE)
        induced.append(ols(PRE, Ox - Ex)["b"] - base)
    induced = np.array(induced)
    for eff in (1.0, 3.0, 6.0):
        d = float(np.interp(eff, induced, ds))
        P("   count drift %+.0f events/decade <- index drift %.2f kt/decade; power %.0f%%" % (
            eff, d, 100 * power(eff, t1["se"], len(PRE) - 2)))
    d65 = float(np.interp(0.65, ds, induced))
    P("   c9dc994's +0.65 kt/decade would add %+.2f events/decade to T1; power %.0f%%" % (
        d65, 100 * power(d65, t1["se"], len(PRE) - 2)))
    P("   minimum count drift detectable with 80%% power: %.2f events/decade" % ((t1["crit"] + 0.86) * t1["se"]))
    P("")

    P("4. ROBUSTNESS (pre-specified)")
    for w in (2, 10):
        stw = strata(basin, minp, w)
        Ow, Ew, _ = excess(season, stw, ex, REF, PRE)
        L("T1 with %d hPa bins" % w, ols(PRE, Ow - Ew))
    for bb in ("atl", "pac"):
        Ob, Eb, _ = excess(season, st5, ex, REF, PRE, mask=(basin == bb))
        L("T1, %s only" % bb, ols(PRE, Ob - Eb))
    L("tracks deeper than 960 hPa, 1979-2000", ols(PRE, [((season == s) & (minp <= 960)).sum() for s in PRE]), "tracks/decade")
    L("mean track minp, 1979-2000", ols(PRE, [minp[season == s].mean() for s in PRE]), "hPa/decade")
    L("tracks per season, 1979-2000", ols(PRE, [(season == s).sum() for s in PRE]), "tracks/decade")
    P("")

    P("6. POST HOC (not in the plan; added after seeing section 1): the LEVEL of the excess")
    P("   T1 tests a ramp within 1979-2000. Section 1 also shows the whole period sits below the")
    P("   depth expectation. Excess per season, rates from 2004-2025, compared with the 2004-2025")
    P("   excess on the same rates (zero mean by construction), Welch t on seasons:")
    d_ref = O3 - E3
    for lab, ss in (("1979-2000", PRE), ("1979-1996", list(range(1979, 1997))), ("1997-2000", list(range(1997, 2001)))):
        Dx = np.array([o - e for s0, o, e in zip(PRE, O, E) if s0 in ss])
        se_w = math.sqrt(Dx.var(ddof=1) / len(Dx) + d_ref.var(ddof=1) / len(d_ref))
        tw = (Dx.mean() - d_ref.mean()) / se_w
        P("   %-10s mean excess %+6.2f events/season (%+.1f%% of expected), se %.2f, t=%+.2f, n=%d" % (
            lab, Dx.mean(), 100 * Dx.mean() / np.mean([e for s0, e in zip(PRE, E) if s0 in ss]), se_w, tw, len(Dx)))
    for bb in ("atl", "pac"):
        Ob, Eb, _ = excess(season, st5, ex, REF, PRE, mask=(basin == bb))
        Or, Er, _ = excess(season, st5, ex, REF, REF, mask=(basin == bb))
        Db, Dr = Ob - Eb, Or - Er
        se_w = math.sqrt(Db.var(ddof=1) / len(Db) + Dr.var(ddof=1) / len(Dr))
        P("   %s 1979-2000 mean excess %+6.2f events/season (%+.1f%%), t=%+.2f" % (bb, Db.mean(), 100 * Db.mean() / Eb.mean(), Db.mean() / se_w))
    P("   What the step does to a whole-record trend (seasons 1979-2025, both basins; 2001-03 kept):")
    ALL = list(range(1979, 2026))
    Oa, Ea, _ = excess(season, st5, ex, REF, ALL)
    L("observed pipeline A events, 1979-2025", ols(ALL, Oa))
    L("expected from depth alone (2004-2025 rates), 1979-2025", ols(ALL, Ea))
    L("excess, 1979-2025", ols(ALL, Oa - Ea))
    P("")

    P("5. DECISION RULE (from the plan), on the T1 95%% CI [%+.2f, %+.2f] events/decade:" % (t1_lo, t1_hi))
    if t1_hi < 2 and (b > 0) == (t1["b"] > 0):
        verdict = "(a) upper bound below +2: counts era-comparable to within the CI"
    elif t1_lo > 0 and b > 0:
        verdict = "(b) lower bound above 0: count drift detected"
    else:
        verdict = "(c) not settled; the CI replaces the +1 to +6 band"
    P("   " + verdict)
    P("   T2 sign: %s; T1 sign: %s" % ("+" if b > 0 else "-", "+" if t1["b"] > 0 else "-"))
    P("")
    P("7. READING")
    for t in READING.strip("\n").split("\n"):
        P(t)
    out = os.path.join(HERE, "drift_counts-result.txt")
    with open(out, "w") as fh:
        fh.write(buf.getvalue())
    print("wrote", out)


READING = """
   Written for the numbers above as committed. A fresh agent recomputed all eleven
   quoted claims from all_tracks.csv.gz with its own code; all matched.

   Pre-registered answer (rule (c)): a RAMP in depth-adjusted counts within 1979-2000 is
   not detectable. T1 = +2.77 events/decade, 95% CI -3.0 to +8.6; T2 agrees in sign
   (+3.6, z = +1.57). Power is 53% at +6 events/decade and 17% at +3; the smallest drift
   detectable with 80% power is about +8. The full track population did not buy power,
   because the noise is season-to-season variation in the excess, not track sampling.
   By basin the ramp is Pacific (+4.5, t = +2.53) and not Atlantic (-1.7); with no
   multiplicity correction and an opposite-signed partner this is not a finding.

   Post hoc (section 6, labelled as such): the LEVEL differs. At the same ERA5 depth,
   1979-2000 has 7.2 fewer gust-based events per season than 2004-2025 rates predict
   (-7.8%, Welch t = -2.77; 1979-1996 -9.2%, t = -3.08), in both basins (-8.0%, -7.5%).
   That offset carries pipeline A's whole-record trend: observed events +2.39 /decade
   (t = +2.35) over 1979-2025, depth-expected events -0.19 (t = -0.28), excess +2.58
   (t = +2.78). Its direction is the one c9dc994 predicted from the per-fix gust drift,
   but it was not the planned test, it has not been replicated, and a real change in gust
   at fixed depth would look the same. The independent buoy and ship comparison is the
   test that can separate those.

   For decision 1: thresholded counts before 2001 can be used for variation within that
   era, but not for levels or trends across 2001. Pipeline A's 1979-2025 rise is not
   supported as a climate signal: the depth-only expectation is flat. The planned test
   cannot exclude a ramp of up to about +8.6 events/decade, and the post hoc level step
   is the size of the whole trend.
"""


if __name__ == "__main__":
    run()
