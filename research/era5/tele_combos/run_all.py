"""Run every analysis in PREREGISTRATION.md. ERA5 PROXY, pipeline A.

usage: run_all.py OUT_DIR
Writes OUT_DIR/results.csv (one row per estimated interaction), OUT_DIR/summary.txt.
"""
import sys
import os
import numpy as np
import pandas as pd
import core

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
TESTS = {
    "T1": ("pac", "ONI", "MJOWP"), "T2": ("pac", "PNA", "MJODL"),
    "T3": ("atl", "NAO", "PNA"), "T4": ("atl", "NAO", "ONI"), "T5": ("atl", "NAO", "SPV"),
}
T = core.load_tracks()
I = core.predictors()
ROWS = []
CORNERS = []


def add(tier, test, spec, outcome, est, se, lo, hi, p, n, ns, extra=None):
    basin, A, B = TESTS[test]
    r = dict(tier=tier, test=test, basin=basin, A=A, B=B, spec=spec, outcome=outcome, gamma=est, se_cr1=se,
             ci_lo=lo, ci_hi=hi, p=p, n_events=n, n_seasons=ns)
    if extra:
        r.update(extra)
    ROWS.append(r)


def run_L(tier, test, spec, sample, s0, s1, now=False, era=False, window="octapr", corners=False):
    basin, A, B = TESTS[test]
    if (B + ("_now" if now else "_lag")) not in I.columns:
        return
    D, E = core.build(T, I, basin, A, B, s0, s1, sample, now=now, era=era, window=window)
    r = core.location_test(E, era=era)
    for k, nm in enumerate(("lon", "lat")):
        add(tier, test, spec, "L:" + nm + ":" + sample, r["gamma"][k], r["se_gamma"][k], r["ci_gamma"][k, 0], r["ci_gamma"][k, 1],
            r["p_lon"] if k == 0 else r["p_lat"], r["n_events"], r["n_seasons"],
            dict(p_joint=r["p_joint"], betaA=r["A"][k], betaB=r["B"][k], y_sd=r["y_sd"][k]))
    if corners:
        for a in (1, -1):
            for b in (1, -1):
                CORNERS.append(dict(test=test, spec=spec, kind="L", a=a, b=b,
                                    lon_total=r["A"][0] * a + r["B"][0] * b + r["gamma"][0] * a * b,
                                    lon_additive=r["A"][0] * a + r["B"][0] * b,
                                    lat_total=r["A"][1] * a + r["B"][1] * b + r["gamma"][1] * a * b,
                                    lat_additive=r["A"][1] * a + r["B"][1] * b))


def run_C(tier, test, spec, sample, s0, s1, now=False, era=False, window="octapr", corners=False, box=None):
    basin, A, B = TESTS[test]
    if (B + ("_now" if now else "_lag")) not in I.columns:
        return
    D, E = core.build(T, I, basin, A, B, s0, s1, sample, now=now, era=era, window=window)
    if box:
        la0, la1, lo0, lo1 = core.BOXES[basin][box]
        E = E[(E.lat >= la0) & (E.lat < la1) & (E.lon >= lo0) & (E.lon < lo1)]
    y = core.daily_counts(D, E)
    r = core.count_test(D, y, era=era)
    add(tier, test, spec, "C:" + sample + (":" + box if box else ""), r["gamma"], r["se_gamma"], r["ci_gamma"][0], r["ci_gamma"][1],
        r["p"], r["n_events"], r["n_seasons"], dict(betaA=r["A"], betaB=r["B"]))
    if corners:
        for a in (1, -1):
            for b in (1, -1):
                CORNERS.append(dict(test=test, spec=spec, kind="C", a=a, b=b,
                                    rr_total=float(np.exp(r["A"] * a + r["B"] * b + r["gamma"] * a * b)),
                                    rr_additive=float(np.exp(r["A"] * a + r["B"] * b))))


def run_S5(test, s0=2004, s1=2025, nboot=1000):
    """HF-equivalent minus all-cyclone interaction on position, season-pairs bootstrap (same seasons)."""
    basin, A, B = TESTS[test]
    if (B + "_lag") not in I.columns:
        return
    D, Eh = core.build(T, I, basin, A, B, s0, s1, "hf")
    _, Ea = core.build(T, I, basin, A, B, s0, s1, "all")
    rng = np.random.default_rng(core.SEED)

    def fit(E):
        cov = core.covariates(E, False)
        X = np.column_stack([E.A, E.B, E.AB, cov.values])
        return X, E[["lon", "lat"]].values.astype(float), E.season.values

    Xh, Yh, gh = fit(Eh.reset_index(drop=True))
    Xa, Ya, ga = fit(Ea.reset_index(drop=True))
    seasons = np.unique(gh)
    ih = {s: np.where(gh == s)[0] for s in seasons}
    ia = {s: np.where(ga == s)[0] for s in seasons}
    bh = np.linalg.lstsq(Xh, Yh, rcond=None)[0][2]
    ba = np.linalg.lstsq(Xa, Ya, rcond=None)[0][2]
    diffs = []
    for _ in range(nboot):
        pick = rng.choice(seasons, size=len(seasons), replace=True)
        jh = np.concatenate([ih[s] for s in pick]); ja = np.concatenate([ia[s] for s in pick])
        diffs.append(np.linalg.lstsq(Xh[jh], Yh[jh], rcond=None)[0][2] - np.linalg.lstsq(Xa[ja], Ya[ja], rcond=None)[0][2])
    diffs = np.array(diffs)
    lo, hi = np.percentile(diffs, [2.5, 97.5], axis=0)
    for k, nm in enumerate(("lon", "lat")):
        add("S5", test, "HF minus all cyclones, 2004-05..2025-26", "L:" + nm + ":hf-minus-all", bh[k] - ba[k], diffs[:, k].std(),
            lo[k], hi[k], np.nan, len(Eh), len(seasons))


# ---------------------------------------------------------------- primary (fixed order T1..T5)
for t in TESTS:
    run_L("P", t, "lagged, Oct-Apr 2004-05..2025-26", "hf", 2004, 2025, corners=True)
    run_C("P", t, "lagged, Oct-Apr 2004-05..2025-26", "hf", 2004, 2025, corners=True)
# ---------------------------------------------------------------- secondary
for t in TESTS:
    run_L("S1", t, "same-time, 2004-05..2025-26", "hf", 2004, 2025, now=True)
    run_C("S1", t, "same-time, 2004-05..2025-26", "hf", 2004, 2025, now=True)
    for bx in core.BOXES[TESTS[t][0]]:
        run_C("S2", t, "box counts, lagged, 2004-05..2025-26", "hf", 2004, 2025, box=bx)
    for sample, tier in (("hf", "S3"), ("all", "S4")):
        run_L(tier, t, "lagged 1979-80..2025-26 + era term", sample, 1979, 2025, era=True)
        run_C(tier, t, "lagged 1979-80..2025-26 + era term", sample, 1979, 2025, era=True)
        run_L(tier, t, "lagged 2001-02..2025-26", sample, 2001, 2025)
        run_C(tier, t, "lagged 2001-02..2025-26", sample, 2001, 2025)
    run_S5(t)
    run_L("S6", t, "lagged, Jun-May 2004-05..2025-26", "hf", 2004, 2025, window="junmay")
    run_C("S6", t, "lagged, Jun-May 2004-05..2025-26", "hf", 2004, 2025, window="junmay")

R = pd.DataFrame(ROWS)
# BH: primary family (joint p for L, once per test) and whole directory
key = R.copy()
key["kind"] = key.outcome.str[0]
key["is_L_rep"] = (key.outcome.str.startswith("L:lon"))
key["p_use"] = np.where(key.kind == "L", key.p_joint, key.p)
rep = key[(key.kind == "C") | key.is_L_rep].copy()
prim = rep[rep.tier == "P"]
R["q_primary_family"] = np.nan
R.loc[prim.index, "q_primary_family"] = core.bh(prim.p_use.values)
rest = rep[rep.tier != "S5"]
R["q_all_family"] = np.nan
R.loc[rest.index, "q_all_family"] = core.bh(rest.p_use.values)
R.to_csv(os.path.join(OUT, "results.csv"), index=False, float_format="%.5g")
pd.DataFrame(CORNERS).to_csv(os.path.join(OUT, "corners.csv"), index=False, float_format="%.4g")
print("rows", len(R), "primary family size", len(prim), "all-family size", len(rest))
