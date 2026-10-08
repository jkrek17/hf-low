"""Power and false-positive rate of the pre-registered threshold rule (PREREGISTRATION.md, "Power").

Outcomes are redrawn on the real fixes; no held-out outcome is used. The fit-season outcomes set the
covariate coefficients of the simulated truth. A threshold is declared when T1, T2, T3 each have
p < 0.01 (an approximation to q < 0.05 over 32 tests), the post-knot slope exceeds the pre-knot slope,
and the knot lies within the 15th-85th fit percentile.

Scenarios (per basin x outcome x variable, 200 simulations):
  hinge   linear trend as fitted plus a hinge at the fit median, odds ratio 2.0 across knot +/- 0.5 SD
  ramp    linear trend as fitted, no hinge (false-positive rate of the threshold rule)
  x1.15   no covariate-adjusted effect except a straight line of odds ratio 1.15 per SD (power of T1)
  null    no effect of x (false-positive rate of T1)
usage: power.py [NSIM]
"""
import sys, itertools, time
import numpy as np, pandas as pd
from multiprocessing import Pool
import jtlib as J

NSIM = int(sys.argv[1]) if len(sys.argv) > 1 else 200
B = 500
SCEN = ("hinge", "ramp", "x1.15", "null")


def job(args):
    basin, outcome, c, scen = args
    D = J.load_all()
    S, y = J.subset(D, outcome, basin)
    cols = J.cov_cols(outcome)
    fm = (S.season <= J.FIT_LAST).values
    Dg = J.Design(S, fm, cols, [c])
    C = Dg.cov(S)
    z = Dg.x(S, c)
    Xl = np.column_stack([C, z])
    b0, _ = J.fit_logit(Xl[fm], y[fm])
    bC, _ = J.fit_logit(C[fm], y[fm])
    k0 = np.median(z[fm])
    if scen == "hinge":
        eta = Xl @ b0 + 2 * np.log(2) * np.maximum(z - k0, 0)
    elif scen == "ramp":
        eta = Xl @ b0
    elif scen == "x1.15":
        eta = C @ bC + np.log(1.15) * z
    else:
        eta = C @ bC
    p = 1 / (1 + np.exp(-eta))
    rng = np.random.default_rng(abs(hash((basin, outcome, c, scen))) % 2**32)
    det = {"T1": 0, "T2": 0, "T3": 0, "threshold": 0}
    for i in range(NSIM):
        ys = (rng.random(len(p)) < p).astype(int)
        r = J.three_tests(S, ys, c, cols, B=B, seed=i)
        ok = [r[t]["p"] < 0.01 for t in ("T1", "T2", "T3")]
        for t, o in zip(("T1", "T2", "T3"), ok):
            det[t] += o
        inside = 15 <= r["knot_pct"] <= 85
        det["threshold"] += all(ok) and r["slope_post"] > r["slope_pre"] and inside
    return dict(basin=basin, outcome=outcome, var=c, scenario=scen, n_sim=NSIM,
                ev_fit=int(y[fm].sum()), ev_test=int(y[~fm].sum()),
                **{k: v / NSIM for k, v in det.items()})


if __name__ == "__main__":
    combos = list(itertools.product(("atl", "pac"), ("BOMB", "HFON"), ("jet250", "trough_up"), SCEN))
    t0 = time.time()
    with Pool(4) as pool:
        rows = []
        for r in pool.imap_unordered(job, combos):
            rows.append(r)
            print(f"{len(rows)}/{len(combos)} {time.time() - t0:.0f}s {r}", flush=True)
    df = pd.DataFrame(rows).sort_values(["basin", "outcome", "var", "scenario"])
    df.to_csv("results/power.csv", index=False)
    with open("results/power.txt", "w") as f:
        f.write(f"Simulation power of the pre-registered rule, {NSIM} simulations per row, {B} bootstrap resamples,\n"
                "p < 0.01 per test as the stand-in for FDR q < 0.05. Held-out outcomes are not used.\n\n")
        f.write(df.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
