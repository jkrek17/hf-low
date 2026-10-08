"""Power for the pre-registered primary test (PREREGISTRATION.md, "Power"). Run before the trough
coefficient is read. Outcomes are redrawn from the covariate-only fitted model plus a planted trough effect.

usage: power.py WORK_DIR   -> results/power.csv, results/power.txt
"""
import sys, os
import numpy as np, pandas as pd
import tblib as T

work = sys.argv[1]
os.makedirs("results", exist_ok=True)
N, R = T.load_tables(work)
std = T.Std(R)
rng = np.random.default_rng(20261008)
PLANT = [1.0, 1.05, 1.10, 1.15, 1.20, 1.25]
NDRAW = 400
rows = []
for d, col in T.DEFS.items():
    S = N[N.elig & N[col].notna()].reset_index(drop=True)
    X = T.design(S, std, col)
    y = S.bomb.values
    X0 = X[:, :-1]
    b0, _ = T.fit_w(X0, y)                 # covariate-only fit; uses the outcome but not the trough
    z = X[:, -1]
    eta0 = X0 @ b0
    for orr in PLANT:
        det = []
        for _ in range(NDRAW):
            p = 1 / (1 + np.exp(-(eta0 + np.log(orr) * z)))
            ys = (rng.random(len(p)) < p).astype(float)
            b, _ = T.fit_w(X, ys)
            se, pv, _ = T.cluster_inference(X, ys, b, S.season.values)
            det.append(b[-1] > 0 and pv < 0.05)
        rows.append(dict(definition=d, planted_OR=orr, power=float(np.mean(det)), draws=NDRAW,
                         fixes=len(S), events=int(y.sum()), seasons=S.season.nunique()))
        print(rows[-1], flush=True)
P = pd.DataFrame(rows)
P.to_csv("results/power.csv", index=False)
with open("results/power.txt", "w") as f:
    f.write("Power of the primary test (CR1 season-clustered Wald, positive OR, p < 0.05), planted OR per SD, "
            f"{NDRAW} draws. Draws are independent given covariates, so power is an upper bound.\n")
    f.write(P.to_string(index=False) + "\n")
    for d in T.DEFS:
        q = P[P.definition == d]
        ok = q[q.power >= 0.8]
        mde = ok.planted_OR.min() if len(ok) else float("nan")
        if len(ok) and ok.planted_OR.min() > q.planted_OR.min():
            i = q.index[q.planted_OR == ok.planted_OR.min()][0]
            a, b_ = q.loc[i - 1], q.loc[i]
            mde = a.planted_OR + (0.8 - a.power) / (b_.power - a.power) * (b_.planted_OR - a.planted_OR)
        f.write(f"{d}: minimum detectable OR per SD at 80% power = {mde:.3f}\n")
print(open("results/power.txt").read())
