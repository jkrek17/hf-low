"""POST HOC (not in the plan; labelled as such in the README). Run after run.py.
PH1: variance to mean ratio of the deseasonalised weekly counts (a Poisson count has 1).
PH2: the same-week primary correlation recomputed on non-overlapping blocks of 2, 3, 5, 6, 10 weeks (all divisors of 30 shown, none picked),
     deseasonalised by block-of-season means, one-sided season-block permutation (10,000), to see whether the correlation
     grows when count noise is averaged down."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run as R

rng = np.random.default_rng(R.SEED + 9)
out = {"PH1": {}, "PH2": []}
for s in ("ARCH", "PRX", "DEP"):
    for b in ("atl", "pac"):
        X, _ = R.block(s, b)
        x = R.deseason_week(X)
        out["PH1"][f"{s}_{b}"] = dict(mean=float(X.mean()), resid_var=float(x.var() * 22 / 21), var_over_mean=float(x.var() * 22 / 21 / X.mean()))
    A, _ = R.block(s, "atl"); P, _ = R.block(s, "pac")
    for m in (1, 2, 3, 5, 6, 10):
        a = R.deseason_week(A.reshape(22, NW // m, m).sum(2)) if (NW := R.NW) else None
        b = R.deseason_week(P.reshape(22, NW // m, m).sum(2))
        r = R.corr(a, b)
        null = R.perm_stats(a, b, rng, nperm=10000)
        lo, hi = R.boot_ci(a, b, rng, nboot=3000)
        out["PH2"].append(dict(sample=s, block_weeks=m, r=r, ci_lo=float(lo), ci_hi=float(hi), p_neg=float((1 + (null <= r).sum()) / 10001)))
json.dump(out, open(os.path.join(R.OUT, "posthoc.json"), "w"), indent=1, default=float)
print(json.dumps(out["PH1"], indent=1))
import pandas as pd
print(pd.DataFrame(out["PH2"]).round(4).to_string())
