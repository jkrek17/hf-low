"""Surprise audit for the post hoc results (S widened; GH negative in S, COLD positive in S): even/odd season consistency, a
within-season label shuffle (placebo), and a check that GRAD/LOW_DEPTH are not simply the wind itself. ERA5 PROXY, post hoc."""
import os, sys
import numpy as np, pandas as pd
sys.argv = ["x", "s2"]
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
exec(open(os.path.join(HERE, "posthoc.py")).read().split("rows = []")[0])     # reuse fit/one/boot/bh and the s2 class labels
out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
rng = np.random.default_rng(7)
P("class factor | OR even seasons | OR odd seasons | within-season shuffle p (2,000)")
for c in ["S", "E", "W"]:
    d = pk[pk.cls == c].copy()
    for f in FACT: d[f + "_z"] = (d[f] - d[f].mean()) / d[f].std()
    for f in FACT:
        e = d[d.season % 2 == 0]; o = d[d.season % 2 == 1]
        real = one(d, f); null = []
        for _ in range(2000):
            dd = d.copy(); dd["hf_ep"] = dd.groupby("season").hf_ep.transform(lambda x: rng.permutation(x.values))
            null.append(one(dd, f))
        null = np.array(null)
        P(c, f, f"| {np.exp(one(e, f)):.2f} (n={len(e)}, hf={int(e.hf_ep.sum())}) | {np.exp(one(o, f)):.2f} (n={len(o)}, hf={int(o.hf_ep.sum())}) | p={(np.abs(null) >= abs(real)).mean():.4f}")
# is GRAD just the wind? correlation of each factor with the maximum box gust itself within class (peak time)
P("\nSpearman correlation of each factor with the box gust (kt), within class:")
for c in ["S", "E", "W"]:
    d = pk[pk.cls == c]; P(c, {f: round(d[f].corr(d.box_gust_kt, method="spearman"), 2) for f in FACT})
open(os.path.join(RES, "audit.txt"), "w").write("\n".join(out)); print("\n".join(out))
