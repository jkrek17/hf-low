"""Extras requested by the claim review: monthly counts, factor SDs in physical units, COLD adjusted for calendar month
(post hoc, widened S). ERA5 PROXY."""
import os, sys
import numpy as np, pandas as pd
sys.argv = ["x", "s2"]
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
exec(open(os.path.join(HERE, "posthoc.py")).read().split("rows = []")[0])
out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
pk["month"] = pk.time.dt.month
three = pk[pk.cls.isin(["S", "E", "W"])]
P("Episodes per month in the three setups (S widened):", three.groupby("month").size().to_dict(), "; by class:", three.groupby(["cls", "month"]).size().unstack().to_dict("index"))
P("HF share by month, three setups:", three.groupby("month").hf_ep.mean().round(2).to_dict())
P("\nWithin-class SD of each factor (physical units: hPa, km, hPa/100 km, hPa, K):")
for c in ["S", "E", "W"]: P(c, pk[pk.cls == c][FACT].std().round(2).to_dict())
P("\nCOLD by month in S, mean K:", pk[pk.cls == "S"].groupby("month").COLD.mean().round(1).to_dict())
P("\nCOLD and GH ORs per SD in S with month dummies and LOW_DEPTH added (all seasons, ridge 0.5, 2,000 season-block bootstraps):")
d = pk[pk.cls == "S"].copy()
for f in FACT: d[f + "_z"] = (d[f] - d[f].mean()) / d[f].std()
for m in (11, 12, 1, 2): d[f"m{m}"] = (d.month == m).astype(float)
dm = [f"m{m}" for m in (11, 12, 1, 2)]
for f in ("COLD", "GH", "LOW_DEPTH", "GRAD"):
    b = boot(d, f, dm); e = one(d, f, dm)
    P(f, f"OR {np.exp(e):.2f} [{np.exp(np.percentile(b, 2.5)):.2f}, {np.exp(np.percentile(b, 97.5)):.2f}]")
P("\nControls pulled:", 'see pull_times.csv')
pt = pd.read_csv(os.path.join(RES, "pull_times.csv")); P(pt.role.value_counts().to_dict())
open(os.path.join(RES, "extras.txt"), "w").write("\n".join(out)); print("\n".join(out))
