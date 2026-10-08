"""Check the daily GH series against PR 58's stage-1 GH at shared 12 UTC times (rule: corr >= 0.999)."""
import os, sys
import numpy as np, pandas as pd
H = os.path.dirname(os.path.abspath(__file__))
g = pd.read_csv(os.path.join(H, "gh_daily.csv"), parse_dates=["date"])
s = pd.read_csv(os.path.join(H, "..", "greenland_jets", "results", "stage1_times.csv"), dtype={"time": str})
s = s[s.time.str[8:10] == "12"].copy()
s["date"] = pd.to_datetime(s.time.str[:8], format="%Y%m%d")
m = s.merge(g, on="date", suffixes=("_pr58", ""))
print("shared 12 UTC times", len(m), "n_pts", g.n_pts.unique())
print("corr", np.corrcoef(m.GH_pr58, m.GH)[0, 1], "max abs diff", (m.GH_pr58 - m.GH).abs().max())
print("missing days", g.GH.isna().sum(), "range", g.GH.min(), g.GH.max())
