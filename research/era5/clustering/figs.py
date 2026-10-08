"""Figures for the clustering README. usage: figs.py RESULTS_DIR OUT_DIR"""
import json
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R, OUT = sys.argv[1], sys.argv[2]
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
coast = json.loads(re.search(r"=\s*(\{.*\})", open(os.path.join(REPO, "docs/data/coastlines.js")).read(), re.S).group(1))["polygons"]
LAB = ["0-24", "24-48", "48-72", "72-96", "4-7 d", "7-14 d", "14-30 d"]
LAB = ["<=24 h", "24-48 h", "48-72 h", "72-96 h", "4-7 d", "7-14 d", "14-30 d"]

# ---- Figure 1: Knox excess by time separation (post hoc breakdown of the pre-registered statistic)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
col = {"arch": "#1f4e79", "r1": "#c26a00"}
nm = {"arch": "archive (OPC warnings)", "r1": "ERA5 proxy, pipeline A"}
for a, b in zip(ax, ("atl", "pac")):
    for k, (t, off) in enumerate((("arch", -0.18), ("r1", 0.18))):
        p = json.load(open(os.path.join(R, f"{t}_posthoc.json")))["basins"][b]["knox_base"]
        a.bar(np.arange(7) + off, p["ratio_by_bin"], 0.34, color=col[t], label=nm[t])
    a.axhline(1, color="k", lw=0.8)
    a.set_xticks(range(7), LAB, rotation=35, ha="right", fontsize=8)
    a.set_title({"atl": "North Atlantic", "pac": "North Pacific"}[b], fontsize=10)
ax[0].set_ylabel("pairs within 1000 km, observed / expected")
ax[0].legend(fontsize=8, frameon=False)
fig.suptitle("HF lows closer together than chance, by time between them (post hoc breakdown)", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "knox_by_lag.png"), dpi=130)

# ---- Figure 2: local excess beyond basin-wide timing, Atlantic and Pacific cells
def draw_coast(a, lon0, lon1):
    for poly in coast:
        xy = np.array(poly)
        x, y = xy[:, 0], xy[:, 1]
        for sh in (0, 360, -360):
            a.plot(x + sh, y, color="#888", lw=0.6)


fig, ax = plt.subplots(2, 2, figsize=(10, 7))
for i, t in enumerate(("arch", "r1")):
    for j, b in enumerate(("atl", "pac")):
        a = ax[i, j]
        d = pd.read_csv(os.path.join(R, f"{t}_{b}_cells.csv"))
        lon = np.where(d.lon > 180, d.lon - 360, d.lon) if b == "atl" else d.lon
        sc = a.scatter(lon, d.lat, c=d.E_shuf, s=520, marker="s", cmap="RdBu_r", vmin=-0.3, vmax=0.3,
                       edgecolor=np.where(d.q_shuf < 0.05, "k", "none"), linewidth=1.8)
        draw_coast(a, 0, 0)
        a.set_xlim((-100, 15) if b == "atl" else (130, 245))
        a.set_ylim(22, 72)
        a.set_title(f"{nm[t]}, {'North Atlantic' if b == 'atl' else 'North Pacific'}", fontsize=9)
        a.tick_params(labelsize=7)
fig.colorbar(sc, ax=ax, shrink=0.7, label="excess weekly dispersion of passages beyond basin-wide timing (E)")
fig.suptitle("Where HF lows arrive in runs (cells with >= 40 passages; black outline = BH q < 0.05)", fontsize=10)
fig.savefig(os.path.join(OUT, "local_cells.png"), dpi=120)
