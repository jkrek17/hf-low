"""Figure for the README: chance that a 00/12 UTC time reaches a 54 kt ERA5 gust, by ingredient band, featured regime vs the others (PROXY)."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in fm.findSystemFonts():
    if "/inter/" in f.lower(): fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": ["Inter", "DejaVu Sans"], "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
                     "axes.edgecolor": "#888", "axes.grid": True, "axes.axisbelow": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.7,
                     "xtick.color": "#444", "ytick.color": "#444", "font.size": 10})
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
Q = pd.read_csv(os.path.join(RES, "q1_table.csv"))
Q = Q[Q.ingredient != "(all times in regime)"]
BANDS = ["bottom half", "50th-80th", "80th-95th", "top 5%"]
PLOTS = [("Z1", "BARRIER", "DEPTH", "North Gulf coast, southeasterly flow"), ("Z2", "BARRIER", "DEPTH", "Panhandle coast, southeasterly flow"),
         ("Z3", "GAP", "DEPTH", "Cook Inlet / Shelikof, northerly flow"),
         ("Z1", "BARRIER", "COASTPILE", None), ("Z2", "BARRIER", "COASTPILE", None), ("Z3", "GAP", "INTERIOR_DP", None)]
XL = {"DEPTH": "how deep the deepest Gulf low is", "COASTPILE": "pressure higher at the coast than offshore", "INTERIOR_DP": "pressure higher inland than at the strait"}
LAB = {"BARRIER": "southeasterly", "GAP": "northerly"}
OFF = {("Z1", "OFFSHORE"): -2.2, ("Z3", "OTHER"): -2.5}
fig, axs = plt.subplots(2, 3, figsize=(11.5, 6.4), sharey=True)
for ax, (zn, rg, ing, title) in zip(axs.flat, PLOTS):
    d = Q[(Q.zone == zn) & (Q.ingredient == ing)]
    for r in d.regime.unique():
        x = d[d.regime == r].set_index("band").reindex(BANDS)
        if r == rg:
            ax.plot(range(4), x.pct_event, color="#d1495b", lw=2.4, marker="o", zorder=3)
            ax.vlines(range(4), x.pct_event_lo, x.pct_event_hi.clip(upper=45), color="#d1495b", lw=1.2, alpha=.55, zorder=2)
            for i_, hi_ in enumerate(x.pct_event_hi):
                if hi_ > 45: ax.annotate("", xy=(i_, 45), xytext=(i_, 42), arrowprops=dict(arrowstyle="->", color="#d1495b", alpha=.6))
            ax.text(3.08, x.pct_event.iloc[-1], LAB[rg], color="#d1495b", va="center", fontsize=9, fontweight="bold")
        else:
            ax.plot(range(4), x.pct_event, color="#9aa5b1", lw=1.2, marker=".", zorder=1)
            ax.text(3.08, x.pct_event.iloc[-1] + OFF.get((zn, r), 0), r.lower(), color="#7b8794", va="center", fontsize=8)
    ax.set_xticks(range(4)); ax.set_xticklabels(["lowest\nhalf", "50-80%", "80-95%", "top\n5%"], fontsize=8.5)
    ax.set_xlim(-0.3, 4.1); ax.set_ylim(0, 45)
    ax.set_xlabel(XL[ing], fontsize=9, color="#555")
    if title: ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", color="#222")
axs[0, 0].set_ylabel("chance a 12-hourly time\nreaches a 54 kt gust (%)", fontsize=9.5)
axs[1, 0].set_ylabel("chance a 12-hourly time\nreaches a 54 kt gust (%)", fontsize=9.5)
fig.text(0.01, 0.005, "ERA5 reanalysis gust at 0.25 degree, a PROXY (not observation); October-April 2004-05 to 2025-26. Bands are percentiles of all times in the same flow.\n"
         "Whiskers: 95% season-resampled range, arrow = runs off the chart (rough in the top bands, only 9-13 times).", fontsize=7.5, color="#666")
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig(os.path.join(RES, "fig_chance_by_ingredient.png"), dpi=150)
print("ok")
