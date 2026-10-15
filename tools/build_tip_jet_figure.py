"""Figure for docs/tip-jets.html: average sea-level pressure and 10 m wind around Cape Farewell for the three tip jet setups,
storm-force and hurricane-force episodes (ERA5 PROXY). Reads the fields pulled by research/era5/tip_jets/stage1_pull.py
(shared folder tar, extracted to FIELDS_DIR), episodes.csv and pull_times.csv from that study. Setups are the widened classes
described in research/era5/tip_jets/README.md.
usage: FIELDS_DIR=<dir with YYYYMMDDHH.npz> EPISODES=<episodes.csv> PULL_TIMES=<pull_times.csv> python3 tools/build_tip_jet_figure.py docs/tip-jets-composites.png
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

FIELDS, EPI, PT = os.environ["FIELDS_DIR"], os.environ["EPISODES"], os.environ["PULL_TIMES"]
OUT = sys.argv[1]
LAT = 90 - 0.25 * (60 + np.arange(120)); LON = -75 + 0.25 * np.arange(280); KT = 1 / 0.514444
names = {f.name for f in font_manager.fontManager.ttflist}
plt.rcParams.update({"font.family": "Inter" if "Inter" in names else "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False,
                     "axes.spines.left": False, "axes.spines.bottom": False, "text.color": "#2b2b2a"})
inb = lambda x, a, b: ((x - a) % 360) < ((b - a) % 360)
pk = pd.read_csv(EPI, parse_dates=["time"])
pk["cls"] = np.where(inb(pk.bearing, 45, 225) & inb(pk.wdir, 225, 315), "S", pk.cls.replace({"S": "other"}))
pt = pd.read_csv(PT, parse_dates=["time"])
load = lambda t: dict(np.load(os.path.join(FIELDS, t.strftime("%Y%m%d%H") + ".npz")))
cm = np.mean([load(t)["msl"] for t in pt[pt.role == "control"].time], 0) / 100
fig, ax = plt.subplots(2, 3, figsize=(13, 7.4), sharex=True, sharey=True)
titles = {"S": "Westerly wind, east/south of cape", "E": "Northerly wind, east of cape", "W": "Easterly wind, west of cape"}
levels = np.arange(-30, 3.1, 3)
for j, c in enumerate(["S", "E", "W"]):
    for i, hf in enumerate([False, True]):
        d = pk[(pk.cls == c) & (pk.hf_ep == hf)]; F = [load(t) for t in d.time]
        m = np.mean([f["msl"] for f in F], 0) / 100; u = np.mean([f["u10"] for f in F], 0); v = np.mean([f["v10"] for f in F], 0)
        a = ax[i, j]; a.set_facecolor("white")
        cs = a.contourf(LON, LAT, np.minimum(m - cm, 0), levels=levels, cmap="Blues_r", extend="min", alpha=0.9)
        a.contour(LON, LAT, m, levels=np.arange(940, 1040, 8), colors="#6b6a66", linewidths=0.6)
        s = 9; a.quiver(LON[::s], LAT[::s], u[::s, ::s], v[::s, ::s], scale=230, width=0.0022, color="#0b3a6b")
        a.plot(-43.9, 59.77, "o", ms=6, mfc="#eb6834", mec="white", zorder=5); a.text(-44.8, 60.7, "Cape Farewell", fontsize=8.5, color="#b44a1c", ha="right", zorder=5)
        k = np.unravel_index(m.argmin(), m.shape)
        a.text(LON[k[1]], LAT[k[0]], f"L {m.min():.0f}", fontsize=10, fontweight="bold", ha="center", va="center", color="#0b0b0b", zorder=6)
        a.set_xlim(-70, -10); a.set_ylim(48, 72); a.grid(True, color="#e8e7e2", lw=0.6); a.tick_params(length=0, labelsize=8.5, colors="#6d6b66")
        a.set_xticks([-60, -40, -20]); a.set_xticklabels(["60W", "40W", "20W"]); a.set_yticks([50, 60, 70]); a.set_yticklabels(["50N", "60N", "70N"])
        a.set_title(f"{'Hurricane force' if hf else 'Storm force'} ({len(d)} episodes)", fontsize=10.5, loc="left", fontweight="bold", pad=6)
        if i == 0: a.annotate(titles[c], (0, 1.16), xycoords="axes fraction", fontsize=11.5, fontweight="bold", color="#0a3f73")
fig.colorbar(cs, ax=ax, shrink=0.55, pad=0.02, label="hPa lower than normal", ticks=[-30, -20, -10, 0])
fig.savefig(OUT, dpi=120, bbox_inches="tight", facecolor="white")
