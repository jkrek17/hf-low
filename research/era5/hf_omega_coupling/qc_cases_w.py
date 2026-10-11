"""Q3: 12 drawn cases (seeded random fixes, 6 with couple_w = 1 and 6 with 0) of the 500 hPa ascent field with centres marked.
usage: python3 qc_cases_w.py    writes results/qc_cases.png, results/qc_cases.csv (re-pulls 12 chunks, about 0.14 GB)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "hf_env_composites"))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from omega_features import *
from omega_features import centres
from common import decode, WB2, WB2_T0, LEVELS, destination, sample
D = pd.read_csv(os.path.join(HERE, "results", "features_w.csv.gz"), dtype={"time": str})
D = D[pd.to_datetime(D.time, format="%Y%m%d%H").dt.month.isin([10, 11, 12, 1, 2, 3, 4])]
rng = np.random.default_rng(20261011)
pick = pd.concat([D[D.couple_w == k].sample(6, random_state=int(rng.integers(1e6))) for k in (1, 0)])
fig, ax = plt.subplots(3, 4, figsize=(16, 12))
rows = []
for a, r in zip(ax.ravel(), pick.itertuples()):
    t = pd.Timestamp(r.time[:4] + "-" + r.time[4:6] + "-" + r.time[6:8] + " " + r.time[8:10] + ":00")
    k = int((np.datetime64(t) - WB2_T0) / np.timedelta64(6, "h"))
    arr, _ = decode(WB2, "vertical_velocity", f"{k // 8}.0.0.0")
    F = time_fields_w(arr[k % 8, LEVELS.index(500)], arr[k % 8, LEVELS.index(700)])
    la, lo = destination(r.lat, r.lon, BRG, BR)
    A5 = smooth(np.nan_to_num(sample(F["a500"], la, lo)))
    cen = centres(A5, np.ones_like(A5, bool), 0.8, 1500.0, 3000.0)
    im = a.pcolormesh(BX, BY, A5, cmap="RdBu_r", vmin=-1.5, vmax=1.5, shading="auto")
    a.contour(BX, BY, A5, [0.8], colors="k", linewidths=0.6)
    a.plot(0, 0, "k+", ms=14, mew=2)
    for i, (x, y, v) in enumerate(cen):
        a.plot(x, y, "o", mfc="none", mec="g" if i == 0 else "m", ms=12, mew=2)
    if np.isfinite(r.heading):
        h = np.radians(r.heading); a.arrow(0, 0, 900 * np.sin(h), 900 * np.cos(h), width=40, color="k")
    a.add_patch(plt.Circle((0, 0), 3000, fill=False, ls=":", color="gray"))
    a.set_title(f"{r.track} {r.time} {r.basin}\nnasc {int(r.nasc)} couple_w {int(r.couple_w)} a500_1000 {r.a500_1000:.2f}", fontsize=9)
    a.set_aspect("equal"); a.set_xlim(-5000, 5000); a.set_ylim(-5000, 5000)
    rows.append(dict(track=r.track, time=r.time, basin=r.basin, nasc=r.nasc, couple_w=r.couple_w, drawn_centres=len(cen)))
fig.colorbar(im, ax=ax, shrink=0.6, label="ascent -omega at 500 hPa (Pa/s); green = strongest centre, magenta = others; ring = 3000 km")
fig.savefig(os.path.join(HERE, "results", "qc_cases.png"), dpi=70)
pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", "qc_cases.csv"), index=False)
print(pd.DataFrame(rows))
