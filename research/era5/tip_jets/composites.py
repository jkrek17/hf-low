"""Composites for the three setup classes (S, E, W), storm force vs hurricane force, and the OPC link.
ERA5 PROXY. Run after analyse.py. Output: results/composites.png, results/composites.txt, results/opc_link.txt
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.argv = sys.argv[:1]
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "tip_jets", "fields")
RES = os.path.join(HERE, "results")
LAT = 90 - 0.25 * (60 + np.arange(120)); LON = -75 + 0.25 * np.arange(280)
KT = 1 / 0.514444
pk = pd.read_csv(os.path.join(RES, "episodes.csv"), parse_dates=["time"])
pt = pd.read_csv(os.path.join(RES, "pull_times.csv"), parse_dates=["time"])


def load(t):
    z = np.load(os.path.join(WORK, t.strftime("%Y%m%d%H") + ".npz")); return {k: z[k] for k in z.files}


_inb = lambda x, a, b: ((x - a) % 360) < ((b - a) % 360)
if os.environ.get("CLASS_MODE") == "s2":      # post hoc widened S, see PREREGISTRATION.md
    pk["cls"] = np.where(_inb(pk.bearing, 45, 225) & _inb(pk.wdir, 225, 315), "S", pk.cls.replace({"S": "other"}))
SUF = "_s2" if os.environ.get("CLASS_MODE") == "s2" else ""
ctl = [load(t) for t in pt[pt.role == "control"].time]
cm = np.mean([c["msl"] for c in ctl], 0) / 100
cu = np.mean([c["u10"] for c in ctl], 0); cv = np.mean([c["v10"] for c in ctl], 0)
fig, ax = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True)
txt = []
for j, c in enumerate(["S", "E", "W"]):
    for i, hf in enumerate([False, True]):
        d = pk[(pk.cls == c) & (pk.hf_ep == hf)]
        if len(d) == 0: continue
        F = [load(t) for t in d.time]
        m = np.mean([f["msl"] for f in F], 0) / 100; u = np.mean([f["u10"] for f in F], 0); v = np.mean([f["v10"] for f in F], 0)
        a = ax[i, j]
        cs = a.contourf(LON, LAT, m - cm, levels=np.arange(-30, 31, 3), cmap="RdBu_r", extend="both")
        a.contour(LON, LAT, m, levels=np.arange(940, 1040, 4), colors="k", linewidths=0.5)
        s = 8; a.quiver(LON[::s], LAT[::s], u[::s, ::s], v[::s, ::s], scale=250, width=0.002)
        a.plot(-43.9, 59.77, "k*", ms=11)
        a.set_title(f"{c}: {'hurricane-force' if hf else 'storm-force'} proxy (n={len(d)})", fontsize=10)
        a.set_xlim(-70, -10); a.set_ylim(48, 72)
        ws = np.hypot(u, v) * KT
        txt.append(f"{c} {'HF' if hf else 'S'}: n={len(d)}; mean-low min MSLP {m.min():.1f} hPa at ({LAT[np.unravel_index(m.argmin(), m.shape)[0]]:.1f}N, "
                   f"{LON[np.unravel_index(m.argmin(), m.shape)[1]]:.1f}E); composite peak wind {ws.max():.0f} kt")
fig.colorbar(cs, ax=ax, shrink=0.7, label="MSLP minus mean of 469 comparison times that were not events (hPa); black = MSLP, arrows = 10 m wind")
fig.suptitle("Cape Farewell box gust events, Nov-Mar 2004-05 to 2025-26 (ERA5 proxy). Star = Cape Farewell")
fig.savefig(os.path.join(RES, "composites" + SUF + ".png"), dpi=110)
open(os.path.join(RES, "composites" + SUF + ".txt"), "w").write("\n".join(txt))
print("\n".join(txt))

# OPC link: any OPC HF/DHF (or S) low fix within 600 nm of the cape box max within +-12 h
op = []
for b in ("Atl",):
    d = pd.read_csv(os.path.join(HERE, "..", "..", "..", "data", "hf_lows", f"HF_Data_-_{b}.csv"))
    d["t"] = pd.to_datetime(d.date.astype(str), format="%Y%m%d%H", errors="coerce"); op.append(d.dropna(subset=["t"]))
op = pd.concat(op)


def nm(la0, lo0, la, lo):
    la0, lo0, la, lo = map(np.radians, (la0, lo0, la, lo))
    return 3440.065 * np.arccos(np.clip(np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0), -1, 1))


lines = ["OPC archive low fixes (Atlantic) within 600 nm of the gust maximum, +-12 h; 'OPC agreed' = parent low in an OPC fix, not a check of the tip jet itself"]
res = []
for _, r in pk.iterrows():
    w = op[(op.t >= r.time - pd.Timedelta(hours=12)) & (op.t <= r.time + pd.Timedelta(hours=12))]
    w = w[nm(r.max_lat, r.max_lon, w.Latitude.values, w.Longitude.values) <= 600]
    res.append((w.Category.isin(["HF", "DHF"]).any(), w.Category.isin(["HF", "DHF", "S"]).any()))
pk["opc_hf"] = [x[0] for x in res]; pk["opc_s_or_hf"] = [x[1] for x in res]
g = pk.groupby(["cls", "hf_ep"]).agg(n=("ep", "size"), opc_hf=("opc_hf", "mean"), opc_storm_or_hf=("opc_s_or_hf", "mean"))
lines.append(g.round(2).to_string())
open(os.path.join(RES, "opc_link" + SUF + ".txt"), "w").write("\n".join(lines))
print("\n".join(lines))
