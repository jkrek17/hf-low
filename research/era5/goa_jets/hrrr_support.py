"""Supporting numbers quoted in README_hrrr.md that the main analysis files do not hold (added after the claim review): cell counts, HRRR-only events near land,
top-5% band counts per season, buoy distance to HRRR land, HRRR domain share of the reference-low box. Writes results/hrrr_support.txt. usage: hrrr_support.py"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hrrr_analyse as A, hrrr_grid as HG, buoy as B
RES = A.RES; out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
la, lo, land = A.load_grid(); dl, tree = A.dland_grid(la, lo, land)
f = pd.read_csv(os.path.join(RES, "hrrr_features.csv"), parse_dates=["time"])
s0 = pd.read_csv(os.path.join(RES, "stage0_zone_gust.csv"), parse_dates=["time"])
e = pd.concat([s0[["time", f"{z}_gust_kt"]].set_axis(["time", "era5_gust"], axis=1).assign(zone=z) for z in A.ZONES])
f = f.merge(e, on=["time", "zone"]); f["regime"] = [A.regime(z, w) for z, w in zip(f.zone, f.wdir)]; f["season"] = [A.season(t) for t in f.time]
f["h54"] = f.hrrr_gust >= 54; f["e54"] = f.era5_gust >= 54
# cell counts
lsm = np.load("/mnt/project-files/era5-cache/lsm_nh_0p25.npy")[80:185, 740:941] < 0.5
ela = 90 - 0.25 * (80 + np.arange(105)); elo = 0.25 * (740 + np.arange(201)) - 360; ELA, ELO = np.meshgrid(ela, elo, indexing="ij")
P("ocean cells in the zone box, HRRR (3 km) vs ERA5 (0.25 degree):")
for zn, (s, n, w, e_) in A.ZONES.items():
    h = int(((la >= s) & (la <= n) & (lo >= w) & (lo <= e_) & ~land).sum()); r = int(((ELA >= s) & (ELA <= n) & (ELO >= w) & (ELO <= e_) & lsm).sum())
    P(f"  {zn}: HRRR {h}, ERA5 {r}, ratio {h / r:.1f}")
# HRRR-only events near land
P("\nHRRR-only events (HRRR zone gust >= 54 kt and ERA5 < 54 kt): share with the HRRR peak cell within 10 / 25 km of land")
for zn in A.ZONES:
    for lab, m in (("all flows", f.zone == zn), ("barrier/gap flow", (f.zone == zn) & f.regime.isin(["BARRIER", "GAP"]))):
        x = f[m & f.h54 & ~f.e54]
        P(f"  {zn} {lab}: n {len(x)}, within 10 km {np.mean(x.g_dland_km <= 10):.3f}, within 25 km {np.mean(x.g_dland_km <= 25):.3f}, median {x.g_dland_km.median():.1f} km")
# 64 kt ratio by flow
h1 = pd.read_csv(os.path.join(RES, "hrrr_h1.csv")); h1["r64"] = h1.hrrr_total64 / h1.era5_total64
P("\nHRRR/ERA5 count of 64 kt times by zone and flow:"); P(h1[["zone", "regime", "hrrr_total64", "era5_total64", "r64"]].round(1).to_string(index=False))
# top-5% band per season
P("\nTop 5% band times per season and events (HRRR 54 kt):")
for zn in A.ZONES:
    ing = "INTERIOR_DP" if zn == "Z3" else "COASTPILE"; rg = "GAP" if zn == "Z3" else "BARRIER"
    d = f[(f.zone == zn) & (f.regime == rg)]; cut = np.nanpercentile(d[ing], 95); x = d[d[ing] >= cut]
    g = x.groupby("season").h54.agg(["size", "sum"]); P(f"  {zn} {rg} {ing}: n {len(x)}; per season (times, events): " + "; ".join(f"{s}: {int(r['size'])}/{int(r['sum'])}" for s, r in g.iterrows()))
# buoys: distance to HRRR land
P("\nBuoy distance to the nearest HRRR land cell (km):")
from hrrr_analyse import xyz
for s_, (bla, blo) in B.ST.items():
    P(f"  {s_}: {tree.query(xyz(np.array([bla]), np.array([blo])))[0][0]:.1f}")
# domain share of the reference-low box
la0 = np.linspace(44, 62, 91); lo0 = np.linspace(-175, -125, 201); LA, LO = np.meshgrid(la0, lo0, indexing="ij"); j, i = HG.ij(LA, LO)
P(f"\nshare of the reference-low box (44-62N, 175-125W) inside the HRRR-Alaska domain: {np.mean((j >= 0) & (j <= 918) & (i >= 0) & (i <= 1298)):.3f}")
open(os.path.join(RES, "hrrr_support.txt"), "w").write("\n".join(out)); print("\n".join(out))
