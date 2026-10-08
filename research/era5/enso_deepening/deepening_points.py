"""Max-deepening point of every track in PR 38's tracker output (ERA5 MSLP proxy; pipeline A detector/linker).

For each fix t with a fix at t+24 h: dP24 = MSL(t) - MSL(t+24 h); B = dP24 * sin 60 / sin(mean lat).
A track's max-deepening window is the one with the highest B (first if tied); its point is the window midpoint.
Also records the fix of minimum pressure. One row per track with >= 5 fixes.

usage: deepening_points.py FIXES_CSV_GZ OUT_CSV_GZ
"""
import sys
import numpy as np, pandas as pd

def main(fixes, out):
    f = pd.read_csv(fixes, dtype={"time": str})
    f["t"] = pd.to_datetime(f.time, format="%Y%m%d%H")
    f = f.sort_values(["tid", "t"]).reset_index(drop=True)
    g = f.groupby("tid")
    for k in ("msl", "lat", "lon", "t"):
        f[k + "2"] = g[k].shift(-4)
    ok = (f.t2 - f.t).dt.total_seconds() == 86400
    w = f[ok].copy()
    w["dp24"] = w.msl - w.msl2
    w["latm"] = (w.lat + w.lat2) / 2
    w["lonm"] = (w.lon + w.lon2) / 2          # Pacific box 120-240E never straddles 0/360
    w["B"] = w.dp24 * np.sin(np.radians(60)) / np.sin(np.radians(w.latm))
    w["tm"] = w.t + pd.Timedelta(hours=12)
    best = w.loc[w.groupby("tid").B.idxmax(), ["tid", "winter", "tm", "latm", "lonm", "B", "dp24"]].copy()
    best["tm"] = best.tm.dt.strftime("%Y%m%d%H")
    mn = f.loc[f.groupby("tid").msl.idxmin(), ["tid", "lat", "lon", "msl", "time"]].rename(
        columns={"lat": "lat_min", "lon": "lon_min", "msl": "msl_min", "time": "t_min"})
    best = best.merge(mn, on="tid")
    best.to_csv(out, index=False)
    print(len(best), "tracks with a deepening window")

if __name__ == "__main__":
    main(*sys.argv[1:3])
