"""Pre-registered tracker check (b): the committed 1979-2003 catalog track points against the pressure-only
re-detection and re-tracking. usage: check_tracker.py LOWS_DIR FIXES_1979_2003.csv.gz"""
import sys, glob
import pandas as pd
lows = pd.concat([pd.read_csv(f, dtype={"time": str}) for f in glob.glob(sys.argv[1] + "/*.csv")])
F = pd.read_csv(sys.argv[2], dtype={"time": str})
C = pd.read_csv("../hf_history/results/era5_hf_catalog_tracks.csv", dtype={"time": str})
C = C[C.time < "2004060100"]
m = C.merge(lows, on="time", how="left", suffixes=("", "_r"))
m = m[m.lat_r.isna() | ((m.lat - m.lat_r).abs() < 0.01) & ((m.lon - m.lon_r).abs() < 0.01)].drop_duplicates(["track", "time"])
found = m.lat_r.notna()
same = found & ((m.msl - m.msl_r).abs() <= 0.2)
print(f"catalog points before 2004-06: {len(C)}; low re-detected at same time and position {found.sum()} ({found.mean():.4%}); "
      f"pressure within 0.2 hPa {same.sum()} ({same.mean():.4%})")
d = C.assign(h=C.time.str[8:]).query("h in ['00','12'] and basin in ['atl','pac']")
j = d.merge(F[["time", "lat", "lon", "msl", "track", "basin"]], on=["time", "lat", "lon"], how="left", suffixes=("", "_n"))
j = j.drop_duplicates(["track", "time"])
print(f"in-domain 00/12 catalog points {len(d)}; found among forecast fixes {j.track_n.notna().sum()} ({j.track_n.notna().mean():.4%}); "
      f"same basin {(j.basin == j.basin_n).sum()}")
g = j.dropna(subset=["track_n"]).groupby("track").track_n.nunique()
print(f"catalog tracks with in-domain points {d.track.nunique()}; mapped to exactly one re-linked track {(g == 1).sum()}; to several {(g > 1).sum()}")
miss = j[j.track_n.isna()]
print("missing in-domain points by year:", miss.time.str[:4].value_counts().sort_index().to_dict())
print((m[~found].time.str[:4].value_counts().sort_index().to_dict()), "<- catalog points whose low was not re-detected, by year")
