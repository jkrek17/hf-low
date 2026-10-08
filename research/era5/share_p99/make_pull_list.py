"""List the 00/12 UTC fixes to pull for RA-28 (no ERA5 access, no pattern index, no HF label).

A fix is listed if its catalog g800 (the maximum owned-cell gust, so an upper bound for any percentile of the same cells) is at least FLOOR kt
and it belongs to (a) a track in the share-mechanism population (first fix in Oct 1 + 210 d of 2004..2025, as chanlib.assign_week), or
(b) the calibration seasons 2021-22..2025-26 (June to May). usage: make_pull_list.py REPO OUT_CSV [FLOOR=60]
"""
import sys, os, numpy as np, pandas as pd
root, out = sys.argv[1], sys.argv[2]
FLOOR = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0
sys.path.insert(0, os.path.join(root, "research/era5/hem_channels"))
import chanlib as C
T = C.assign_week(C.load_tracks(first=2003, last=2025))
pop = set(T.track[T.basin.isin(C.BASINS)])
F = pd.read_csv(C.FIXES, dtype={"time": str})
F = F[F.basin.isin(C.BASINS) & (F.g800 >= FLOOR)]
cal = F.season.between(2021, 2025)
L = F[F.track.isin(pop) | cal]
L[["track", "time", "lat", "lon", "g800"]].to_csv(out, index=False)
print("fixes", len(L), "times", L.time.nunique(), "tracks", L.track.nunique(), "floor", FLOOR)
print("population tracks", len(pop), "| est GB at 5.36 MB per time:", round(L.time.nunique() * 5.36e-3, 1))
