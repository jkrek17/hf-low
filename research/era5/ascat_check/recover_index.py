"""R0 (pre-registered, outcome-free): recompute the pipeline A index from the cache at 00/12 UTC sample tracks.
usage: recover_index.py out.csv [max_tracks]"""
import sys, numpy as np, pandas as pd
import era5side as E

out = sys.argv[1]
mx = int(sys.argv[2]) if len(sys.argv) > 2 else None
S = E.sample()
idx = E.cache_index()
ocean = E.lsm_ocean()
S = S[S.h1900.isin(idx)]
if mx:
    S = S.sample(mx, random_state=0)
S["fr"] = [idx[h] for h in S.h1900]
S = S.sort_values("fr")
rec = np.full(len(S), np.nan)
for n, (k, row) in enumerate(S.iterrows()):
    f, r_ = row.fr
    g = E.field(f, r_)
    R, C = E.disc_cells(row.peak_lat, row.peak_lon)
    v = g[R, C] * ocean[R, C]
    rec[n] = v.max() if v.size else np.nan
S["g800_recomputed"] = rec
S[["track", "basin", "season", "peak_time", "peak_lat", "peak_lon", "gust800_kt", "g800_recomputed"]].to_csv(out, index=False)
d = S.g800_recomputed - S.gust800_kt
print("n", len(S), "frac within 1kt", float((d.abs() <= 1).mean()), "frac >= stored-1", float((d >= -1).mean()),
      "median diff", float(d.median()))
