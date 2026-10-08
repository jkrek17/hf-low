"""Pilot sample (no pattern index, no HF label): 60 random 00/12 fixes with catalog g800 60-85 kt, Oct-Apr 2004-05 on.
Sets only the fix-level floor of the main pull (how far below 71.7 kt the percentile cut can fall). Seed 20261008."""
import sys, pandas as pd, numpy as np
root = sys.argv[1]
F = pd.read_csv(root + "/research/era5/intensity/results/fixes_2004.csv.gz", dtype={"time": str})
F["m"] = F.time.str[4:6].astype(int)
F = F[F.m.isin([10, 11, 12, 1, 2, 3, 4]) & F.g800.between(60, 85)]
S = F.sample(60, random_state=20261008)[["track", "time", "lat", "lon", "g800"]]
S.to_csv(sys.argv[2], index=False)
print(len(S), S.time.nunique())
