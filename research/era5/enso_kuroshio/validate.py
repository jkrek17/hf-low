"""Pre-registered tracker validation: fraction of pipeline A Pacific tracks (peak_time in
the pulled windows) whose peak fix is a fix of a new track. Also writes the HF label map.

usage: validate.py ALL_TRACKS OUT_DIR
"""
import sys, os
import numpy as np, pandas as pd

A = pd.read_csv(sys.argv[1]); out = sys.argv[2]
F = pd.read_csv(os.path.join(out, "fixes.csv.gz"), dtype={"time": str})
A["pt"] = A.peak_time.astype(str)
md = A.pt.str[4:8].astype(int)
inwin = ((md >= 1203) | (md <= 328)) & A.pt.str[:4].astype(int).between(1979, 2026)
# winter of the peak: Dec of s .. Mar of s+1
A["winter"] = np.where(A.pt.str[4:6].astype(int) == 12, A.pt.str[:4].astype(int), A.pt.str[:4].astype(int) - 1)
B = A[inwin & (A.basin == "pac") & (A.winter <= 2025)].copy()
key = F.assign(k=F.time + "_" + F.lat.round(2).astype(str) + "_" + F.lon.round(2).astype(str)).set_index("k").tid
B["k"] = B.pt + "_" + B.peak_lat.round(2).astype(str) + "_" + B.peak_lon.round(2).astype(str)
B["tid"] = B.k.map(key[~key.index.duplicated()])
rate = B.tid.notna().mean()
print(f"pipeline A Pacific tracks with peak in window: {len(B)}; recovered by exact peak fix: {B.tid.notna().sum()} ({rate:.4f})")
print(B.groupby(pd.cut(B.winter, [1978, 1990, 2000, 2010, 2025])).tid.apply(lambda s: s.notna().mean()).round(4).to_string())
B[["track", "winter", "tid", "gust800_kt", "minp", "peak_time", "peak_lat", "peak_lon"]].to_csv(os.path.join(out, "pipelineA_match.csv"), index=False)
