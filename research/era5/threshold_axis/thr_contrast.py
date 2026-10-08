"""RA-23 paired differences between thresholds (68, 71.7, 75 kt) in the PR 14 daily and PR 64 weekly setups. ERA5 PROXY, pipeline A.
Reuses sustained_hf/contrast.py (estimation, season-block bootstrap and permutation, clustered SE) with classes defined by peak gust.
Outcomes: all, t68, t717, t75 (peak gust index >= 68, 71.7, 75 kt). Pairs: t75-t68, t75-t717, t717-t68 (log RR differences).
usage: thr_contrast.py daily|weekly OUT_CSV NPERM NBOOT
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sustained_hf"))
import contrast as K  # noqa: E402

THR = {"t68": 68.0, "t717": 71.7, "t75": 75.0}
K.CLASSES.clear()
K.CLASSES["all"] = lambda hf, n: np.ones(len(hf), bool)
for nm, t in THR.items():
    K.CLASSES[nm] = (lambda hf, n, t=t: n >= t)
K.PAIRS[:] = [("t75", "t68"), ("t75", "t717"), ("t717", "t68")]

design, out, nperm, nboot = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(K.SEED + 23)
rows = []
names = ["t68", "t717", "t75", "all"]
if design == "daily":
    T = K.S.load_tracks(os.path.join(K.ERA, "hf_history/results/all_tracks.csv.gz"))
    T["n_hf"] = T.gust800_kt
    T["hf"] = 1
    I = K.S.indices(K.CPC, K.REPO)
    for basin, idx in (("atl", "NAO"), ("pac", "PNA")):
        M = K.Daily(T, I, basin, idx, names)
        df = K.analyse(M, names, nperm, nboot, rng)
        df.insert(0, "basin", basin); df.insert(1, "design", f"daily {idx}"); df["n_hf"] = str(M.n); rows.append(df)
else:
    Tfull = K.C.load_tracks(first=2003, last=2025)
    Tfull["n_hf"] = Tfull.gust800_kt
    Tfull["hf"] = 1
    for basin in K.C.BASINS:
        M, n = K.proxy_weekly(basin, Tfull)
        df = K.analyse(M, names, nperm, nboot, rng)
        df.insert(0, "basin", basin); df.insert(1, "design", "weekly"); df["n_hf"] = str(n); rows.append(df)
pd.concat(rows).to_csv(out, index=False, float_format="%.6g")
