"""Test 2: in-scope HF storms, their weekly pattern index and the HF-onset anchors still to be extracted (labels only).
usage: select_storms.py REPO_ROOT. Writes results/hf_index.csv (every HF storm: week, index, tercile) and results/times_pp.csv (onsets not yet extracted by hf_vs_storm)."""
import sys, os, datetime as dt, numpy as np, pandas as pd
root = sys.argv[1]
E5 = os.path.join(root, "research/era5"); OUT = os.path.join(E5, "hf_pattern_phase/results"); os.makedirs(OUT, exist_ok=True)
S = pd.read_csv(os.path.join(E5, "hf_vs_storm/results/storms.csv")); S = S[S.grp == "HF"].copy()
d0 = pd.to_datetime(S.start.astype(str).str[:8], format="%Y%m%d")
oct1 = pd.to_datetime(S.season.astype(str) + "-10-01")
S["week"] = ((d0 - oct1).dt.days // 7)
S["in_scope"] = (S.week >= 0) & (S.week < 30) & S.onset_time.notna()
rows = []
for b in ("atl", "pac"):
    I = pd.read_csv(os.path.join(E5, f"hemispheric/results/oos_index_{b}.csv"))
    I["z"] = (I.idx - I.idx.mean()) / I.idx.std(ddof=0)
    lo, hi = np.percentile(I.z, [100 / 3, 200 / 3])
    m = S.basin == b
    X = S[m].merge(I[["season", "week", "z"]], on=["season", "week"], how="left")
    X["tercile"] = np.where(X.z <= lo, "bottom", np.where(X.z >= hi, "top", "middle"))
    X.loc[X.z.isna() | ~X.in_scope, "tercile"] = "out"
    X["cut_lo"], X["cut_hi"] = lo, hi
    rows.append(X)
H = pd.concat(rows)
H[["track", "basin", "season", "mon", "week", "z", "tercile", "cut_lo", "cut_hi", "in_scope", "onset_time", "onset_lat", "onset_lon", "hd_onset", "peak_time"]].to_csv(os.path.join(OUT, "hf_index.csv"), index=False)
print(H.groupby(["basin", "tercile"]).size().unstack())
done = pd.read_csv(os.path.join(E5, "hf_vs_storm/results/times_storm.csv"))
done = set(done[done.anchor == "HF_onset"].track)
T = H[(H.tercile != "out") & ~H.track.isin(done)]
T = pd.DataFrame(dict(track=T.track, anchor="HF_onset", time=T.onset_time.astype("int64"), lat=T.onset_lat, lon=T.onset_lon, heading=T.hd_onset, basin=T.basin))
T.to_csv(os.path.join(OUT, "times_pp.csv"), index=False)
print("to extract:", len(T), "storms,", T.time.nunique(), "distinct times; already extracted:", len(done))
