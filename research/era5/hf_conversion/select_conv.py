"""Test 1 stage B selection (labels and positions only): random N pairs per basin from results/pairs.csv; times t0, t0-12, t0-24 h per member.
usage: select_conv.py REPO_ROOT N. Position at a lag = the 00/12 UTC fix on the same track (heading = heading at t0). Writes results/times_conv.csv, results/conv_sample.csv."""
import sys, os, numpy as np, pandas as pd
root, N = sys.argv[1], int(sys.argv[2]); E5 = os.path.join(root, "research/era5"); OUT = os.path.join(E5, "hf_conversion/results")
P = pd.read_csv(os.path.join(OUT, "pairs.csv"))
keep = []
for b, g in P.groupby("basin"):
    pr = g.pair.drop_duplicates().sample(frac=1, random_state=20261012).iloc[:N]; keep.append(g[g.pair.isin(pr)])
S = pd.concat(keep)
F = pd.read_csv(os.path.join(E5, "intensity/results/fixes_2004.csv.gz"), usecols=["track", "time", "lat", "lon"]).set_index(["track", "time"])
rows = []
for r in S.itertuples():
    t0 = pd.to_datetime(str(r.time), format="%Y%m%d%H")
    for L in (0, 12, 24):
        tl = int((t0 - pd.Timedelta(hours=L)).strftime("%Y%m%d%H"))
        if (r.track, tl) in F.index:
            f = F.loc[(r.track, tl)]
            rows.append(dict(track=r.track, anchor=f"g{r.grp}_p{r.pair}_L{L}", time=tl, lat=float(f.lat), lon=float(f.lon), heading=r.heading, basin=r.basin))
T = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, "conv_sample.csv"), index=False); T.to_csv(os.path.join(OUT, "times_conv.csv"), index=False)
print(f"N={N}: pairs {S.pair.nunique()} (members {len(S)}), lag rows {len(T)}, distinct times {T.time.nunique()}")
