"""Select the RA-27 sample (no field is read): pipeline A tracks of seasons 2016-17..2021-22 with a 00/12 UTC
in-domain gust index >= FLOOR kt, and the 00/12 UTC fixes of those tracks.
Run: python3 -I research/era5/second_analysis/select.py <repo_root> [floor_kt]
"""
import sys, os, pandas as pd, numpy as np
root = sys.argv[1]; floor = float(sys.argv[2]) if len(sys.argv) > 2 else 55.0
F = pd.read_csv(os.path.join(root, "research/era5/hf_probability_tracks/results/fix_probs_2004.csv.gz"))
T = pd.read_csv(os.path.join(root, "research/era5/hf_history/results/all_tracks.csv.gz"), usecols=["track", "gust800_kt", "season"])
F = F[(F.season >= 2016) & (F.season <= 2021)]
idx = F.groupby("track").g800.max().rename("idx1200")
sel = idx[idx >= floor].index
S = F[F.track.isin(sel)].merge(idx, left_on="track", right_index=True)
out = os.path.join(root, "research/era5/second_analysis/results")
S[["track", "time", "basin", "season", "lat", "lon", "msl", "g800", "idx1200"]].to_csv(os.path.join(out, f"sample_fixes_{int(floor)}.csv"), index=False)
print(f"floor {floor}: tracks {len(sel)}  fixes {len(S)}  unique times {S.time.nunique()}  tracks>=71.7: {(idx>=71.7).sum()}  by basin {S.drop_duplicates('track').basin.value_counts().to_dict()}")
