"""RA-28 track labels at 00/12 UTC fixes: M (maximum owned-cell gust), P99, P98, each at its calibrated cut (results/cuts.json).

A track's index is the maximum over its in-domain 00/12 UTC fixes in intensity/results/fixes_2004.csv.gz; for the percentiles, fixes with
catalog g800 below the pull floor were not pulled and count as 0. Tracks with no 00/12 fix are not in the table (the caller fills 0).
ERA5 PROXY, pipeline A.
"""
import json, os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def labels(repo):
    cuts = json.load(open(os.path.join(HERE, "results", "cuts.json")))["cuts_kt"]
    F = pd.read_csv(os.path.join(repo, "research/era5/intensity/results/fixes_2004.csv.gz"), usecols=["track", "time", "g800"], dtype={"time": str})
    X = pd.read_csv(os.path.join(HERE, "results", "fix_p99.csv.gz"), dtype={"time": str}, usecols=["track", "time", "p98", "p99"])
    F = F.merge(X, on=["track", "time"], how="left")
    g = pd.DataFrame({"M": F.groupby("track").g800.max(), "P99": F.assign(v=F.p99.fillna(0.0)).groupby("track").v.max(),
                      "P98": F.assign(v=F.p98.fillna(0.0)).groupby("track").v.max()})
    return {k: (g[k] >= cuts[k]).astype(int) for k in g}, cuts
