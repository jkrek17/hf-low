"""Threshold axis for RA-23 (ERA5 proxy, pipeline A). Plan: PREREGISTRATION.md.

Writes, in work/ (ignored), copies of all_tracks.csv.gz in which a track's gust800_kt is 99.0 if its peak gust index is at or above T
and 0.0 otherwise, for T = 68, 71.7, 75 kt, so the unchanged PR 14 / 63 / 64 scripts (which cut at 71.7) see HF_T. Demoted tracks stay cyclones.
T = 71.7 reproduces the original flags exactly (control). usage: labels.py WORKDIR
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
THRESHOLDS = {"t68": 68.0, "t717": 71.7, "t75": 75.0}

if __name__ == "__main__":
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    T = pd.read_csv(os.path.join(ERA, "hf_history/results/all_tracks.csv.gz"), dtype={"start": str, "peak_time": str})
    orig = T.gust800_kt >= 71.7
    for name, thr in THRESHOLDS.items():
        x = T.copy()
        x["gust800_kt"] = np.where(T.gust800_kt >= thr, 99.0, 0.0)
        x.to_csv(os.path.join(work, f"all_tracks_{name}.csv.gz"), index=False)
        if name == "t717":
            assert ((x.gust800_kt >= 71.7) == orig).all()
    w = T[T.season.between(2004, 2025)]
    for b in ("atl", "pac"):
        x = w[w.basin == b]
        print(b, {n: int((x.gust800_kt >= t).sum()) for n, t in THRESHOLDS.items()}, "tracks", len(x))
