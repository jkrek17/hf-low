"""Link the MSLP-only lows of each winter into tracks with pipeline A's linker and
write one row per track (genesis = first fix) plus a fix table.

usage: tracks.py WORK_DIR OUT_DIR
"""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_history"))
import track as A            # pipeline A linker (link, MIN_LEN, basin)
from common import LAT, LON

def main(work, out):
    lsm = np.load(os.path.join(work, "lsm.npy"))
    trk, fix = [], []
    for f in sorted(glob.glob(os.path.join(work, "lows_msl", "*.csv"))):
        s = int(os.path.basename(f)[:4])
        df = pd.read_csv(f, dtype={"time": str})
        df = df[(df.lat >= 20) & (df.lat <= 75)]
        df = A.link(df)
        n = df.groupby("track").size()
        df = df[df.track.isin(n[n >= A.MIN_LEN].index)].sort_values(["track", "time"]).copy()
        df["winter"] = s
        df["tid"] = s * 100000 + df.track
        fix.append(df)
        g = df.groupby("tid")
        first = g.head(1).set_index("tid")
        t = pd.DataFrame({"winter": s, "t0": first.time, "lat0": first.lat, "lon0": first.lon,
                          "msl0": first.msl, "n_fix": g.size(), "minp": g.msl.min(),
                          "t_end": g.time.last()})
        ri = np.clip(np.round((90 - t.lat0.values) / 0.25).astype(int), 0, 720)
        ci = np.round(t.lon0.values / 0.25).astype(int) % 1440
        t["ocean0"] = lsm[ri, ci] < 0.5
        trk.append(t.reset_index())
    T = pd.concat(trk, ignore_index=True)
    F = pd.concat(fix, ignore_index=True)[["tid", "winter", "time", "lat", "lon", "msl"]]
    os.makedirs(out, exist_ok=True)
    T.to_csv(os.path.join(out, "tracks.csv.gz"), index=False)
    F.to_csv(os.path.join(out, "fixes.csv.gz"), index=False)
    print(len(T), "tracks", len(F), "fixes")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
