"""Link detections into tracks. M: whole 1979-2022 record in one pass. V: per season (no linking across seasons).
Writes WORK/out/{M,V}_fixes.csv.gz (track,time,lat,lon,val) for tracks of >= 4 fixes (24 h).
usage: track.py M|V
"""
import os, sys, glob
import pandas as pd
import common as C

WORK = os.environ.get("ERA5_WORK", os.path.join(os.path.dirname(os.path.abspath(__file__)), "work"))
tr = sys.argv[1]
parts = []
off = 0
files = sorted(glob.glob(f"{WORK}/{tr}/*.csv"), key=lambda f: int(os.path.basename(f)[:-4]))
if tr == "M":
    d = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates(["time", "lat", "lon"]).sort_values("time")
    df = C.link(d.values.tolist(), "M")
    out = df[df.track >= 0]
else:
    outs = []
    for f in files:
        d = pd.read_csv(f).sort_values("time")
        df = C.link(d.values.tolist(), "V")
        df = df[df.track >= 0].copy()
        df["track"] += off
        off = df.track.max() + 1 if len(df) else off
        outs.append(df)
    out = pd.concat(outs)
os.makedirs(f"{WORK}/out", exist_ok=True)
out.to_csv(f"{WORK}/out/{tr}_fixes.csv.gz", index=False)
n = out.groupby("track").size()
print(tr, "tracks", len(n), "fixes", len(out), "median fixes", n.median())
