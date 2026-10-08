"""RA-30 step 1: central pressure for tracker V fixes (V tracks carry vorticity only).
For every V fix at 00/12 UTC: minimum 1.5 degree ERA5 (WeatherBench2) MSLP within 500 km of the fix (ERA5 proxy).
usage: python3 vpressure.py V_FIXES.csv.gz OUT.csv.gz         (pulls MSLP 2-day chunks, about 1 GB; resumable cache in work/)
"""
import os, sys, datetime as dt
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "second_tracker"))
import common as C

V, OUT = sys.argv[1], sys.argv[2]
WORK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work", "mslp")
os.makedirs(WORK, exist_ok=True)
d = pd.read_csv(V)
d["time"] = d.time.astype(np.int64)
d["hh"] = d.time % 100
d = d[d.hh.isin([0, 12])].copy()
print("V fixes at 00/12 UTC:", len(d), flush=True)
cks = {}
for t in d.time.unique():
    c, k = C.chunk_of(dt.datetime.strptime(str(t), "%Y%m%d%H"))
    cks.setdefault(c, []).append((int(t), k))
print("chunks", len(cks), flush=True)


def get(c):
    f = f"{WORK}/{c}.npy"
    if os.path.exists(f):
        return np.load(f)
    a = C.fetch("mean_sea_level_pressure", c)
    np.save(f + ".tmp.npy", a.astype(np.float32))
    os.replace(f + ".tmp.npy", f)
    return a


cl = sorted(cks)
with ThreadPoolExecutor(8) as ex:
    arrs = dict(zip(cl, ex.map(get, cl)))
print("pulled bytes", C.BYTES[0], flush=True)
LATD, LON = C.LATD, C.LON
byt = {t: g for t, g in d.groupby("time")}
out = []
for c in cl:
    for t, k in cks[c]:
        g = byt[t]
        f = arrs[c][k]                                    # [lat, lon]
        for tr, la, lo in zip(g.track.values, g.lat.values, g.lon.values):
            i0 = int(round((la - LATD[0]) / C.H)); j0 = int(round(lo / C.H)) % 240
            di = 4; dj = int(np.ceil(4.6 / (C.H * max(np.cos(np.radians(la)), 0.2))))
            ii = np.arange(max(i0 - di, 0), min(i0 + di + 1, len(LATD)))
            jj = (np.arange(j0 - dj, j0 + dj + 1)) % 240
            sub = f[np.ix_(ii, jj)]
            dd = C.dist(la, lo, LATD[ii][:, None], LON[jj][None, :])
            m = dd <= 500.0
            out.append((tr, t, float(sub[m].min())))
r = pd.DataFrame(out, columns=["track", "time", "mslp500"])
r.to_csv(OUT, index=False)
print("rows", len(r), "mslp range", r.mslp500.min(), r.mslp500.max())
