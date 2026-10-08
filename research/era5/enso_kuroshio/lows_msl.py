"""MSLP-only low detection for Dec-Mar, 6-hourly, 20-75N: pipeline A's detector
(research/era5/hf_history/extract.py::lows_at) with the gust part removed.

One CSV per winter in WORK/lows_msl/, resumable. Columns: time,lat,lon,msl.
usage: lows_msl.py FIRST_WINTER LAST_WINTER [PROCS]
"""
import os, sys
import numpy as np, pandas as pd
from scipy import ndimage
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
from common import get, tindex

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
R0, R1 = 60, 281
LAT = 90 - 0.25 * np.arange(R0, R1)
LON = 0.25 * np.arange(1440)
RE = 6371.0
LA2, LO2 = np.meshgrid(np.radians(LAT), np.radians(LON), indexing="ij")
XYZ = np.stack([np.cos(LA2) * np.cos(LO2), np.cos(LA2) * np.sin(LO2), np.sin(LA2)], -1)


def lows_from(msl_full):
    msl = msl_full[R0:R1] / 100.0
    sm = ndimage.gaussian_filter(msl, 2, mode=["nearest", "wrap"])
    mn = ndimage.minimum_filter(sm, size=(21, 21), mode=["nearest", "wrap"])
    r, c = np.where((sm == mn) & (sm < 1010.0))
    keep = (r > 2) & (r < len(LAT) - 3)
    r, c = r[keep], c[keep]
    pts = []
    for i, j in zip(r, c):
        rr = slice(max(i - 2, 0), i + 3)
        cols = np.arange(j - 2, j + 3) % 1440
        box = msl[rr][:, cols]
        a, b = np.unravel_index(np.argmin(box), box.shape)
        pts.append((box[a, b], max(i - 2, 0) + a, cols[b]))
    pts.sort()
    sel = []
    for p, i, j in pts:
        v = XYZ[i, j]
        if all(np.arccos(np.clip(v @ XYZ[a, b], -1, 1)) * RE > 400 for _, a, b in sel):
            sel.append((p, i, j))
    return [(LAT[i], LON[j], round(float(p), 1)) for p, i, j in sel]


def winter(s):
    path = f"{WORK}/lows_msl/{s}.csv"
    if os.path.exists(path):
        return s
    times = np.arange(np.datetime64(f"{s}-12-01T00"), np.datetime64(f"{s + 1}-04-01T00"), np.timedelta64(6, "h"))
    with ThreadPoolExecutor(6) as ex:
        fields = ex.map(lambda t: get("mean_sea_level_pressure", tindex(t)), times)
        rows = []
        for t, f in zip(times, fields):
            stamp = str(t).replace("-", "").replace("T", "")[:10]
            for la, lo, p in lows_from(f):
                rows.append((stamp, la, lo, p))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(rows, columns=["time", "lat", "lon", "msl"]).to_csv(path + ".tmp", index=False)
    os.replace(path + ".tmp", path)
    return s


if __name__ == "__main__":
    a, b = int(sys.argv[1]), int(sys.argv[2])
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    with Pool(n) as p:
        for s in p.imap_unordered(winter, range(b, a - 1, -1)):
            print("done", s, flush=True)
