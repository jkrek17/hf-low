"""Stage 0 (no pull): how often does each Gulf of Alaska zone see a gust of each class?

Reads only the shared ERA5 gust cache (PROXY, 00/12 UTC, 2004-05..2025-26). Writes results/stage0_zone_gust.csv: one row
per October-April time with, for each zone, the largest ocean-cell gust and where it is. No wind direction, pressure or
other ingredient, and no outcome other than the zone gust, is looked at here.
usage: stage0_counts.py [CACHE_DIR]
"""
import os, sys, glob
import numpy as np, pandas as pd, numcodecs
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = sys.argv[1] if len(sys.argv) > 1 else "/mnt/project-files/era5-cache"
# zone: (south, north, west, east) in degrees; west/east as degrees west (positive)
ZONES = {"Z1": (56, 61, 150, 138), "Z2": (53, 58, 140, 130), "Z3": (55, 60, 158, 148)}
codec = numcodecs.Blosc()
lsm = np.load(os.path.join(CACHE, "lsm_nh_0p25.npy"))
T0 = np.datetime64("1900-01-01T00")


def window(z):
    s, n, w, e = ZONES[z]
    j0, j1 = int((90 - n) / .25), int((90 - s) / .25)
    k0, k1 = int((360 - w) / .25), int((360 - e) / .25)
    return j0, j1, k0, k1


def one(f):
    z = np.load(f)
    t = T0 + z["hours_since_1900"].astype("timedelta64[h]")
    a = np.frombuffer(codec.decode(z["data"].tobytes()), np.uint8).reshape(z["shape"])
    m = t.astype("datetime64[M]").astype(int) % 12
    keep = np.where(np.isin(m, [9, 10, 11, 0, 1, 2, 3]))[0]
    rows = []
    for i in keep:
        r = [str(t[i])]
        for zn in ZONES:
            j0, j1, k0, k1 = window(zn)
            g = a[i, j0:j1, k0:k1].astype(float) / 2
            g = np.where(lsm[j0:j1, k0:k1] < 0.5, g, -1)
            jj, kk = np.unravel_index(np.argmax(g), g.shape)
            r += [g[jj, kk], 90 - 0.25 * (j0 + jj), (k0 + kk) * 0.25 - 360]
        rows.append(r)
    return rows


if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(CACHE, "gust_nh_12h", "gust_*.npz")))
    with Pool(4) as p:
        out = [r for rr in p.imap(one, files) for r in rr]
    cols = ["time"] + [f"{zn}_{c}" for zn in ZONES for c in ("gust_kt", "lat", "lon")]
    df = pd.DataFrame(out, columns=cols).sort_values("time")
    df.to_csv(os.path.join(HERE, "results", "stage0_zone_gust.csv"), index=False)
    print(len(df), "Oct-Apr times")
    ev = pd.Series(False, index=df.index)
    for zn in ZONES:
        print(zn, {thr: int((df[f"{zn}_gust_kt"] >= thr).sum()) for thr in (38, 50, 54, 64, 71.7)})
        ev |= df[f"{zn}_gust_kt"] >= 54
    print("union >=54 in any zone:", int(ev.sum()))
