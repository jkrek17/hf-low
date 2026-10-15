"""Stage 0 (no pull): how often does the Cape Farewell box see a gust of each class?

Reads only the shared gust cache (ERA5 PROXY, 00/12 UTC, 2004-05..2025-26). Writes
results/stage0_box_gust.csv: one row per Nov-Mar time with the box maximum gust, where it is, and the
ocean-cell count. No ingredient and no outcome other than the box gust is looked at here.
usage: stage0_counts.py [CACHE_DIR]
"""
import os, sys, glob
import numpy as np, pandas as pd, numcodecs

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = sys.argv[1] if len(sys.argv) > 1 else "/mnt/project-files/era5-cache"
CF = (59.77, -43.9)            # Cape Farewell (Nunap Isua)
RE = 6371.0
codec = numcodecs.Blosc()
# cache grid: rows 90N..0N (0.25), cols lon 0..359.75
lsm = np.load(os.path.join(CACHE, "lsm_nh_0p25.npy"))
lat = 90 - 0.25 * np.arange(361)
lon = 0.25 * np.arange(1440)
LA, LO = np.meshgrid(lat, lon, indexing="ij")
lo180 = np.where(LO > 180, LO - 360, LO)
d = RE * np.arccos(np.clip(np.sin(np.radians(LA)) * np.sin(np.radians(CF[0])) +
        np.cos(np.radians(LA)) * np.cos(np.radians(CF[0])) * np.cos(np.radians(lo180 - CF[1])), -1, 1))
ocean = lsm[:361] < 0.5
box = ocean & (d <= 400)
rows = np.where(box.any(1))[0]; cols = np.where(box.any(0))[0]
sub = box[np.ix_(rows, cols)]
subd = d[np.ix_(rows, cols)]; sublat = LA[np.ix_(rows, cols)]; sublon = lo180[np.ix_(rows, cols)]
T0 = np.datetime64("1900-01-01T00")
out = []
for f in sorted(glob.glob(os.path.join(CACHE, "gust_nh_12h", "gust_*.npz"))):
    z = np.load(f)
    t = T0 + z["hours_since_1900"].astype("timedelta64[h]")
    a = np.frombuffer(codec.decode(z["data"].tobytes()), np.uint8).reshape(z["shape"])
    m = (t.astype("datetime64[M]").astype(int) % 12)  # 0=Jan
    keep = np.isin(m, [10, 11, 0, 1, 2])
    for i in np.where(keep)[0]:
        w = a[i][np.ix_(rows, cols)].astype(float) / 2
        w = np.where(sub, w, -1)
        k = np.unravel_index(np.argmax(w), w.shape)
        out.append((str(t[i]), w[k], sublat[k], sublon[k], subd[k]))
df = pd.DataFrame(out, columns=["time", "box_gust_kt", "max_lat", "max_lon", "max_dist_cf_km"])
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
df.to_csv(os.path.join(HERE, "results", "stage0_box_gust.csv"), index=False)
print(len(df), "Nov-Mar times")
for thr in (38, 48, 54, 60, 64, 71.7):
    print(f">= {thr}: {(df.box_gust_kt >= thr).sum()} times ({(df.box_gust_kt >= thr).mean():.3%})")
