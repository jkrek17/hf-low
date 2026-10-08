"""Pipeline A low detection on ERA5 surface pressure only (no gust), 6-hourly, one CSV per month.

ERA5 proxy. Same detection as `../hf_history/extract.py` (lightly smoothed MSLP minima below 1010 hPa,
20-75 N, refined to the raw 5x5 minimum, 400 km de-duplication) with the gust step removed, so the 1979-2003
lows can be found without the gust field the project bars before 2004.

usage: extract_lows.py OUT_DIR YYYYMM [YYYYMM ...]       (resumable: skips months already written)
       extract_lows.py OUT_DIR SEASONS Y0 Y1             (June Y0 .. May Y1+1, newest month first)
Env: NPROC (default 4). Bytes streamed are printed at the end (about 2.2 MB per time).
"""
import os, sys, time
import numpy as np, numcodecs, gcsfs
from scipy import ndimage
from multiprocessing import Pool

ROOT = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
T0 = np.datetime64("1900-01-01T00")
R0, R1 = 60, 281                    # rows for 75N .. 20N
LAT = 90 - 0.25 * np.arange(R0, R1)
LON = 0.25 * np.arange(1440)
RE = 6371.0
_fs = None
_codec = numcodecs.Blosc(cname="lz4", clevel=5, shuffle=1)
LA2, LO2 = np.meshgrid(np.radians(LAT), np.radians(LON), indexing="ij")
XYZ = np.stack([np.cos(LA2) * np.cos(LO2), np.cos(LA2) * np.sin(LO2), np.sin(LA2)], -1)
NBYTES = 0


def fs():
    global _fs
    if _fs is None:
        _fs = gcsfs.GCSFileSystem(token="anon")
    return _fs


def get(var, ti):
    global NBYTES
    for k in range(6):
        try:
            b = fs().cat(f"{ROOT}/{var}/{ti}.0.0")
            NBYTES += len(b)
            return np.frombuffer(_codec.decode(b), "<f4").reshape(721, 1440)[R0:R1]
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {var} {ti}")


def lows_at(ti):
    msl = get("mean_sea_level_pressure", ti) / 100.0
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


def month(args):
    out, (y, m) = args
    path = f"{out}/{y}{m:02d}.csv"
    if os.path.exists(path):
        return (y, m), 0
    a = np.datetime64(f"{y}-{m:02d}-01T00")
    b = np.datetime64(f"{y + (m == 12)}-{m % 12 + 1:02d}-01T00")
    rows, t = [], a
    global NBYTES
    n0 = NBYTES
    while t < b:
        ti = int((t - T0) / np.timedelta64(1, "h"))
        stamp = str(t).replace("-", "").replace("T", "")[:10]
        for r in lows_at(ti):
            rows.append((stamp,) + r)
        t += np.timedelta64(6, "h")
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write("time,lat,lon,msl\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    os.replace(tmp, path)
    return (y, m), NBYTES - n0


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    if sys.argv[2] == "SEASONS":
        y0, y1 = int(sys.argv[3]), int(sys.argv[4])
        months = [(y, m) for y in range(y0, y1 + 2) for m in range(1, 13) if (y0, 6) <= (y, m) <= (y1 + 1, 5)]
        months.sort(key=lambda ym: (-ym[0], -ym[1]))
    else:
        months = [(int(s[:4]), int(s[4:])) for s in sys.argv[2:]]
    t = time.time()
    tot = 0
    with Pool(int(os.environ.get("NPROC", 4))) as p:
        for i, (ym, nb) in enumerate(p.imap_unordered(month, [(out, ym) for ym in months])):
            tot += nb
            print(f"{ym} {i + 1}/{len(months)} {time.time() - t:.0f}s  streamed {tot / 1e9:.1f} GB", flush=True)
