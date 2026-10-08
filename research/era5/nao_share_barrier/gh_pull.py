"""Daily Greenland-high index GH for the NAO-share study (ERA5 PROXY).

GH = mean MSLP (hPa) over Greenland ice-sheet land points 64-78N, 55W-30W, as PR 58
(greenland_jets/stage1.py): ARCO-ERA5 0.25 degree, the land component containing
72N 40W, land-sea mask >= 0.5. One value per day at 12 UTC, 21 Sep to 28 Apr of the
seasons 2004-05 to 2025-26 (lag window needs 10 days before 1 Oct).

Pull: one 2.2 MB MSLP chunk per day, about 4,850 days, about 11 GB streamed, nothing
stored except the 12-value-per-day output. Resumable.

usage: gh_pull.py [NPROC]   writes gh_daily.csv next to this file (small, committed)
"""
import os, sys, time, urllib.request
import numpy as np, pandas as pd, numcodecs
from scipy import ndimage
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
T0 = np.datetime64("1900-01-01T00")
_codec = numcodecs.Blosc()
GHM = None


def fetch(var, ti):
    url = f"{ARCO}/{var}/{ti}.0.0"
    for k in range(7):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            return np.frombuffer(_codec.decode(raw), "<f4").reshape(721, 1440)
        except Exception:
            if k == 6:
                raise
            time.sleep(2 ** k)


def init():
    global GHM
    lsm = fetch("land_sea_mask", 876576)            # 2000-01-01 00Z, as pipeline A
    lab, _ = ndimage.label(lsm >= 0.5)
    gl = lab == lab[int(round((90 - 72.0) / 0.25)), int(round(320.0 / 0.25))]
    assert gl.sum() > 1000
    la, lo = np.meshgrid(90 - 0.25 * np.arange(721), 0.25 * np.arange(1440), indexing="ij")
    GHM = gl & (la >= 64.0) & (la <= 78.0) & (lo >= 305.0) & (lo <= 330.0)


def one(day):
    path = os.path.join(WORK, "gh", f"{day}.txt")
    if os.path.exists(path):
        return day
    t = np.datetime64(f"{day[:4]}-{day[4:6]}-{day[6:8]}T12")
    msl = fetch("mean_sea_level_pressure", int((t - T0) / np.timedelta64(1, "h"))) / 100.0
    with open(path + ".tmp", "w") as f:
        f.write(f"{float(msl[GHM].mean()):.3f} {int(GHM.sum())}\n")
    os.replace(path + ".tmp", path)
    return day


if __name__ == "__main__":
    os.makedirs(os.path.join(WORK, "gh"), exist_ok=True)
    days = []
    for s in range(2004, 2026):
        days += [d.strftime("%Y%m%d") for d in pd.date_range(f"{s}-09-21", f"{s + 1}-04-28")]
    print(len(days), "days", flush=True)
    t0 = time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8, initializer=init) as p:
        for i, d in enumerate(p.imap_unordered(one, days, chunksize=4)):
            if i % 250 == 0:
                print(d, i + 1, f"{time.time() - t0:.0f}s", flush=True)
    rows = []
    for d in days:
        v, n = open(os.path.join(WORK, "gh", f"{d}.txt")).read().split()
        rows.append((pd.Timestamp(d).strftime("%Y-%m-%d"), float(v), int(n)))
    pd.DataFrame(rows, columns=["date", "GH", "n_pts"]).to_csv(os.path.join(HERE, "gh_daily.csv"), index=False)
    print("done", flush=True)
