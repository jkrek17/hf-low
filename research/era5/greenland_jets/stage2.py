"""Stage 2: the remaining ingredients at the cases and the control sample.

ERA5 proxy. For each time in results/sample_times.csv (all cases plus the seeded
control sample, PREREGISTRATION.md) it records

  STAB     mean (SST - 2 m temperature) in K over ocean points 60-70N, 40W-10W within
           600 km of Greenland (ARCO-ERA5, 0.25 degree)
  GBI      mean 500 hPa height (m, cos-latitude weighted) over 60-80N, 80W-20W, from
           WeatherBench2's 1.5 degree store; only before 2023-01-10 (else NaN)

NAO, the low's motion and its depth are taken from the CPC file and the fix table in
analyse.py. Output (ignored, resumable): $ERA5_WORK/greenland_jets/stage2/<stamp>.csv.

usage: stage2.py [NPROC]
"""
import os, sys, time, urllib.request, json
import numpy as np, pandas as pd, numcodecs
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
OUT = os.path.join(WORK, "stage2")
sys.path.insert(0, HERE)
import stage1 as s1

WB2 = "https://storage.googleapis.com/weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
WB2_T0 = np.datetime64("1959-01-01T00")
WB2_END = np.datetime64("2023-01-10T00")
WB2_LEVELS = (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000)
G = 9.80665
_codec = numcodecs.Blosc()
_meta = None
STABM = None


def get(url):
    for k in range(7):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception:
            if k == 6:
                raise
            time.sleep(2 ** k)


def init():
    global STABM
    s1.init()
    ocean = s1.S["lsm"] < 0.5
    LA2, LO2 = np.meshgrid(s1.LA_FULL, s1.LON, indexing="ij")
    STABM = ocean & (LA2 >= 60) & (LA2 <= 70) & (LO2 >= 320) & (LO2 <= 350) & (s1.S["dgl"] <= 600.0)
    os.makedirs(OUT, exist_ok=True)


def gbi(t):
    global _meta
    if t >= WB2_END:
        return np.nan, 0
    if _meta is None:
        _meta = json.loads(get(f"{WB2}/.zmetadata"))["metadata"]
    za = _meta["geopotential/.zarray"]
    k6 = int((t - WB2_T0) / np.timedelta64(6, "h"))
    chunk, j = divmod(k6, za["chunks"][0])
    raw = get(f"{WB2}/geopotential/{chunk}.0.0.0")
    codec = numcodecs.get_codec(za["compressor"])
    a = np.frombuffer(codec.decode(raw), za["dtype"]).reshape(za["chunks"])
    z = a[j, WB2_LEVELS.index(500)] / G                      # (lon, lat)
    lat = -90 + 1.5 * np.arange(121)
    lon = 1.5 * np.arange(240)
    li = (lat >= 60) & (lat <= 80)
    lo = (lon >= 280) & (lon <= 340)
    w = np.cos(np.radians(lat[li]))[None, :] * np.ones((lo.sum(), 1))
    return float((z[lo][:, li] * w).sum() / w.sum()), len(raw)


def one(stamp):
    path = f"{OUT}/{stamp}.csv"
    if os.path.exists(path):
        return stamp
    t = np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}")
    ti = int((t - s1.T0) / np.timedelta64(1, "h"))
    sst, n1 = s1.fetch("sea_surface_temperature", ti)
    t2m, n2 = s1.fetch("2m_temperature", ti)
    d = (sst - t2m)[STABM]
    stab = float(np.nanmean(d))
    nsst = int(np.isfinite(d).sum())
    gb, n3 = gbi(t)
    pd.DataFrame([dict(time=stamp, STAB=round(stab, 3), n_stab_pts=nsst, n_stab_total=int(STABM.sum()),
                       GBI=round(gb, 2) if gb == gb else np.nan, bytes=n1 + n2 + n3)]).to_csv(path + ".tmp", index=False)
    os.replace(path + ".tmp", path)
    return stamp


if __name__ == "__main__":
    tg = [str(x) for x in pd.read_csv(os.path.join(HERE, "results", "sample_times.csv")).time]
    print(len(tg), "times", flush=True)
    t0 = time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8, initializer=init) as p:
        for i, s in enumerate(p.imap_unordered(one, tg)):
            if i % 100 == 0:
                print(s, i + 1, f"{time.time() - t0:.0f}s", flush=True)
