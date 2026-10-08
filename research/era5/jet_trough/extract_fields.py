"""Z500, U250 and V250 at 00 and 12 UTC on the 5.625 degree grid, Northern Hemisphere 14-87 N.

ERA5 proxy. Outcome-free: nothing here reads a track, a fix or an outcome.

Sources, as `../hemispheric/fields.py` (whose readers and conservative regrid are reused):
  <= 2023-01-09  WeatherBench2 64 x 32 conservative regrid of ERA5, 6-hourly
  >= 2023-01-10  ARCO-ERA5 hourly 0.25 degree, regridded onto the same cells; only the blosc
                 blocks of the one level wanted are fetched (HTTP range request).

usage: extract_fields.py OUT_DIR           (resumable: one npz per season)
Writes season npz files with time (hours since 1900-01-01), z500 (m), u250, v250 (m/s) of shape
[ntime, 14, 64]; rows are cell centres 14.06 .. 87.19 N, longitudes 0, 5.625, ...
Volume: WB2 about 0.5 MB per 100-time chunk per variable; ARCO about 1-2 MB per level per time.
Expected total under 10 GB, nothing but OUT_DIR (ignored) is written.
"""
import os, sys, time, datetime as dt
import numpy as np, numcodecs
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hemispheric"))
import fields as F

NH = np.arange(18, 32)                       # cell centres 14.06 .. 87.19 N
F.NH = NH
F.WLAT = F._lat_weights()[NH]
G = F.G
H0_WB2 = dt.datetime(1959, 1, 1)
H0_ARCO = dt.datetime(1900, 1, 1)
NAMES = {"z500": ("geopotential", 500), "u250": ("u_component_of_wind", 250),
         "v250": ("v_component_of_wind", 250)}
CH = 100                                     # WB2 64 x 32 time chunk


def wb2_chunk(var, c):
    name, lev = NAMES[var]
    m = F.meta()
    za = m[name + "/.zarray"]
    codec = numcodecs.get_codec(za["compressor"])
    raw = F.retry(F.fs().cat, f"{F.WB2}/{name}/{c}.0.0.0")
    a = np.frombuffer(codec.decode(raw), za["dtype"]).reshape(za["chunks"])
    nt = min(za["chunks"][0], za["shape"][0] - c * za["chunks"][0])
    a = a[:nt, F.WB2_LEVELS.index(lev)].transpose(0, 2, 1)[:, NH]     # time, lat, lon
    return (a / G if var == "z500" else a).astype(np.float32), len(raw)


def wb2_times(hours):
    """hours: int array of hours since 1900-01-01 (all 00/12 UTC). -> dict var [n, 14, 64], bytes."""
    off = int((H0_WB2 - H0_ARCO).total_seconds() // 3600)
    idx = (np.asarray(hours) - off) // 6
    out = {v: np.full((len(idx), len(NH), 64), np.nan, np.float32) for v in NAMES}
    chunks = sorted({int(i) // CH for i in idx})
    nb = 0
    def job(vc):
        v, c = vc
        a, n = wb2_chunk(v, c)
        return v, c, a, n
    with ThreadPoolExecutor(8) as ex:
        for v, c, a, n in ex.map(job, [(v, c) for v in NAMES for c in chunks]):
            nb += n
            for j, i in enumerate(idx):
                if int(i) // CH == c:
                    out[v][j] = a[int(i) - c * CH]
    return out, nb


def arco_one(h):
    res, nb = {}, 0
    for v, (name, lev) in NAMES.items():
        a, n = F.arco_level(name, int(h), lev)
        res[v] = F.regrid(a / G if v == "z500" else a)
        nb += n
    return res, nb


def arco_times(hours):
    out = {v: np.full((len(hours), len(NH), 64), np.nan, np.float32) for v in NAMES}
    nb = 0
    with ThreadPoolExecutor(8) as ex:
        for j, (r, n) in enumerate(ex.map(arco_one, hours)):
            nb += n
            for v in NAMES:
                out[v][j] = r[v]
    return out, nb


def season_hours(s):
    t0 = dt.datetime(s, 6, 1)
    t1 = dt.datetime(s + 1, 6, 1)
    n = int((t1 - t0).total_seconds() // 3600 // 12)
    return np.array([int((t0 - H0_ARCO).total_seconds() // 3600) + 12 * k for k in range(n)])


if __name__ == "__main__":
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    cut = int((dt.datetime(2023, 1, 10) - H0_ARCO).total_seconds() // 3600)
    tot = 0
    for s in range(2004, 2026):
        path = f"{outdir}/fields_{s}.npz"
        if os.path.exists(path):
            continue
        t0 = time.time()
        hrs = season_hours(s)
        if s == 2025:                           # fixes end 2026-05-31 12 UTC
            pass
        w = hrs < cut
        parts, nb = [], 0
        if w.any():
            o, n = wb2_times(hrs[w]); parts.append((np.where(w)[0], o)); nb += n
        if (~w).any():
            o, n = arco_times(hrs[~w]); parts.append((np.where(~w)[0], o)); nb += n
        arr = {v: np.full((len(hrs), len(NH), 64), np.nan, np.float32) for v in NAMES}
        for ii, o in parts:
            for v in NAMES:
                arr[v][ii] = o[v]
        tmp = path + ".tmp.npz"
        np.savez_compressed(tmp, time=hrs, **arr)
        os.replace(tmp, path)
        tot += nb
        print(f"{s} {len(hrs)} times  {nb / 1e9:.2f} GB  {time.time() - t0:.0f}s  nan z500 {np.isnan(arr['z500']).mean():.4f}",
              flush=True)
    print(f"streamed about {tot / 1e9:.1f} GB this run")
