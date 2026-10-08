"""Copy of hemispheric/fields.py for the HF vs storm-force composites (hf_vs_storm): adds column water vapour (tcwv),
seasons 2004-05 to 2025-26 only, and daily dates 20 August to 31 May so lags of 10 days reach back from September anchors.
Original docstring follows.

Hemispheric daily (12 UTC) fields for the data-driven pattern search (ERA5 proxy).

Fields: Z500 (m), U250 (m/s), MSLP (hPa), SST (K, ocean only) on WeatherBench2's
64 x 32 conservative grid (5.625 deg), NH rows only (cell centres 25.3125 .. 87.1875 N,
12 rows) x 64 longitudes.

Sources
  <= 2023-01-09  WeatherBench2 1959-2023_01_10-6h-64x32_equiangular_conservative (6-hourly, 12 UTC taken)
  >= 2023-01-10  ARCO-ERA5 hourly 0.25 degree, regridded here onto the same cells by area-weighted
                 (conservative) overlap. The 3-D ARCO chunks hold all 37 levels (80 MB); one level
                 is about 9 of 294 blosc blocks, so only those blocks are fetched (HTTP range request)
                 and decoded (blosc blocks are independent). Checked against full-chunk decode and
                 against WB2 where the two overlap (check_overlap).

usage: fields.py OUT_DIR [overlap]     (resumable: one npz per season in OUT_DIR)
Volume: about 15 GB streamed for 1979-2025 (WB2 ~7 GB, ARCO ~8 GB). Nothing but OUT_DIR (ignored) is written.
"""
import os, sys, json, struct, time, datetime as dt
import numpy as np, numcodecs, gcsfs
from concurrent.futures import ThreadPoolExecutor

WB2 = "weatherbench2/datasets/era5/1959-2023_01_10-6h-64x32_equiangular_conservative.zarr"
ARCO = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
WB2_LEVELS = (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000)
ARCO_LEVELS = (1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400,
               450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000)
G = 9.80665
LAT = -87.1875 + 5.625 * np.arange(32)
LON = 5.625 * np.arange(64)
NH = np.arange(20, 32)              # cell centres 25.3125 .. 87.1875 N
WB2_T0 = dt.datetime(1959, 1, 1)
WB2_LAST = dt.date(2023, 1, 9)
ARCO_T0 = dt.datetime(1900, 1, 1)
VARS = ("z500", "u250", "mslp", "sst", "tcwv")

_fs = None
_meta = None


def fs():
    global _fs
    if _fs is None:
        _fs = gcsfs.GCSFileSystem(token="anon")
    return _fs


def retry(f, *a, **k):
    for i in range(7):
        try:
            return f(*a, **k)
        except Exception as e:
            err = e
            time.sleep(2 ** i)
    raise err


def meta():
    global _meta
    if _meta is None:
        _meta = json.loads(retry(fs().cat, WB2 + "/.zmetadata"))["metadata"]
    return _meta


# ---------------------------------------------------------------- WB2
_wb2_cache = {}


def wb2_chunk(var, c):
    key = (var, c)
    m = meta()
    name = {"z500": "geopotential", "u250": "u_component_of_wind", "mslp": "mean_sea_level_pressure",
            "sst": "sea_surface_temperature", "tcwv": "total_column_water_vapour"}[var]
    za = m[name + "/.zarray"]
    codec = numcodecs.get_codec(za["compressor"])
    three = var in ("z500", "u250")
    raw = retry(fs().cat, f"{WB2}/{name}/{c}.0.0" + (".0" if three else ""))
    a = np.frombuffer(codec.decode(raw), za["dtype"]).reshape(za["chunks"])
    nt = min(za["chunks"][0], za["shape"][0] - c * za["chunks"][0])
    a = a[:nt]
    if three:
        a = a[:, WB2_LEVELS.index(500 if var == "z500" else 250)]
    a = a.transpose(0, 2, 1)[:, NH]                     # time, lat, lon
    if var == "z500":
        a = a / G
    if var == "mslp":
        a = a / 100.0
    return a.astype(np.float32), len(raw)


def wb2_days(dates, nth=8):
    """dict var -> [ndate, 12, 64] at 12 UTC (NaN where WB2 is NaN)."""
    idx = [((dt.datetime(d.year, d.month, d.day) - WB2_T0).days * 4 + 2) for d in dates]
    chunks = sorted({i // 100 for i in idx})
    out = {v: np.full((len(dates), len(NH), 64), np.nan, np.float32) for v in VARS}
    nb = [0]
    def job(vc):
        v, c = vc
        a, n = wb2_chunk(v, c)
        return v, c, a, n
    with ThreadPoolExecutor(nth) as ex:
        for v, c, a, n in ex.map(job, [(v, c) for v in VARS for c in chunks]):
            nb[0] += n
            for j, i in enumerate(idx):
                if i // 100 == c:
                    out[v][j] = a[i - c * 100]
    return out, nb[0]


# ---------------------------------------------------------------- regrid 0.25 -> 5.625 (conservative)
def _lat_weights():
    s = np.arange(721)
    a = np.clip(90 - 0.25 * s - 0.125, -90, 90)          # source cell lower edge
    b = np.clip(90 - 0.25 * s + 0.125, -90, 90)
    c = LAT
    e0 = np.maximum(c - 2.8125, -90)
    e1 = np.minimum(c + 2.8125, 90)
    W = np.zeros((32, 721))
    for k in range(32):
        lo = np.maximum(a, e0[k]); hi = np.minimum(b, e1[k])
        ov = np.maximum(0, np.sin(np.radians(hi)) - np.sin(np.radians(lo)))
        W[k] = ov
    return W / W.sum(1, keepdims=True)


def _lon_weights():
    s = np.arange(1440)
    a = 0.25 * s - 0.125
    b = 0.25 * s + 0.125
    W = np.zeros((64, 1440))
    for j in range(64):
        for sh in (-360, 0, 360):
            lo = np.maximum(a, 5.625 * j - 2.8125 + sh); hi = np.minimum(b, 5.625 * j + 2.8125 + sh)
            W[j] += np.maximum(0, hi - lo)
    return W / W.sum(1, keepdims=True)


WLAT = _lat_weights()[NH]
WLON = _lon_weights()


def regrid(a):
    """[721, 1440] (lat 90..-90) -> [12, 64]; NaN if any contributing source cell is NaN."""
    nan = np.isnan(a)
    z = WLAT @ np.where(nan, 0.0, a) @ WLON.T
    if nan.any():
        f = WLAT @ nan.astype(float) @ WLON.T
        z = np.where(f > 1e-6, np.nan, z)
    return z


# ---------------------------------------------------------------- ARCO
def arco_time(d):
    return int((dt.datetime(d.year, d.month, d.day, 12) - ARCO_T0).total_seconds() // 3600)


def arco_level(var, ti, level):
    """One level of a 37-level ARCO chunk via blosc blocks. Returns [721,1440] float32."""
    key = f"{ARCO}/{var}/{ti}.0.0.0"
    h = retry(fs().cat_file, key, start=0, end=16)
    ver, vl, fl, ts, nbytes, bs, cb = struct.unpack("<BBBBIII", h)
    nblk = -(-nbytes // bs)
    bst = np.frombuffer(retry(fs().cat_file, key, start=16, end=16 + 4 * nblk), "<u4")
    n = 721 * 1440
    L = ARCO_LEVELS.index(level)
    lo, hi = L * n * 4, (L + 1) * n * 4
    b0, b1 = lo // bs, (hi - 1) // bs
    ends = np.append(bst, cb)
    body = retry(fs().cat_file, key, start=int(bst[b0]), end=int(ends[b1 + 1]))
    nb = b1 - b0 + 1
    newb = (bst[b0:b1 + 1].astype("<u4") - bst[b0]) + 16 + 4 * nb
    last = (b1 == nblk - 1)
    nbytes_new = nb * bs if not last else (nb - 1) * bs + (nbytes - (nblk - 1) * bs)
    hdr = struct.pack("<BBBBIII", ver, vl, fl, ts, nbytes_new, bs, 16 + 4 * nb + len(body))
    buf = hdr + newb.astype("<u4").tobytes() + body
    arr = np.frombuffer(numcodecs.Blosc().decode(buf), "<f4")
    off = lo - b0 * bs
    return arr[off // 4: off // 4 + n].reshape(721, 1440), len(body) + 1200


def arco_full(var, ti):
    raw = retry(fs().cat, f"{ARCO}/{var}/{ti}.0.0")
    return np.frombuffer(numcodecs.Blosc().decode(raw), "<f4").reshape(721, 1440), len(raw)


def arco_day(d):
    ti = arco_time(d)
    nb = 0
    z, n = arco_level("geopotential", ti, 500); nb += n
    u, n = arco_level("u_component_of_wind", ti, 250); nb += n
    p, n = arco_full("mean_sea_level_pressure", ti); nb += n
    s, n = arco_full("sea_surface_temperature", ti); nb += n
    w, n = arco_full("total_column_water_vapour", ti); nb += n
    return (regrid(z / G), regrid(u), regrid(p / 100.0), regrid(s), regrid(w)), nb


def arco_days(dates, nth=8):
    out = {v: np.full((len(dates), len(NH), 64), np.nan, np.float32) for v in VARS}
    nb = 0
    with ThreadPoolExecutor(nth) as ex:
        for j, (r, n) in enumerate(ex.map(arco_day, dates)):
            nb += n
            for v, a in zip(VARS, r):
                out[v][j] = a
    return out, nb


# ---------------------------------------------------------------- driver
def season_dates(s):
    d0 = dt.date(s, 8, 20)
    return [d0 + dt.timedelta(days=i) for i in range((dt.date(s + 1, 5, 31) - d0).days + 1)]


def season(s):
    D = season_dates(s)
    w = [d for d in D if d <= WB2_LAST]
    a = [d for d in D if d > WB2_LAST]
    parts, nb = [], 0
    if w:
        o, n = wb2_days(w); parts.append(o); nb += n
    if a:
        o, n = arco_days(a); parts.append(o); nb += n
    out = {v: np.concatenate([p[v] for p in parts]) for v in VARS}
    return D, out, nb


def check_overlap(outdir):
    """ARCO regrid vs WB2 on 12 UTC of 12 days in Dec 2022 - Jan 2023 (all before 2023-01-10), and the block decode vs full decode."""
    ds = [dt.date(2022, 12, 1) + dt.timedelta(days=3 * i) for i in range(12)]
    w, _ = wb2_days(ds)
    a, nb = arco_days(ds)
    lines = []
    for v in VARS:
        d = a[v] - w[v]
        mw, ma = np.isnan(w[v]), np.isnan(a[v])
        lines.append(f"{v}: rms diff {np.sqrt(np.nanmean(d**2)):.4g}  wb2 std {np.nanstd(w[v]):.4g}  "
                     f"corr {np.corrcoef(a[v][~mw & ~ma], w[v][~mw & ~ma])[0, 1]:.5f}  nan-mask mismatch {(mw != ma).sum()} of {mw.size}")
    ti = arco_time(ds[0])
    full = np.frombuffer(numcodecs.Blosc().decode(retry(fs().cat, f"{ARCO}/geopotential/{ti}.0.0.0")), "<f4").reshape(37, 721, 1440)
    blk, n = arco_level("geopotential", ti, 500)
    lines.append(f"block decode == full decode (z500): {np.array_equal(full[ARCO_LEVELS.index(500)], blk)}; "
                 f"bytes {n/1e6:.2f} MB of 80 MB")
    txt = "\n".join(lines) + f"\nARCO bytes for 12 days: {nb/1e6:.1f} MB\n"
    open(os.path.join(outdir, "overlap_check.txt"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    if len(sys.argv) > 2 and sys.argv[2] == "overlap":
        check_overlap(out); sys.exit()
    tot = 0
    for s in range(2004, 2026):
        fn = os.path.join(out, f"season_{s}.npz")
        if os.path.exists(fn):
            continue
        t = time.time()
        D, o, nb = season(s)
        tot += nb
        np.savez_compressed(fn, dates=np.array([d.isoformat() for d in D]), **o)
        print(s, f"{nb/1e9:.2f} GB  {time.time()-t:.0f}s  running {tot/1e9:.1f} GB", flush=True)
