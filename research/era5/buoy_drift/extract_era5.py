"""ERA5 10 m wind speed and instantaneous 10 m gust at the moored-buoy sites.

ERA5 here is a PROXY record (see STATUS.md); this script only samples it.

Source: Google's public ARCO-ERA5 Zarr stores (anonymous HTTPS):
  10 m wind speed      ar/1959-2022-6h-1440x721.zarr/10m_wind_speed        (6-hourly)
  instantaneous gust   ar/full_37-1h-0p25deg-chunk-1.zarr-v3/instantaneous_10m_wind_gust
  land-sea mask        ar/full_37-1h-0p25deg-chunk-1.zarr-v3/land_sea_mask (2000-01-01 00Z)
These are the same gust variable pipelines A and B use. The wind speed is the
instantaneous analysis speed sqrt(u10^2 + v10^2) at the synoptic hour.

Sampling: 00, 06, 12 and 18 UTC, October through March, 1979-01-01 to
2004-12-31 (the span of the ISD buoy record). About 18,900 times.

Volume: each Zarr chunk is one global hour (721 x 1440 float32, blosc-lz4,
about 3.2 MB). Blosc compresses in independent 512 KiB blocks of 91 latitude
rows, so only the two blocks that hold 22.0-67.25 N are fetched (HTTP range
requests), with the 48-byte header. That is about 0.8 MB per variable per
time, so about 30 GB in all, against 120 GB for whole chunks. Nothing is
stored but a 3 x 3 grid-point neighbourhood around each buoy.

    python3 research/era5/buoy_drift/extract_era5.py [--procs 24]
Writes work/era5/<yyyymm>.npz (resumable: existing months are skipped).
"""
import argparse
import concurrent.futures as cf
import datetime as dt
import os
import struct
import time
import urllib.request

import lz4.block
import numpy as np

from fetch_buoys import STATIONS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "work", "era5")
G = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/"
WS = G + "1959-2022-6h-1440x721.zarr/10m_wind_speed/{i}.0.0"          # i = hours since 1959 / 6
GU = G + "full_37-1h-0p25deg-chunk-1.zarr-v3/instantaneous_10m_wind_gust/{i}.0.0"  # i = hours since 1900
LSM = G + "full_37-1h-0p25deg-chunk-1.zarr-v3/land_sea_mask/876576.0.0"
NLON = 1440
BLOCKS = (1, 2)                     # rows 91-181 (67.25-45.0 N) and 182-272 (44.75-22.0 N)
ROWS_PER_BLOCK = 91                 # 524288 bytes / 4 / 1440, not an integer; handled by flat index


def get(url, rng=None, tries=6):
    for k in range(tries):
        try:
            req = urllib.request.Request(url)
            if rng:
                req.add_header("Range", "bytes=%d-%d" % rng)
            return urllib.request.urlopen(req, timeout=120).read()
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 ** k)


def decode_block(buf, flags, ts, bsize):
    """Decompress one blosc-1 block (lz4, byte shuffle, split streams)."""
    split = not (flags & 0x10)
    ns = ts if split else 1
    neb = bsize // ns
    out = bytearray()
    p = 0
    for _ in range(ns):
        cs = struct.unpack("<i", buf[p:p + 4])[0]
        p += 4
        out += buf[p:p + cs] if cs == neb else lz4.block.decompress(bytes(buf[p:p + cs]), uncompressed_size=neb)
        p += cs
    a = np.frombuffer(bytes(out), np.uint8)
    if flags & 0x1:
        a = a.reshape(ts, -1).T.reshape(-1)
    return a.view("<f4")


def read_blocks(url, flat_idx):
    """Values at flat indices (row*1440+col) of one chunk, fetching only BLOCKS."""
    h = get(url, (0, 47))
    flags, ts = h[2], h[3]
    nbytes, bs, cbytes = struct.unpack("<iii", h[4:16])
    if flags & 0x2 or bs != 524288 or ts != 4:
        raise RuntimeError("unexpected blosc layout in " + url)
    st = np.frombuffer(h[16:48], "<i4")
    b0, b1 = BLOCKS
    raw = get(url, (int(st[b0]), int(st[b1 + 1]) - 1))
    vals = {}
    for b in BLOCKS:
        blk = decode_block(raw[st[b] - st[b0]: st[b + 1] - st[b0]], flags, ts, bs)
        lo = b * bs // 4
        for j in flat_idx:
            if lo <= j < lo + len(blk):
                vals[j] = blk[j - lo]
    return np.array([vals[j] for j in flat_idx], np.float32)


def neighbourhoods():
    """Flat indices of the 3 x 3 grid points around each buoy (nominal ISD position)."""
    import csv, gzip
    idx = []
    for sid in STATIONS:
        with gzip.open(os.path.join(HERE, "work", "isd", sid + ".csv.gz"), "rt") as f:
            r = next(csv.DictReader(f))
        lat, lon = float(r["LATITUDE"]), float(r["LONGITUDE"]) % 360
        r0, c0 = int(round((90 - lat) * 4)), int(round(lon * 4))
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                idx.append((r0 + dr) * NLON + (c0 + dc) % NLON)
    return np.array(idx)


def times_for_month(y, m):
    t = dt.datetime(y, m, 1)
    out = []
    while t.month == m:
        out.append(t)
        t += dt.timedelta(hours=6)
    return out


def do_month(args):
    y, m, idx = args
    path = os.path.join(OUT, "%04d%02d.npz" % (y, m))
    if os.path.exists(path):
        return path, 0
    ts = times_for_month(y, m)
    ws = np.empty((len(ts), len(idx)), np.float32)
    gu = np.empty_like(ws)
    for k, t in enumerate(ts):
        h1900 = int((t - dt.datetime(1900, 1, 1)).total_seconds() // 3600)
        h1959 = int((t - dt.datetime(1959, 1, 1)).total_seconds() // 3600)
        ws[k] = read_blocks(WS.format(i=h1959 // 6), idx)
        gu[k] = read_blocks(GU.format(i=h1900), idx)
    np.savez_compressed(path + ".tmp.npz", time=np.array(ts, "datetime64[h]"), ws=ws, gust=gu, idx=idx)
    os.replace(path + ".tmp.npz", path)
    return path, len(ts)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=24)
    ap.add_argument("--test", action="store_true", help="one month only")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    idx = neighbourhoods()
    lsm_path = os.path.join(OUT, "lsm.npz")
    if not os.path.exists(lsm_path):
        from numcodecs import Blosc
        lsm = np.frombuffer(Blosc().decode(get(LSM)), "<f4")
        np.savez_compressed(lsm_path, lsm=lsm[idx], idx=idx)
    months = [(y, m) for y in range(2004, 1978, -1) for m in (12, 11, 10, 3, 2, 1)]
    if a.test:
        months = months[:1]
    t0 = time.time()
    done = 0
    with cf.ThreadPoolExecutor(a.procs) as ex:
        for path, n in ex.map(do_month, [(y, m, idx) for y, m in months]):
            done += 1
            print("%s %d  [%d/%d, %.0f s]" % (os.path.basename(path), n, done, len(months), time.time() - t0), flush=True)
