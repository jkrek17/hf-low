"""Compact all-month ERA5 10 m gust cache, 12-hourly (00 and 12 UTC), 1 June 2004 to 31 May 2026 (seasons 2004-05 to 2025-26).

ERA5 PROXY (reanalysis instantaneous gust), not observation. Streams ARCO-ERA5 instantaneous_10m_wind_gust (one global field per
time, ~3.2 MB): 16,070 times, ~51 GB streamed. Keeps the northern hemisphere (rows 90N..0N, 361 x 1440) as uint8 in 0.5 kt steps
(value = round(gust_kt * 2), clipped at 255 = 127.5 kt), zstd-compressed, 100 times per file. Nothing else is kept. Later gust
questions read this cache instead of pulling again.

usage: python3 -I extract_cache.py OUT_DIR [NPROC]       (resumable: finished chunks are skipped)
"""
import os, sys, time, urllib.request, datetime
import numpy as np, numcodecs
from multiprocessing import Pool

ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
T0 = datetime.datetime(1900, 1, 1)
KT = 1 / 0.514444
R1 = 361
SEASONS = list(range(2004, 2026))
_in = numcodecs.Blosc()
_out = numcodecs.Blosc(cname="zstd", clevel=5, shuffle=numcodecs.Blosc.BITSHUFFLE)


def times():
    out = []
    t = datetime.datetime(2004, 6, 1)
    while t < datetime.datetime(2026, 6, 1):
        for h in (0, 12):
            out.append(t + datetime.timedelta(hours=h))
        t += datetime.timedelta(days=1)
    return out


def fetch(url):
    for k in range(6):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            return np.frombuffer(_in.decode(raw), "<f4").reshape(721, 1440), len(raw)
        except Exception:
            if k == 5:
                raise
            time.sleep(2 ** k)


def chunk(args):
    ci, items, out = args
    path = os.path.join(out, f"gust_{ci:04d}.npz")
    if os.path.exists(path):
        return ci, 0
    arr = np.zeros((len(items), R1, 1440), np.uint8)
    nbytes = 0
    for i, t in enumerate(items):
        ti = int((t - T0).total_seconds() // 3600)
        g, b = fetch(f"{ARCO}/instantaneous_10m_wind_gust/{ti}.0.0")
        nbytes += b
        arr[i] = np.clip(np.rint(g[:R1] * KT * 2), 0, 255).astype(np.uint8)
    tv = np.array([int((t - T0).total_seconds() // 3600) for t in items], np.int64)
    comp = np.frombuffer(_out.encode(arr), np.uint8)
    np.savez(path + ".tmp.npz", hours_since_1900=tv, shape=np.array(arr.shape), data=comp)
    os.replace(path + ".tmp.npz", path)
    return ci, nbytes


if __name__ == "__main__":
    out = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    os.makedirs(out, exist_ok=True)
    T = times()
    chunks = [(i, T[i * 100:(i + 1) * 100], out) for i in range((len(T) + 99) // 100)]
    print(f"{len(T)} times, {len(chunks)} chunks, about {len(T) * 3.2 / 1000:.1f} GB", flush=True)
    t0 = time.time()
    with Pool(nproc) as p:
        for ci, nb in p.imap_unordered(chunk, chunks):
            print(f"chunk {ci} done {nb / 1e9:.2f} GB  t={time.time() - t0:.0f}s", flush=True)
