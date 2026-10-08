"""Gridded ERA5 10 m wind-gust climatology, October-April, 12-hourly (00 and 12 UTC), seasons 2004-05 to 2025-26.

ERA5 PROXY (reanalysis gust), not observation. Streams ARCO-ERA5 instantaneous_10m_wind_gust, one global field per time
(3.2 MB each, measured on 3 test times): 9,328 times, about 30 GB. Under the 50 GB gate (CLAUDE.md); no existing cache covers it
(hf_structure pulled only HF-strength times). Fields are aggregated on the fly; nothing is stored per time.

Per grid cell it accumulates, over the sampled times: n, sum, max, counts of gust >= THRESH kt, and the maximum within each
season (for the mean annual maximum). Sampling is 12-hourly, so maxima and exceedance frequencies are for the sampled
times, not the true maxima.

usage: ERA5_WORK=<dir> python3 extract.py [NPROC]       (resumable: finished chunks are skipped)
"""
import os, sys, time, urllib.request, datetime
import numpy as np, numcodecs
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")))
ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
T0 = datetime.datetime(1900, 1, 1)
KT = 1 / 0.514444
THRESH = (34.0, 50.0, 64.0, 71.7, 80.0, 90.0)
R0, R1 = 0, 361                       # rows 90N..0N, whole northern hemisphere
SEASONS = list(range(2004, 2026))
_codec = numcodecs.Blosc()


def times():
    out = []
    for s in SEASONS:
        t = datetime.datetime(s, 10, 1)
        while t < datetime.datetime(s + 1, 5, 1):
            for h in (0, 12):
                out.append((s, t + datetime.timedelta(hours=h)))
            t += datetime.timedelta(days=1)
    return out


def fetch(url):
    for k in range(6):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            return np.frombuffer(_codec.decode(raw), "<f4").reshape(721, 1440), len(raw)
        except Exception:
            if k == 5:
                raise
            time.sleep(2 ** k)


def chunk(args):
    ci, items = args
    path = os.path.join(WORK, f"chunk_{ci:04d}.npz")
    if os.path.exists(path):
        return ci, 0
    shp = (R1 - R0, 1440)
    n = 0
    sm = np.zeros(shp)
    mx = np.zeros(shp, np.float32)
    ex = np.zeros((len(THRESH),) + shp, np.uint32)
    smax = {}
    nbytes = 0
    for s, t in items:
        ti = int((t - T0).total_seconds() // 3600)
        g, b = fetch(f"{ARCO}/instantaneous_10m_wind_gust/{ti}.0.0")
        nbytes += b
        g = g[R0:R1] * KT
        n += 1
        sm += g
        mx = np.maximum(mx, g)
        for k, th in enumerate(THRESH):
            ex[k] += g >= th
        smax[s] = np.maximum(smax.get(s, 0), g)
    np.savez(path + ".tmp.npz", n=n, sm=sm, mx=mx, ex=ex, seasons=np.array(sorted(smax)),
             smax=np.stack([smax[s] for s in sorted(smax)]).astype(np.float16), nbytes=nbytes)
    os.replace(path + ".tmp.npz", path)
    return ci, nbytes


if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True)
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    T = times()
    # chunk by season-contiguous blocks of ~100 times so each chunk holds at most 2 seasons
    chunks = [(i, T[i * 100:(i + 1) * 100]) for i in range((len(T) + 99) // 100)]
    print(f"{len(T)} times, {len(chunks)} chunks, about {len(T) * 3.2 / 1000:.1f} GB", flush=True)
    t0 = time.time()
    with Pool(nproc) as p:
        for ci, nb in p.imap_unordered(chunk, chunks):
            print(f"chunk {ci} done {nb / 1e9:.2f} GB  t={time.time() - t0:.0f}s", flush=True)
