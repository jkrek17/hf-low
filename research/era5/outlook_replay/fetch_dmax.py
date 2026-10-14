"""Pull GFS GUST:surface (00 UTC cycles, 1 Oct-30 Apr, 2021-22..2025-26, valid every 6 h f006..f096) by .idx byte ranges
from the NOAA open-data bucket and save, per cycle, the forecast-DAY maximum gust on a 0.5 degree ocean grid.
Day d (1..4) = valid hours 24(d-1)+6 .. 24d (four 6-hourly samples). Grid: 110 rows x 720 cols; row k = lat 74.75-0.5k, col c = lon 0.125+0.5c
(2x2 maximum of the 0.25 degree field, land cells set to 0 first). Stored uint8 in 0.5 kt steps, 255 = a sample of that day was missing.
GFS GUST is in m/s in the file; converted to kt. Decoder: grib2_t53.py. Resumable: a cycle with an npz is skipped.
usage: fetch_dmax.py WORKDIR [NPROC]"""
import os, sys, time, struct, urllib.request, urllib.error
import numpy as np
from datetime import date, timedelta
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import grib2_t53 as G
BUCKET = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"; UA = "hf-low-outlook-replay/0.1"
LEADS = list(range(6, 97, 6)); KT = 1.9438445
WORK = sys.argv[1]; NPROC = int(sys.argv[2]) if len(sys.argv) > 2 else 12
LSM = np.load("/mnt/project-files/era5-cache/lsm_nh_0p25.npy")[60:280] < 0.5     # 220 rows 75N..20.25N, ocean True
LSM = np.concatenate([LSM, LSM[:, :0]], axis=1)


class Missing(Exception):
    pass


def get(u, rng=None):
    h = {"User-Agent": UA}; want = None
    if rng:
        h["Range"] = rng; a, b = rng.split("=")[1].split("-"); want = (int(b) - int(a) + 1) if b else None
    err = None; n404 = 0
    for k in range(8):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=120) as r:
                b = r.read(); cl = r.headers.get("Content-Length")
            if (cl is not None and len(b) != int(cl)) or (want is not None and len(b) != want):
                raise IOError("short read")
            return b
        except urllib.error.HTTPError as e:
            if e.code == 404:
                n404 += 1
                if n404 >= 3: raise Missing(u)
            err = e
        except Exception as e:
            err = e
        time.sleep(min(2 ** k, 30))
    raise RuntimeError(f"fetch failed {u}: {err}")


def walk(u, i):
    off = 0
    for k in range(i + 1):
        h = get(u, f"bytes={off}-{off + 15}")
        if h[:4] != b"GRIB": raise RuntimeError("walk lost sync")
        L = struct.unpack(">Q", h[8:16])[0]
        if k == i: return get(u, f"bytes={off}-{off + L - 1}")
        off += L


def gust_field(day, lead):
    u = f"{BUCKET}/gfs.{day}/00/atmos/gfs.t00z.pgrb2.0p25.f{lead:03d}"
    idx = get(u + ".idx"); rows = []
    for line in idx.decode().splitlines():
        p = line.split(":")
        if len(p) >= 6 and p[1].isdigit(): rows.append((int(p[1]), p[3], p[4]))
    for i, (off, v, l) in enumerate(rows):
        if v == "GUST" and l == "surface":
            end = rows[i + 1][0] - 1
            blob = get(u, f"bytes={off}-{end}")
            try:
                return G.decode_message(blob)[60:280] * KT * 1.0, len(blob) + len(idx)
            except (AssertionError, struct.error, IndexError, KeyError):
                blob = walk(u, i)       # stale .idx (offsets a few hundred bytes off, seen Nov 2022): walk GRIB section-0 lengths
                return G.decode_message(blob)[60:280] * KT * 1.0, len(blob) + len(idx)
    raise Missing(u + " no GUST")


def job(day):
    out = os.path.join(WORK, "cycles", f"{day}.npz")
    if os.path.exists(out): return day, 0, "skip"
    samples = {}; nb = 0
    for lead in LEADS:
        try:
            g, n = gust_field(day, lead); nb += n
        except Missing:
            continue
        g = np.where(LSM, g, 0.0)
        samples[lead] = g.reshape(110, 2, 720, 2).max(axis=(1, 3))
    if not samples: return day, nb, "missing"
    d = np.full((4, 110, 720), 255, np.uint8)
    for k in range(4):
        ls = [24 * k + 6 * j for j in (1, 2, 3, 4)]
        if all(l in samples for l in ls):
            m = np.max([samples[l] for l in ls], axis=0)
            d[k] = np.clip(np.round(m * 2), 0, 254).astype(np.uint8)
    np.savez_compressed(out + ".tmp.npz", d=d, nmiss=np.array([l for l in LEADS if l not in samples]))
    os.replace(out + ".tmp.npz", out)
    return day, nb, "ok"


if __name__ == "__main__":
    os.makedirs(os.path.join(WORK, "cycles"), exist_ok=True)
    days = []
    for y in range(2021, 2026):
        d = date(y, 10, 1)
        while d <= date(y + 1, 4, 30):
            days.append(d.strftime("%Y%m%d")); d += timedelta(days=1)
        # October-April only: the loop above spans Oct y .. Apr y+1 continuously (it includes May-Sep gaps)
    days = [x for x in days if int(x[4:6]) in (10, 11, 12, 1, 2, 3, 4)]
    if os.environ.get("ONLY"): days = os.environ["ONLY"].split(",")
    t0 = time.time(); tot = 0
    with Pool(NPROC) as p, open(os.path.join(WORK, "bytes.log"), "a") as lg:
        for k, (day, nb, st) in enumerate(p.imap_unordered(job, days)):
            tot += nb; lg.write(f"{day} {st} {nb}\n"); lg.flush()
            if k % 25 == 0: print(k, "of", len(days), f"{tot / 1e9:.2f} GB {time.time() - t0:.0f}s", flush=True)
    print("done", len(days), f"{tot / 1e9:.2f} GB")
