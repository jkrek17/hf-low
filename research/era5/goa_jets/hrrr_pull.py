"""HRRR-Alaska pull (PREREGISTRATION_hrrr.md): GUST surface, UGRD/VGRD 10 m, MSLMA from the 3 hour forecast of the cycle 3 h before each 00/12 UTC time, Oct-Apr 2018-19..2025-26.

By .idx byte ranges from noaa-hrrr-bdp-pds; winds rotated to east/north; cropped to the box 48-67N, 160-128W. One npz per valid time in $ERA5_WORK/goa_jets/hrrr/:
gust (kt, float16), u, v (kt, float16), msl (0.1 hPa, uint16). land.npy, crop.npz (lat, lon) written once. MODEL values, not observation. Resumable; bytes.log.
usage: hrrr_pull.py probe | hrrr_pull.py run [NPROC]
"""
import os, sys, time, struct, urllib.request, urllib.error
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import grib2_hrrr as G, hrrr_grid as HG
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "goa_jets", "hrrr")
BUCKET = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"; UA = "hf-low-goa-jets-hrrr/0.1"
KT = 1.9438445
WANT = {("GUST", "surface"): "gust", ("UGRD", "10 m above ground"): "ug", ("VGRD", "10 m above ground"): "vg", ("MSLMA", "mean sea level"): "msl"}
LA, LO, ALPHA = HG.grid()
_m = (LA >= 48) & (LA <= 67) & (LO >= -160) & (LO <= -128)
J0, J1 = np.where(_m.any(axis=1))[0][[0, -1]]; I0, I1 = np.where(_m.any(axis=0))[0][[0, -1]]
SJ, SI = slice(J0, J1 + 1), slice(I0, I1 + 1)


class Missing(Exception):
    pass


def get(u, rng=None):
    h = {"User-Agent": UA}; want = None
    if rng:
        h["Range"] = rng; a, b = rng.split("=")[1].split("-"); want = int(b) - int(a) + 1
    err = None; n404 = 0
    for k in range(8):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=120) as r:
                b = r.read(); cl = r.headers.get("Content-Length")
            if (cl is not None and len(b) != int(cl)) or (want is not None and len(b) != want):
                raise IOError("short read")
            return b
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
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


def url(valid):
    cyc = valid - pd.Timedelta(hours=3)
    return f"{BUCKET}/hrrr.{cyc:%Y%m%d}/alaska/hrrr.t{cyc:%H}z.wrfsfcf03.ak.grib2"


def fields(valid, extra=()):
    u = url(valid); idx = get(u + ".idx"); rows = []
    for line in idx.decode().splitlines():
        p = line.split(":")
        if len(p) >= 6 and p[1].isdigit(): rows.append((int(p[1]), p[3], p[4], p[5]))
    want = dict(WANT); want.update({e: e[0].lower() for e in extra})
    out = {}; nb = len(idx)
    for i, (off, v, l, f) in enumerate(rows):
        if (v, l) in want and want[(v, l)] not in out:
            assert "3 hour fcst" in f, f
            end = rows[i + 1][0] - 1 if i + 1 < len(rows) else ""
            blob = get(u, f"bytes={off}-{end}"); nb += len(blob)
            try:
                a = G.decode_message(blob)
            except (AssertionError, struct.error, IndexError, KeyError) as e:
                if "packing" in str(e) or "grid" in str(e): raise
                blob = walk(u, i); a = G.decode_message(blob)
            assert a.shape == (919, 1299)
            out[want[(v, l)]] = a
    if not all(k in out for k in ("gust", "ug", "vg", "msl")): raise Missing(u + " records " + str(sorted(out)))
    return out, nb


def job(valid):
    fn = os.path.join(WORK, f"{valid:%Y%m%d%H}.npz")
    if os.path.exists(fn): return valid, 0, "skip"
    try:
        o, nb = fields(valid)
    except Missing:
        return valid, 0, "missing"
    ue, ve = HG.to_earth(o["ug"], o["vg"], ALPHA)
    np.savez_compressed(fn + ".tmp.npz", gust=(o["gust"] * KT)[SJ, SI].astype(np.float16), u=(ue * KT)[SJ, SI].astype(np.float16), v=(ve * KT)[SJ, SI].astype(np.float16),
                        msl=np.round(o["msl"][SJ, SI] / 10).astype(np.uint16))
    os.replace(fn + ".tmp.npz", fn)
    return valid, nb, "ok"


def setup():
    os.makedirs(WORK, exist_ok=True)
    np.savez_compressed(os.path.join(WORK, "crop.npz"), lat=LA[SJ, SI].astype(np.float32), lon=LO[SJ, SI].astype(np.float32), j0=J0, i0=I0)
    if not os.path.exists(os.path.join(WORK, "land.npy")):
        o, _ = fields(pd.Timestamp("2023-01-16 00:00"), extra=(("LAND", "surface"),))
        np.save(os.path.join(WORK, "land.npy"), o["land"][SJ, SI].astype(np.float16))


def vtimes():
    d = pd.read_csv(os.path.join(HERE, "results", "stage0_zone_gust.csv"), parse_dates=["time"])
    return list(d.time[(d.time >= "2018-10-01") & (d.time < "2026-05-01")])


if __name__ == "__main__":
    setup(); T = vtimes()
    print(len(T), "valid times; crop", SJ, SI, "->", (J1 - J0 + 1, I1 - I0 + 1))
    if sys.argv[1] == "probe":
        t0 = time.time(); print(job(T[0]), f"{time.time() - t0:.1f}s"); sys.exit()
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    tot = 0; t0 = time.time(); miss = 0
    with Pool(n) as p, open(os.path.join(WORK, "bytes.log"), "a") as lg:
        for k, (v, nb, st) in enumerate(p.imap_unordered(job, T)):
            tot += nb; miss += st == "missing"; lg.write(f"{v:%Y%m%d%H} {st} {nb}\n"); lg.flush()
            if k % 100 == 0: print(k, "of", len(T), f"{tot / 1e9:.2f} GB missing {miss} {time.time() - t0:.0f}s", flush=True)
    print("done", len(T), f"{tot / 1e9:.2f} GB missing {miss}")
