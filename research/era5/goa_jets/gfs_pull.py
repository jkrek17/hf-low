"""Follow-up stage 2 (PREREGISTRATION_followup.md): GFS 0.25 degree forecasts at 24 and 48 h lead valid at every 00/12 UTC, October-April, seasons 2021-22..2025-26.

Records PRMSL (mean sea level), UGRD and VGRD (10 m), GUST (surface) by .idx byte ranges from the NOAA open-data bucket; cropped to rows 80:185, cols 740:941
(70N..43.75N, 175W..125W, the same window as the ERA5 study). One npz per (lead, valid time) in $ERA5_WORK/goa_jets/gfs/ with msl (hPa), u, v, gust (kt), float32.
GFS FORECAST (model), not observation. Resumable; writes bytes.log. usage: gfs_pull.py probe | gfs_pull.py run [NPROC]
"""
import os, sys, time, struct, urllib.request, urllib.error
import numpy as np, pandas as pd
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "outlook_replay"))
import grib2_t53 as G
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "goa_jets", "gfs")
BUCKET = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"; UA = "hf-low-goa-jets/0.1"
KT = 1.9438445; LEADS = (24, 48)
WANT = {("PRMSL", "mean sea level"): "msl", ("UGRD", "10 m above ground"): "u", ("VGRD", "10 m above ground"): "v", ("GUST", "surface"): "gust"}
ROW, COL = slice(80, 185), slice(740, 941)


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


def fields(valid, lead):
    cyc = valid - pd.Timedelta(hours=lead)
    u = f"{BUCKET}/gfs.{cyc:%Y%m%d}/{cyc:%H}/atmos/gfs.t{cyc:%H}z.pgrb2.0p25.f{lead:03d}"
    idx = get(u + ".idx"); rows = []
    for line in idx.decode().splitlines():
        p = line.split(":")
        if len(p) >= 6 and p[1].isdigit(): rows.append((int(p[1]), p[3], p[4]))
    out = {}; nb = len(idx)
    for i, (off, v, l) in enumerate(rows):
        if (v, l) in WANT and WANT[(v, l)] not in out:
            end = rows[i + 1][0] - 1
            blob = get(u, f"bytes={off}-{end}"); nb += len(blob)
            try:
                a = G.decode_message(blob)
            except (AssertionError, struct.error, IndexError, KeyError):
                blob = walk(u, i); a = G.decode_message(blob)
            assert a.shape == (721, 1440)
            out[WANT[(v, l)]] = a[ROW, COL]
    if len(out) != 4: raise Missing(u + " records " + str(sorted(out)))
    return dict(msl=(out["msl"] / 100).astype(np.float32), u=(out["u"] * KT).astype(np.float32), v=(out["v"] * KT).astype(np.float32),
                gust=(out["gust"] * KT).astype(np.float32)), nb


def job(a):
    valid, lead = a
    fn = os.path.join(WORK, f"{lead}_{valid:%Y%m%d%H}.npz")
    if os.path.exists(fn): return valid, lead, 0, "skip"
    try:
        d, nb = fields(valid, lead)
    except Missing:
        return valid, lead, 0, "missing"
    np.savez_compressed(fn + ".tmp.npz", **d); os.replace(fn + ".tmp.npz", fn)
    return valid, lead, nb, "ok"


def vtimes():
    d = pd.read_csv(os.path.join(HERE, "results", "stage0_zone_gust.csv"), parse_dates=["time"])
    t = d.time[(d.time >= "2021-10-01") & (d.time < "2026-05-01")]
    return list(t)


if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True)
    T = vtimes(); jobs = [(t, l) for t in T for l in LEADS]
    print(len(T), "valid times,", len(jobs), "files; first", T[0], "last", T[-1])
    if sys.argv[1] == "probe":
        t0 = time.time(); r = job(jobs[0]); print(r, f"{time.time() - t0:.1f}s"); sys.exit()
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    tot = 0; t0 = time.time(); miss = 0
    with Pool(n) as p, open(os.path.join(WORK, "bytes.log"), "a") as lg:
        for k, (v, l, nb, st) in enumerate(p.imap_unordered(job, jobs)):
            tot += nb; miss += st == "missing"; lg.write(f"{v:%Y%m%d%H} {l} {st} {nb}\n"); lg.flush()
            if k % 100 == 0: print(k, "of", len(jobs), f"{tot / 1e9:.2f} GB missing {miss} {time.time() - t0:.0f}s", flush=True)
    print("done", len(jobs), f"{tot / 1e9:.2f} GB missing {miss}")
