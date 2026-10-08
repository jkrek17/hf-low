"""Stage A GFS pull: PRMSL + GUST:surface by .idx byte ranges (method of wind-particles/server.py), processed on the fly.

Per cycle (00/12 UTC, f000..f048 step 6, 1 Oct .. 30 Apr, seasons 2021..2025) writes, in GFS_WORK (/home/claude/gfs_work):
  lows/<YYYYMMDDCC>.csv   cycle,lead,lat,lon,msl,g500,g800,g1200      (a_tracks.lows_in, kt)
  gust/<YYYYMMDDCC>.csv   cycle,lead,basin,gmax_kt,g999_kt           (ocean points in the basin box)
  missing/<YYYYMMDDCC>.txt  leads whose files are absent (404 on the .idx or the message)
Grids are never stored. Resumable: a cycle with a lows file is skipped. bytes.log: one line per file, bytes downloaded.
usage: a_fetch.py [--day YYYYMMDD] [--workers N]      (no --day: all seasons)
"""
import os, sys, time, argparse, urllib.request, urllib.error, tempfile
import numpy as np
from datetime import date, timedelta
from multiprocessing import Pool
from eccodes import codes_get, codes_get_values, codes_new_from_message, codes_release

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a_tracks as A

WORK = A.WORK
BUCKET = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
UA = "hf-low-gfs-test/0.1"
LEADS = list(range(0, 49, 6))
KT = A.KT


def url(day, cc, lead):
    return f"{BUCKET}/gfs.{day}/{cc}/atmos/gfs.t{cc}z.pgrb2.0p25.f{lead:03d}"


class Missing(Exception):
    pass


def get(u, rng=None):
    """Returns bytes; raises Missing only on a repeated HTTP 404 from the bucket. Everything else (403/5xx from the
    proxy or S3 under load, short reads, resets) is retried with backoff and, if it persists, raises (cycle not written)."""
    h = {"User-Agent": UA}
    want = None
    if rng:
        h["Range"] = rng
        a, b = rng.split("=")[1].split("-")
        want = (int(b) - int(a) + 1) if b else None
    err = None
    n404 = 0
    for k in range(10):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=120) as r:
                b = r.read()
                cl = r.headers.get("Content-Length")
            if (cl is not None and len(b) != int(cl)) or (want is not None and len(b) != want):
                raise IOError(f"short read {len(b)} of {cl}/{want}")
            return b
        except urllib.error.HTTPError as e:
            if e.code == 404:
                n404 += 1
                if n404 >= 6:
                    raise Missing(u)
            err = e
        except Exception as e:
            err = e
        time.sleep(min(2 ** (k + (1 if n404 else 0)), 30))
    raise RuntimeError(f"fetch failed {u}: {err}")


def idx_rows(text):
    rows = []
    for line in text.splitlines():
        p = line.strip().split(":")
        if len(p) >= 6 and p[1].isdigit():
            rows.append((int(p[1]), p[3], p[4]))
    return rows


def message_walk(u, rows, var, level):
    """Fallback for a stale .idx (seen late Nov 2022: offsets a few hundred bytes off): walk the GRIB section-0 length
    fields from byte 0 to the message at the idx position of var:level and read that message exactly."""
    import struct
    i = next(k for k, (o, v, l) in enumerate(rows) if v == var and l == level)
    off = 0
    for k in range(i + 1):
        h = get(u, f"bytes={off}-{off + 15}")
        if h[:4] != b"GRIB":
            raise RuntimeError("walk lost sync")
        L = struct.unpack(">Q", h[8:16])[0]
        if k == i:
            return get(u, f"bytes={off}-{off + L - 1}")
        off += L


def message(u, rows, var, level):
    for i, (off, v, l) in enumerate(rows):
        if v == var and l == level:
            end = rows[i + 1][0] - 1 if i + 1 < len(rows) else None
            rng = f"bytes={off}-{end}" if end is not None else f"bytes={off}-"
            return get(u, rng)
    raise Missing(f"{var}:{level} not in idx of {u}")


def decode(blob, short):
    g = codes_new_from_message(blob)
    try:
        assert codes_get(g, "shortName") == short, codes_get(g, "shortName")
        ni, nj = int(codes_get(g, "Ni")), int(codes_get(g, "Nj"))
        assert (ni, nj) == (1440, 721), (ni, nj)
        v = np.array(codes_get_values(g), dtype=np.float32).reshape(nj, ni)
        assert not int(codes_get(g, "iScansNegatively")) and not int(codes_get(g, "jScansPositively"))
        assert abs(float(codes_get(g, "latitudeOfFirstGridPointInDegrees")) - 90) < 1e-6
        assert abs(float(codes_get(g, "longitudeOfFirstGridPointInDegrees"))) < 1e-6
        return v
    finally:
        codes_release(g)


def land_mask(day, cc):
    p = os.path.join(WORK, "land.npy")
    if os.path.exists(p):
        return np.load(p) == 0
    u = url(day, cc, 0)
    b = get(u + ".idx")
    rows = idx_rows(b.decode())
    blob = message(u, rows, "LAND", "surface")
    g = codes_new_from_message(blob)
    try:
        v = np.array(codes_get_values(g), dtype=np.float32).reshape(721, 1440)
    finally:
        codes_release(g)
    np.save(p + ".tmp.npy", v)
    os.replace(p + ".tmp.npy", p)
    with open(os.path.join(WORK, "bytes.log"), "a") as f:
        f.write(f"LAND {day}{cc} f000 {len(b) + len(blob)}\n")
    return v == 0


def box_masks():
    la = (90 - 0.25 * np.arange(721))[:, None]
    lo = (0.25 * np.arange(1440))[None, :]
    return {"atl": (la >= 30) & (la <= 67) & ((lo >= 262) | (lo <= 10)),
            "pac": (la >= 27) & (la <= 67) & (lo >= 135) & (lo <= 240)}


def do_cycle(args):
    try:
        return _cycle(args)
    except Exception as e:
        return f"{args[0]}{args[1]}", f"ERROR {e!r}"


def _cycle(args):
    day, cc = args
    cyc = f"{day}{cc}"
    lo_p = os.path.join(WORK, "lows", cyc + ".csv")
    if os.path.exists(lo_p):
        return cyc, "skip"
    ocean = land_mask(day, cc)
    boxes = {k: m & ocean for k, m in box_masks().items()}
    lows, gst, miss = [], [], []
    for lead in LEADS:
        u = url(day, cc, lead)
        try:
            for att in range(4):
                ib = get(u + ".idx")
                rows = idx_rows(ib.decode())
                if len(rows) < 500:                     # truncated index: retry, do not call it missing
                    time.sleep(2 ** att)
                    continue
                pb = message(u, rows, "PRMSL", "mean sea level")
                gb = message(u, rows, "GUST", "surface")
                try:
                    try:
                        msl = decode(pb, "prmsl") / 100.0
                    except Exception:
                        pb = message_walk(u, rows, "PRMSL", "mean sea level")
                        msl = decode(pb, "prmsl") / 100.0
                    try:
                        gust = decode(gb, "gust")
                    except Exception:
                        gb = message_walk(u, rows, "GUST", "surface")
                        gust = decode(gb, "gust")
                    break
                except Exception:
                    time.sleep(2 ** att)
            else:
                raise RuntimeError(f"unusable {u}")
        except Missing:
            miss.append(lead)
            continue
        for r in A.lows_in(msl, gust, ocean):
            lows.append((cyc, lead) + r)
        for b, m in boxes.items():
            x = gust[m] * KT
            gst.append((cyc, lead, b, round(float(x.max()), 2), round(float(np.percentile(x, 99.9)), 2)))
        with open(os.path.join(WORK, "bytes.log"), "a") as f:
            f.write(f"{cyc} f{lead:03d} {len(ib) + len(pb) + len(gb)}\n")
    for d in ("lows", "gust", "missing"):
        os.makedirs(os.path.join(WORK, d), exist_ok=True)
    with open(os.path.join(WORK, "gust", cyc + ".csv"), "w") as f:
        f.write("cycle,lead,basin,gmax_kt,g999_kt\n")
        f.writelines(",".join(map(str, r)) + "\n" for r in gst)
    with open(os.path.join(WORK, "missing", cyc + ".txt"), "w") as f:
        f.write(" ".join(map(str, miss)) + "\n")
    tmp = lo_p + ".tmp"
    with open(tmp, "w") as f:                        # written last: marks the cycle done
        f.write("cycle,lead,lat,lon,msl,g500,g800,g1200\n")
        f.writelines(",".join(map(str, r)) + "\n" for r in lows)
    os.replace(tmp, lo_p)
    return cyc, f"missing {miss}" if miss else "ok"


def days():
    for y in range(2021, 2026):
        d = date(y, 10, 1)
        while d <= date(y + 1, 4, 30):
            if d.month >= 10 or d.month <= 4:
                yield d.strftime("%Y%m%d")
            d += timedelta(1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--day")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    for d in ("lows", "gust", "missing"):
        os.makedirs(os.path.join(WORK, d), exist_ok=True)
    ds = [a.day] if a.day else sorted(set(days()))
    jobs = [(d, cc) for d in ds for cc in ("00", "12")]
    t0 = time.time()
    with Pool(min(a.workers, 4)) as p:
        for i, (c, s) in enumerate(p.imap_unordered(do_cycle, jobs)):
            if s != "skip" or i % 200 == 0:
                print(f"{c} {s} {i + 1}/{len(jobs)} {time.time() - t0:.0f}s", flush=True)
