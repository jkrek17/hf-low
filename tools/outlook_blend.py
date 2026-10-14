"""GFE-style blend of 10 m sustained wind for the Outlook (second layer of the live page, built by tools/outlook_blend_live.py).

Recipe, ported from jkrek17/awips-tools marine_wind/pipeline/build.py (read-only reference):
    base   GFS, IFS (ECMWF) and AIFS 10 m speeds, weighted at each point by member support: the share of the
           GEFS + ENS members within about 5 kt of each model (Gaussian kernel, SUPPORT_KT wide)
    wind   base + min(k x member spread, 5 kt); k = 0.25, 0.15 at +24 h and earlier
    p34/p48/p64  share of GEFS + ENS + AIFS-ENS members at 34/48/64 kt or more; p48 and p64 recalibrated with
           data/outlook_live/gfe_calibration.json (the GFE page's hindcast table). These are chances that the wind
           occurs, judged against buoys and ships, NOT the "OPC agreed" rates of the gust shading.
Grids are the GFE page's: natl 18-67N 95W-10E, npac 18-67N 135E-115W, 0.25 degrees. Sustained 10 m wind, not gust.
Decoding uses tools/grib2_t53.py (GFS, GEFS) and tools/grib2_ccsds.py (ECMWF, through the system libaec), because
eccodes cannot be installed here.

    python3 tools/outlook_blend.py --cycle 2026101400 --steps 24 --out OUTDIR    (keeps the 0.25 degree steps; the daily job is tools/outlook_blend_live.py)
"""
import argparse, base64, json, os, random, sys, time, zlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np, requests

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import grib2_t53 as T53, grib2_ccsds as CC

KT = 1.943844
GFS = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
GEFS = "https://noaa-gefs-pds.s3.amazonaws.com"
EC = "https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com"
GEFS_MEMBERS = ["c00"] + ["p%02d" % i for i in range(1, 31)]
CAP_KT, SUPPORT_KT, MIN_MEMBERS = 5.0, 5.0, 20
BLEND_FIELDS = ("wind", "mean", "sigma", "boost", "gfs", "ecmwf", "aifs", "p34", "p48", "p64", "p48raw", "p64raw")
STEPS = list(range(12, 97, 12))
BASINS = {"natl": dict(north=67.0, south=18.0, west=-95.0, east=10.0),
          "npac": dict(north=67.0, south=18.0, west=135.0, east=245.0)}
for B in BASINS.values():
    B["iy"] = slice(int((90 - B["north"]) * 4), int((90 - B["south"]) * 4) + 1)
    B["ix"] = (int(round((B["west"] % 360) * 4)) + np.arange(int(round((B["east"] - B["west"]) * 4)) + 1)) % 1440
CALIB = json.load(open(os.path.join(HERE, "..", "data", "outlook_live", "gfe_calibration.json")))
_S = requests.Session(); _S.mount("https://", requests.adapters.HTTPAdapter(pool_connections=16, pool_maxsize=32))


class Missing(Exception):
    pass


def get(url, rng=None, tries=8):
    """GET with optional byte range. 403/404 -> Missing (not posted). 503 SlowDown retried fast, other errors back off."""
    h = {"Range": "bytes=%d-%s" % (rng[0], "" if rng[1] is None else rng[1])} if rng else None
    fails = slow = 0
    while fails < tries and slow < 80:
        try:
            r = _S.get(url, headers=h, timeout=180)
            if r.status_code in (200, 206): return r
            if r.status_code in (403, 404): raise Missing(url)
            if r.status_code == 503: slow += 1; time.sleep(0.4 + 1.2 * random.random()); continue
        except requests.RequestException: pass
        fails += 1; time.sleep(min(2 ** fails, 30) * (0.5 + random.random()))
    raise RuntimeError("failed: " + url)


def crop(v):
    """Global (721, 1440) 0..359.75 north-up -> {basin: crop} in kt-ready float32."""
    return {k: v[B["iy"]][:, B["ix"]].copy() for k, B in BASINS.items()}


def spd(u, v):
    return {k: (np.hypot(u[k], v[k]) * KT).astype(np.float32) for k in u}


def ncep_uv(url):
    L = [l.split(":") for l in get(url + ".idx").text.splitlines() if l.strip()]
    out = []
    for key in ("UGRD:10 m above ground", "VGRD:10 m above ground"):
        i = next(j for j, f in enumerate(L) if f[3] + ":" + f[4] == key)
        end = int(L[i + 1][1]) - 1 if i + 1 < len(L) else None
        out.append(crop(T53.decode_message(get(url, (int(L[i][1]), end)).content)))
    return spd(*out)


def gfs(ymd, hh, step):
    return ncep_uv("%s/gfs.%s/%s/atmos/gfs.t%sz.pgrb2.0p25.f%03d" % (GFS, ymd, hh, hh, step))


def gefs(ymd, hh, step, m):
    return ncep_uv("%s/gefs.%s/%s/atmos/pgrb2sp25/ge%s.t%sz.pgrb2s.0p25.f%03d" % (GEFS, ymd, hh, m, hh, step))


def ec_base(model, stream, typ, ymd, hh, step):
    return "%s/%s/%sz/%s/0p25/%s/%s%s0000-%dh-%s-%s" % (EC, ymd, hh, model, stream, ymd, hh, step, stream, typ)


def ec_records(base):
    return [json.loads(l) for l in get(base + ".index").text.splitlines() if l.strip().startswith("{")]


def ec_fields(base, recs):
    """{member: {'10u','10v'}} for the sfc 10u/10v records of one file."""
    jobs = [r for r in recs if r.get("levtype") == "sfc" and r["param"] in ("10u", "10v")]
    def one(r):
        return int(r.get("number", 0)), r["param"], crop(CC.decode_message(get(base + ".grib2", (r["_offset"], r["_offset"] + r["_length"] - 1)).content))
    mem = {}
    with ThreadPoolExecutor(6) as ex:
        for m, p, a in ex.map(one, jobs): mem.setdefault(m, {})[p] = a
    return [spd(d["10u"], d["10v"]) for _, d in sorted(mem.items()) if len(d) == 2]


def ec_det(model, ymd, hh, step):
    base = ec_base(model, "oper", "fc", ymd, hh, step)
    r = ec_fields(base, ec_records(base))
    if len(r) != 1: raise Missing(base)
    return r[0]


def ec_ens(model, ymd, hh, step):
    out = []
    for typ in (("ef",) if model == "ifs" else ("cf", "pf")):
        base = ec_base(model, "enfo", typ, ymd, hh, step)
        out += ec_fields(base, ec_records(base))
    return out


def k_for(step):
    return 0.15 if step <= 24 else 0.25


def support_base(speeds, members):
    w = np.array([np.exp(-0.5 * ((members - s[None]) / SUPPORT_KT) ** 2).mean(0) for s in speeds]) + 1e-6
    return (w * speeds).sum(0) / w.sum(0)


def calibrate(pct, thr, basin, step):
    g = "le48" if step <= 48 else ("72" if step <= 84 else "ge96")
    t = CALIB["thresholds"][str(thr)][basin][g]
    return np.interp(pct, t["raw"], t["cal"])


def blend_step(ymd, hh, step):
    """{basin: {'wind','mean','sigma','p34','p48','p64', ...}} plus counts, for one valid time."""
    with ThreadPoolExecutor(3) as ex:
        det_f = [ex.submit(gfs, ymd, hh, step), ex.submit(ec_det, "ifs", ymd, hh, step), ex.submit(ec_det, "aifs-single", ymd, hh, step)]
        ens_f = ex.submit(ec_ens, "ifs", ymd, hh, step)
        aens_f = ex.submit(ec_ens, "aifs-ens", ymd, hh, step)
        with ThreadPoolExecutor(8) as gx:
            gefs_f = [gx.submit(gefs, ymd, hh, step, m) for m in GEFS_MEMBERS]
            dets = [f.result() for f in det_f]
            ge = [f.result() for f in gefs_f]
        en = ens_f.result()
        try: aens = aens_f.result()
        except Missing: aens = []
    mem = ge + en
    if len(mem) < MIN_MEMBERS: raise RuntimeError("only %d members at +%dh" % (len(mem), step))
    out = {}
    for b in BASINS:
        speeds = np.array([d[b] for d in dets]); members = np.array([m[b] for m in mem])
        pmem = members if not aens else np.concatenate([members, np.array([m[b] for m in aens])])
        mean = support_base(speeds, members); sigma = members.std(0)
        boost = np.minimum(k_for(step) * sigma, CAP_KT)
        p34 = 100 * (pmem >= 34).mean(0)
        out[b] = dict(wind=mean + boost, mean=mean, sigma=sigma, boost=boost, gfs=speeds[0], ecmwf=speeds[1], aifs=speeds[2],
                      p34=p34, p48=calibrate(100 * (pmem >= 48).mean(0), 48, b, step), p64=calibrate(100 * (pmem >= 64).mean(0), 64, b, step),
                      p48raw=100 * (pmem >= 48).mean(0), p64raw=100 * (pmem >= 64).mean(0))
    n = dict(step=step, members=len(mem), gefs=len(ge), ens=len(en), aifs_ens=len(aens), prob_members=len(mem) + len(aens))
    return out, n


DAY_STEPS = {1: (12, 24), 2: (36, 48), 3: (60, 72), 4: (84, 96)}
DAY_FIELDS = ("wind", "gfs", "ecmwf", "aifs", "sigma", "p34", "p48", "p64")
SCALE = {"wind": 2, "gfs": 2, "ecmwf": 2, "aifs": 2, "sigma": 4, "p34": 1, "p48": 1, "p64": 1}   # stored = round(value x scale), uint8
CELL = 4            # 0.25 degree cells per page cell (1 degree)
LAT_MIN = 20.0      # the page's map shows 20N and north


def run_cycle(ymd, hh, steps=STEPS, cache=None, log=print):
    """{step: {basin: {field: 0.25 degree array}}}. With cache (a folder) each step is saved as fNNN.npz and reused."""
    res = {}
    for step in steps:
        f = os.path.join(cache, "f%03d.npz" % step) if cache else None
        if f and os.path.exists(f):
            z = np.load(f); res[step] = {b: {k: z["%s_%s" % (b, k)].astype(np.float32) for k in BLEND_FIELDS} for b in BASINS}; continue
        t = time.time(); out, n = blend_step(ymd, hh, step)
        log("+%dh %s in %.0fs; max wind natl %.1f npac %.1f kt" % (step, n, time.time() - t, out["natl"]["wind"].max(), out["npac"]["wind"].max()))
        res[step] = out
        if f:
            os.makedirs(cache, exist_ok=True)
            np.savez_compressed(f, **{"%s_%s" % (b, k): v.astype(np.float16 if k.startswith("p") else np.float32) for b, d in out.items() for k, v in d.items()})
            json.dump(n, open(os.path.join(cache, "f%03d.json" % step), "w"))
    return res


def pool(a, n=CELL):
    """Max of each n x n block (edges padded by repeating the last row/column)."""
    a = np.pad(a, ((0, -a.shape[0] % n), (0, -a.shape[1] % n)), mode="edge")
    return a.reshape(a.shape[0] // n, n, a.shape[1] // n, n).max(axis=(1, 3))


def day_grids(results):
    """Per basin and forecast day 1-4 the higher of the day's two 12-hourly valid times, as the highest value in each 1 degree cell, packed to uint8
    (see SCALE). Returns {basin: {'lat0','lon0','res','ny','nx','days': {day: {field: uint8 array [ny, nx]}}}}; a day is left out when one of its steps is missing."""
    out = {}
    for b, B in BASINS.items():
        ny = int((B["north"] - LAT_MIN) * 4 // CELL); days = {}
        for d, steps in DAY_STEPS.items():
            if not all(s in results for s in steps): continue
            rec = {}
            for k in DAY_FIELDS:
                a = pool(np.max([results[s][b][k].astype(np.float32) for s in steps], axis=0))[:ny]
                rec[k] = np.clip(np.rint(a * SCALE[k]), 0, 255).astype(np.uint8)
            days[d] = rec
        if days:
            out[b] = dict(lat0=B["north"], lon0=B["west"], res=1.0, ny=ny, nx=next(iter(days.values()))["wind"].shape[1], days=days)
    return out


def pack_live(results, run, made):
    """The page's data object: per basin (atl = natl, pac = npac) and day, one base64 deflate blob of the eight fields back to back (field-major)."""
    g = day_grids(results); P = {"run": run, "made": made, "res": 1.0, "lat0": 67.0, "fields": list(DAY_FIELDS), "scale": [SCALE[k] for k in DAY_FIELDS], "basins": {}}
    for b, key in (("natl", "atl"), ("npac", "pac")):
        v = g.get(b)
        if v is None: continue
        P["basins"][key] = dict(lon0=v["lon0"], ny=v["ny"], nx=v["nx"],
            days=[(base64.b64encode(zlib.compress(b"".join(v["days"][d][k].tobytes() for k in DAY_FIELDS), 9)).decode() if d in v["days"] else None) for d in (1, 2, 3, 4)])
    return P


def ready(ymd, hh, step=max(STEPS)):
    """True when every source has its last step posted (HEAD on the index of each; GEFS control stands for the members)."""
    urls = [GFS + "/gfs.%s/%s/atmos/gfs.t%sz.pgrb2.0p25.f%03d.idx" % (ymd, hh, hh, step),
            GEFS + "/gefs.%s/%s/atmos/pgrb2sp25/gep30.t%sz.pgrb2s.0p25.f%03d.idx" % (ymd, hh, hh, step),
            ec_base("ifs", "oper", "fc", ymd, hh, step) + ".index", ec_base("aifs-single", "oper", "fc", ymd, hh, step) + ".index",
            ec_base("ifs", "enfo", "ef", ymd, hh, step) + ".index", ec_base("aifs-ens", "enfo", "cf", ymd, hh, step) + ".index",
            ec_base("aifs-ens", "enfo", "pf", ymd, hh, step) + ".index"]
    for u in urls:
        for _ in range(6):
            try:
                r = _S.head(u, timeout=60)
            except requests.RequestException:
                time.sleep(2); continue
            if r.status_code == 200: break
            if r.status_code in (403, 404): return False
            time.sleep(1)
        else:
            return False
    return True


def main():
    ap = argparse.ArgumentParser(description="Pull one cycle and keep the 0.25 degree steps (testing); the daily job uses tools/outlook_blend_live.py")
    ap.add_argument("--cycle", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--steps", default=",".join(map(str, STEPS)))
    a = ap.parse_args(); run_cycle(a.cycle[:8], a.cycle[8:], [int(s) for s in a.steps.split(",")], cache=a.out)


if __name__ == "__main__":
    main()
