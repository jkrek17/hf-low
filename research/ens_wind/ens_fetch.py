"""Pull 10 m sustained wind for GEFS, ECMWF ENS and AIFS ENS members, keep 1-degree cell maxima over ocean for the two basins.

Range requests for the two 10 m records only (the GRIB readers are the committed ones from tools/). Resumable: one npz per
(model, run date, lead). Work directory is outside the repository (ENS_WORK, default /tmp/ens_work). Plan: PREREGISTRATION.md."""
import os, sys, json, time, random, datetime as dt
from concurrent.futures import ThreadPoolExecutor
import numpy as np, requests

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import grib2_t53 as T53, grib2_ccsds as CC

WORK = os.environ.get("ENS_WORK", "/tmp/ens_work")
KT = 1.943844
GFS = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
GEFS = "https://noaa-gefs-pds.s3.amazonaws.com"
EC = "https://ecmwf-forecasts.s3.eu-central-1.amazonaws.com"
GEFS_MEMBERS = ["c00"] + ["p%02d" % i for i in range(1, 31)]
# whole-degree boxes (pipeline A basins); lon in 0..360
BOX = {"natl": dict(s=30, n=67, w=262, e=370), "npac": dict(s=27, n=67, w=135, e=240)}
LEADS = (24, 48, 72, 96)
_S = requests.Session(); _S.mount("https://", requests.adapters.HTTPAdapter(pool_connections=16, pool_maxsize=32))
_LAND = None


class Missing(Exception):
    pass


def get(url, rng=None, tries=8):
    h = {"Range": "bytes=%d-%s" % (rng[0], "" if rng[1] is None else rng[1])} if rng else None
    fails = slow = 0
    while fails < tries and slow < 80:
        try:
            r = _S.get(url, headers=h, timeout=180)
            if r.status_code in (200, 206): return r
            if r.status_code in (403, 404):
                fails += 3 if fails < 3 else 99    # transient proxy 403/404 happen: a few retries, then Missing
                time.sleep(1 + random.random()); continue
            if r.status_code == 503: slow += 1; time.sleep(0.4 + 1.2 * random.random()); continue
        except requests.RequestException: pass
        fails += 1; time.sleep(min(2 ** fails, 30) * (0.5 + random.random()))
    raise Missing(url)


def land_mask():
    """0.25 degree land flag (1 = land) from one GFS f000 file, cached."""
    global _LAND
    if _LAND is None:
        p = os.path.join(WORK, "land.npy")
        if not os.path.exists(p):
            url = GFS + "/gfs.20251115/00/atmos/gfs.t00z.pgrb2.0p25.f000"
            L = [l.split(":") for l in get(url + ".idx").text.splitlines() if l.strip()]
            i = next(j for j, f in enumerate(L) if f[3] == "LAND" and f[4] == "surface")
            m = T53.decode_message(get(url, (int(L[i][1]), int(L[i + 1][1]) - 1)).content)
            os.makedirs(WORK, exist_ok=True); np.save(p + ".tmp.npy", (m > 0.5)); os.replace(p + ".tmp.npy", p)
        _LAND = np.load(p)
    return _LAND


def cells(speed_kt):
    """global (721,1440) kt -> {basin: (ny,nx) uint8 of 2*kt}, ocean-only 4x4 max, north-up."""
    v = np.where(land_mask(), 0.0, speed_kt)
    out = {}
    for k, b in BOX.items():
        iy0 = (90 - b["n"]) * 4; ny = b["n"] - b["s"]
        ix = (np.arange((b["e"] - b["w"]) * 4) + b["w"] * 4) % 1440
        sub = v[iy0:iy0 + ny * 4][:, ix]
        out[k] = np.clip(np.round(sub.reshape(ny, 4, -1, 4).max(axis=(1, 3)) * 2), 0, 255).astype(np.uint8)
    return out


def spd(u, v):
    return cells(np.hypot(u, v) * KT)


def ncep_uv(url):
    L = [l.split(":") for l in get(url + ".idx").text.splitlines() if l.strip()]
    f = []
    for key in ("UGRD:10 m above ground", "VGRD:10 m above ground"):
        i = next(j for j, x in enumerate(L) if x[3] + ":" + x[4] == key)
        end = int(L[i + 1][1]) - 1 if i + 1 < len(L) else None
        r = get(url, (int(L[i][1]), end)); LOG.append(len(r.content))
        f.append(T53.decode_message(r.content))
    return spd(*f)


def gefs(ymd, step, m):
    return ncep_uv("%s/gefs.%s/00/atmos/pgrb2sp25/ge%s.t00z.pgrb2s.0p25.f%03d" % (GEFS, ymd, m, step))


def ec_base(model, stream, typ, ymd, step):
    return "%s/%s/00z/%s/0p25/%s/%s000000-%dh-%s-%s" % (EC, ymd, model, stream, ymd, step, stream, typ)


def ec_fields(base):
    recs = [json.loads(l) for l in get(base + ".index").text.splitlines() if l.strip().startswith("{")]
    jobs = [r for r in recs if r.get("levtype") == "sfc" and r["param"] in ("10u", "10v")]
    def one(r):
        c = get(base + ".grib2", (r["_offset"], r["_offset"] + r["_length"] - 1)).content; LOG.append(len(c))
        return int(r.get("number", 0)), r["param"], CC.decode_message(c)
    mem = {}
    with ThreadPoolExecutor(6) as ex:
        for m, p, a in ex.map(one, jobs): mem.setdefault(m, {})[p] = a
    return [spd(d["10u"], d["10v"]) for _, d in sorted(mem.items()) if len(d) == 2]


def ec_ens(model, ymd, step):
    out = []
    for typ in (("ef",) if model == "ifs" else ("cf", "pf")):
        out += ec_fields(ec_base(model, "enfo", typ, ymd, step))
    return out


LOG = []   # bytes of every response in this process


def task(model, ymd, step):
    """-> path of npz (members, ny, nx) per basin; None if the run is not posted."""
    p = os.path.join(WORK, model, "%s_%03d.npz" % (ymd, step))
    if os.path.exists(p): return p
    os.makedirs(os.path.dirname(p), exist_ok=True); land_mask()
    try:
        if model == "gefs":
            with ThreadPoolExecutor(8) as ex: mem = list(ex.map(lambda m: gefs(ymd, step, m), GEFS_MEMBERS))
        else:
            mem = ec_ens("ifs" if model == "ens" else "aifs-ens", ymd, step)
    except Missing as e:
        print("MISSING", model, ymd, step, e, flush=True); return None
    arr = {k: np.stack([m[k] for m in mem]) for k in BOX}
    np.savez_compressed(p + ".tmp.npz", **arr); os.replace(p + ".tmp.npz", p)
    with open(os.path.join(WORK, "bytes.log"), "a") as f: f.write("%s %s %d %d\n" % (model, ymd, step, sum(LOG)))
    LOG.clear()
    return p


def run_dates(season):
    """Every 5th day from 15 Oct to 15 Apr; season = start year."""
    d = dt.date(season, 10, 15); end = dt.date(season + 1, 4, 15); out = []
    while d <= end: out.append(d.strftime("%Y%m%d")); d += dt.timedelta(days=5)
    return out


def plan(model, seasons, leads):
    return [(model, d, s) for y in seasons for d in run_dates(y) for s in leads]


if __name__ == "__main__":
    import argparse, multiprocessing as mp
    a = argparse.ArgumentParser()
    a.add_argument("model", choices=["gefs", "ens", "aifs-ens"]); a.add_argument("--seasons", type=int, nargs="+")
    a.add_argument("--leads", type=int, nargs="+", default=list(LEADS)); a.add_argument("--procs", type=int, default=3)
    a.add_argument("--limit", type=int, default=0); a.add_argument("--count", action="store_true")
    o = a.parse_args()
    P = plan(o.model, o.seasons, o.leads)
    if o.limit: P = P[:o.limit]
    print(len(P), "tasks", flush=True)
    if o.count: sys.exit()
    land_mask()
    t0 = time.time()
    with mp.Pool(o.procs) as pool:
        for i, r in enumerate(pool.starmap(task, P, chunksize=1)):
            pass
    print("done", round(time.time() - t0), "s", flush=True)
