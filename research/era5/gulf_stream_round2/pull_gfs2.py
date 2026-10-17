"""Round 2 extra GFS fields (2 m dewpoint, surface pressure, surface momentum flux UFLX/VFLX 18-24 h average) for the same times as pull_gfs.py. Original docstring: GFS 0.25 degree 24 h forecasts valid 00 and 12 UTC (init 24 h earlier), Oct-Apr winters 2021-22..2025-26: 10 m U/V, 2 m T, surface T,
fetched by HTTP range from the NOAA GFS AWS bucket, cropped to 25-45N, 85-50W, int16 (m/s*100, (K-200)*100). Resumable per winter. About 7.1 GB.
Usage: python3 pull_gfs.py OUT [THREADS]"""
import os, sys, time, datetime as dt, tempfile, concurrent.futures as cf, requests, numpy as np, eccodes
OUT = sys.argv[1]; NT = int(sys.argv[2]) if len(sys.argv) > 2 else 8
KEEP = {("DPT", "2 m above ground"): "dpt", ("PRES", "surface"): "pres", ("UFLX", "surface"): "uflx", ("VFLX", "surface"): "vflx"}
S = requests.Session(); S.mount("https://", requests.adapters.HTTPAdapter(pool_maxsize=NT * 2))
def http(url, rng=None, tries=12):
    for k in range(tries):
        try:
            r = S.get(url, headers={"Range": f"bytes={rng}"} if rng else {}, timeout=120)
            if r.status_code in (403, 404): return None
            r.raise_for_status(); return r.content
        except Exception:
            time.sleep(min(2 ** k * 0.5, 20))
    return None
def one(valid):
    init = valid - dt.timedelta(hours=24)
    base = f"https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.{init:%Y%m%d}/{init:%H}/atmos/gfs.t{init:%H}z.pgrb2.0p25.f024"
    ix = http(base + ".idx")
    if ix is None: return valid, None
    L = [l.split(":") for l in ix.decode().strip().split("\n")]; got = {}
    for i, r in enumerate(L):
        k = (r[3], r[4])
        if k in KEEP:
            end = int(L[i + 1][1]) - 1 if i + 1 < len(L) else ""
            b = http(base, f"{r[1]}-{end}")
            if b is None: return valid, None
            with tempfile.NamedTemporaryFile(suffix=".grb") as f:
                f.write(b); f.flush()
                with open(f.name, "rb") as fh:
                    g = eccodes.codes_grib_new_from_file(fh)
                    v = eccodes.codes_get_values(g).reshape(int(eccodes.codes_get(g, "Nj")), int(eccodes.codes_get(g, "Ni")))
                    la1 = eccodes.codes_get(g, "latitudeOfFirstGridPointInDegrees"); eccodes.codes_release(g)
            if la1 < 0: v = v[::-1]
            got[KEEP[k]] = v[180:261, 1100:1241]                    # 45N..25N rows, 275E..310E
    if len(got) != 4: return valid, None
    enc = lambda a, o, s: np.clip(np.rint((a - o) * s), -32768, 32767).astype(np.int16)
    return valid, np.stack([enc(got["dpt"], 200, 100), enc(got["pres"], 80000, 0.1), enc(got["uflx"], 0, 1000), enc(got["vflx"], 0, 1000)])
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for w in range(2021, 2026):
        fn = f"{OUT}/gfs2_f024_{w}.npz"
        if os.path.exists(fn): continue
        vs = []; t = dt.datetime(w, 10, 1)
        while t < dt.datetime(w + 1, 5, 1): vs += [t, t + dt.timedelta(hours=12)]; t += dt.timedelta(days=1)
        with cf.ThreadPoolExecutor(NT) as ex: res = list(ex.map(one, vs))
        ok = [(v, a) for v, a in res if a is not None]
        np.savez_compressed(fn + ".tmp.npz", valid=np.array([int(v.strftime("%Y%m%d%H")) for v, _ in ok]), data=np.stack([a for _, a in ok]))
        os.rename(fn + ".tmp.npz", fn); print(w, len(ok), "of", len(vs), flush=True)
