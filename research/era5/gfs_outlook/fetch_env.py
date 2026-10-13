"""Pull GFS 0.25 degree upper-air fields for leads 0 and 24 h of every cycle in ../gfs_test/results/gfs_tracks_2021..2025.csv.gz and
compute the three environment terms of intensity/env.py at every GFS forecast fix. Plan: PREREGISTRATION.md (ERA5 pipeline A is a proxy;
this is the GFS side of the A/B models).

Per file: .idx, then byte ranges for UGRD, VGRD (250/500/850 mb), TMP (500/850 mb), HGT (500/850 mb), 10 records, method of gfs_test/a_fetch.py.
Fields are decimated to the 1.5 degree grid with intensity/env.py coarsen() (the conservative 7-point cell average ERA5 pipeline A uses for
the 0.25 degree ARCO store), then env.derived() gives ws250, vadv500, eady, and env.py's own masks/wmean give, at each fix:
  jet250  max 250 hPa wind speed within 1000 km (kt)       vadv500  mean 500 hPa absolute vorticity advection within 500 km (1e-9 s-2)
  eady    mean 850-500 hPa Eady growth rate within 500 km (day-1)
(env.py is imported with gcsfs and numcodecs stubbed: its zarr readers are not used here; GFS HGT in gpm is the ERA5 z/g of env.py.)
Work area GFS_WORK (default /home/claude/gfs_work): env/<cycle>.csv (written last = cycle done, resumable), env_missing/<cycle>.txt,
env_bytes.log (one line per file: cycle fLLL bytes). Grids are never stored. Hard stop when cumulative env_bytes.log exceeds 40 GB.
usage: fetch_env.py [--workers N<=4] [--limit N] [--cycle YYYYMMDDCC]
"""
import os, sys, time, argparse, types, glob
import numpy as np, pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
GT = os.path.join(HERE, "..", "gfs_test")
sys.path.insert(0, GT)
sys.path.insert(0, os.path.join(HERE, "..", "intensity"))
for _m in ("gcsfs", "numcodecs"):
    sys.modules.setdefault(_m, types.ModuleType(_m))
import a_fetch as AF                    # get, idx_rows, message, message_walk, decode (decode asserts the 721 x 1440 grid)
import env as E                         # intensity/env.py
import hart
from eccodes import codes_get, codes_get_values, codes_new_from_message, codes_release

WORK = os.environ.get("GFS_WORK", "/home/claude/gfs_work")
LEADS = (0, 24)
HARD_STOP = 40e9
RECS = [("UGRD", "250 mb", "u250"), ("VGRD", "250 mb", "v250"), ("UGRD", "500 mb", "u500"), ("VGRD", "500 mb", "v500"),
        ("UGRD", "850 mb", "u850"), ("VGRD", "850 mb", "v850"), ("TMP", "500 mb", "t500"), ("TMP", "850 mb", "t850"),
        ("HGT", "500 mb", "z500"), ("HGT", "850 mb", "z850")]
LOG = os.path.join(WORK, "env_bytes.log")


def decode_any(blob):
    """GRIB message -> float32 [721,1440], lat 90..-90, lon 0..359.75 (checked as a_fetch.decode, minus the shortName check)."""
    g = codes_new_from_message(blob)
    try:
        ni, nj = int(codes_get(g, "Ni")), int(codes_get(g, "Nj"))
        assert (ni, nj) == (1440, 721), (ni, nj)
        v = np.array(codes_get_values(g), dtype=np.float32).reshape(nj, ni)
        assert not int(codes_get(g, "iScansNegatively")) and not int(codes_get(g, "jScansPositively"))
        assert abs(float(codes_get(g, "latitudeOfFirstGridPointInDegrees")) - 90) < 1e-6
        assert abs(float(codes_get(g, "longitudeOfFirstGridPointInDegrees"))) < 1e-6
        return v, codes_get(g, "shortName"), int(codes_get(g, "level"))
    finally:
        codes_release(g)


def read_rec(u, rows, var, level, lev_num):
    """Decoded field, bytes used. A stale .idx (offsets off) gives a bad or wrong-variable message: fall back to the section-0 walk."""
    blob = AF.message(u, rows, var, level)
    nbytes = len(blob)
    try:
        v, sn, lv = decode_any(blob)
        assert lv == lev_num, (sn, lv)
    except Exception:
        blob = AF.message_walk(u, rows, var, level)
        nbytes += len(blob)
        v, sn, lv = decode_any(blob)
        assert lv == lev_num, (sn, lv)
    return v, nbytes


def fix_features(F, D, la, lo):
    d = hart.great_circle_km(E.LAT2, E.LON2, la, lo)
    m5, m10 = d <= 500, d <= 1000
    return (float(np.nanmax(np.where(m10, D["ws250"], np.nan))), E.wmean(D["vadv500"], m5), E.wmean(D["eady"], m5))


def coarse_fields(raw):
    F = {}
    for k, v in raw.items():
        f = E.coarsen(v)[E.R0:E.R1]
        F[k] = f
    return F


def derived_three(F):
    """env.derived() with zero stand-ins for the fields these three terms do not use."""
    G = dict(F)
    z = np.zeros_like(F["z500"])
    for k in ("z925", "z700", "z600", "z400", "z300", "u300", "v300", "sst", "lhf", "shf"):
        G.setdefault(k, z)
    return E.derived(G)


def cumulative():
    try:
        return sum(int(l.split()[-1]) for l in open(LOG) if l.strip())
    except FileNotFoundError:
        return 0


def do_cycle(args):
    try:
        return _cycle(args)
    except Exception as e:
        return args[0], f"ERROR {e!r}"


def _cycle(args):
    cyc, fixes = args
    out_p = os.path.join(WORK, "env", cyc + ".csv")
    if os.path.exists(out_p):
        return cyc, "skip"
    if cumulative() > HARD_STOP:
        return cyc, "HALT"
    day, cc = cyc[:8], cyc[8:10]
    rows_out, miss = [], []
    for lead in LEADS:
        fx = fixes[fixes.lead == lead]
        if not len(fx):
            continue
        u = AF.url(day, cc, lead)
        try:
            for att in range(4):
                ib = AF.get(u + ".idx")
                rows = AF.idx_rows(ib.decode())
                if len(rows) < 500:
                    time.sleep(2 ** att)
                    continue
                raw, nb = {}, len(ib)
                try:
                    for var, level, name in RECS:
                        raw[name], b = read_rec(u, rows, var, level, int(level.split()[0]))
                        nb += b
                    break
                except AF.Missing:
                    raise
                except Exception:
                    time.sleep(2 ** att)
            else:
                raise RuntimeError(f"unusable {u}")
        except AF.Missing:
            miss.append(lead)
            continue
        with open(LOG, "a") as f:
            f.write(f"{cyc} f{lead:03d} {nb}\n")
        F = coarse_fields(raw)
        D = derived_three(F)
        for r in fx.itertuples():
            j, v, e = fix_features(F, D, r.lat, r.lon)
            rows_out.append((cyc, lead, r.track, r.lat, r.lon, round(j, 3), round(v, 4), round(e, 4)))
    os.makedirs(os.path.join(WORK, "env_missing"), exist_ok=True)
    with open(os.path.join(WORK, "env_missing", cyc + ".txt"), "w") as f:
        f.write(" ".join(map(str, miss)) + "\n")
    tmp = out_p + ".tmp"
    with open(tmp, "w") as f:
        f.write("cycle,lead,track,lat,lon,jet250,vadv500,eady\n")
        f.writelines(",".join(map(str, r)) + "\n" for r in rows_out)
    os.replace(tmp, out_p)
    return cyc, f"missing {miss}" if miss else "ok"


def cycles():
    T = pd.concat([pd.read_csv(os.path.join(GT, "results", f"gfs_tracks_{s}.csv.gz"), dtype={"cycle": str}) for s in range(2021, 2026)],
                  ignore_index=True)
    T = T[T.lead.isin(LEADS)]
    return {c: g[["lead", "track", "lat", "lon"]] for c, g in T.groupby("cycle")}, T.cycle.nunique()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--cycle")
    a = ap.parse_args()
    os.makedirs(os.path.join(WORK, "env"), exist_ok=True)
    C, n = cycles()
    jobs = sorted(C.items())
    if a.cycle:
        jobs = [j for j in jobs if j[0] == a.cycle]
    if a.limit:
        jobs = jobs[:a.limit]
    print(f"cycles with fixes at lead 0/24: {n}; jobs {len(jobs)}; cumulative bytes so far {cumulative()}", flush=True)
    t0 = time.time()
    halted = False
    with Pool(min(a.workers, 4)) as p:
        for i, (c, s) in enumerate(p.imap_unordered(do_cycle, jobs)):
            if s == "HALT":
                halted = True
            if s != "skip" or i % 200 == 0:
                print(f"{c} {s} {i + 1}/{len(jobs)} {time.time() - t0:.0f}s cum {cumulative() / 1e9:.2f} GB", flush=True)
    print("HALTED: cumulative bytes exceeded 40 GB" if halted else "done", flush=True)
