"""Part 2 pull: GFS 0.25 degree fields for the Hart terms and the remaining environment terms of the FULL model, leads 0 and 24 h of every cycle in
../gfs_test/results/gfs_tracks_2021..2025.csv.gz, evaluated at every GFS forecast fix with intensity/env.py's own functions. Plan: PREREGISTRATION_FULL.md
and its Deviations. ERA5 pipeline A is a proxy; this is the GFS side of the FULL model.

Per file (byte-range records, method of gfs_test/a_fetch.py; stale .idx -> message_walk):
  HGT 925 850 700 600 500 400 300 mb, UGRD VGRD 300 mb, TMP 500 mb, PWAT (entire atmosphere), TMP surface, ICEC surface; lead 24 also LHTFL, SHTFL surface
  (18-24 h average); LAND surface once globally. HGT 850/500 and TMP 500 were pulled in part 1 but only a per-fix table was kept: re-pulled (Deviations 1).
Fields -> env.coarsen (1.5 degree) -> env.derived / env.features:
  B VTL VTU div300 sst sstgrad sst_t500 flux tcwv at each fix (heading from the previous 6 h of the track). SST = TMP surface over ocean (LAND == 0, ICEC < 0.15), NaN elsewhere.
  flux = LHTFL + SHTFL, upward-positive (no flux at lead 0: NaN). jet250 / vadv500 / eady come from part 1's table (their inputs are zero stand-ins here and discarded).
Work area GFS_WORK (default: gfs_work in the system temp folder): envfull/<cycle>.csv (written last = done, resumable), envfull_missing/<cycle>.txt, envfull_bytes.log (cycle fLLL bytes; .idx included).
Grids are never stored. HARD STOP when cumulative envfull_bytes.log exceeds 46 GB (Deviations 1).
usage: fetch_env_full.py [--workers N<=4] [--limit N] [--cycle YYYYMMDDCC]    ;    fetch_env_full.py --assemble OUT.csv.gz
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
import a_fetch as AF
import env as E
from eccodes import codes_get, codes_get_values, codes_new_from_message, codes_release

WORK = os.environ.get("GFS_WORK", os.path.join(__import__("tempfile").gettempdir(), "gfs_work"))
AF.WORK = WORK
LEADS = (0, 24)
HARD_STOP = 46e9
ATM = "entire atmosphere (considered as a single layer)"
# (var, level text in .idx, key, level number or None, accepted shortNames)
RECS = [("HGT", f"{p} mb", f"z{p}", p, ("gh",)) for p in (925, 850, 700, 600, 500, 400, 300)] + [
    ("UGRD", "300 mb", "u300", 300, ("u",)), ("VGRD", "300 mb", "v300", 300, ("v",)), ("TMP", "500 mb", "t500", 500, ("t",)),
    ("PWAT", ATM, "tcwv", None, ("pwat",)), ("TMP", "surface", "tsfc", None, ("t",)), ("ICEC", "surface", "icec", None, ("ci",))]
FLUX_RECS = [("LHTFL", "surface", "lhtfl", None, ("avg_slhtf", "slhtf", "lhtfl")), ("SHTFL", "surface", "shtfl", None, ("avg_ishf", "avg_sshf", "sshf", "shtfl"))]
LOG = os.path.join(WORK, "envfull_bytes.log")
COLS = ["B", "VTL", "VTU", "div300", "sst", "sstgrad", "sst_t500", "flux", "tcwv"]


def decode_any(blob):
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


def ok_msg(v, sn, lv, lev_num, names):
    assert sn in names, (sn, names)
    if lev_num is not None:
        assert lv == lev_num, (sn, lv)
    assert np.isfinite(v).all()


def read_rec(u, rows, var, level, lev_num, names):
    blob = AF.message(u, rows, var, level)
    nbytes = len(blob)
    try:
        v, sn, lv = decode_any(blob)
        ok_msg(v, sn, lv, lev_num, names)
    except Exception:
        blob = AF.message_walk(u, rows, var, level)
        nbytes += len(blob)
        v, sn, lv = decode_any(blob)
        ok_msg(v, sn, lv, lev_num, names)
    return v, nbytes


def land_ocean():
    p = os.path.join(WORK, "land.npy")
    if not os.path.exists(p):
        AF.land_mask("20240115", "00")          # writes land.npy (and a LAND line in bytes.log of WORK)
    return np.load(p) == 0


def cumulative():
    try:
        return sum(int(l.split()[-1]) for l in open(LOG) if l.strip())
    except FileNotFoundError:
        return 0


def fields(raw, ocean, lead):
    """raw 0.25 degree dict -> coarse dict for env.derived / env.features (zero stand-ins for the part 1 terms' inputs)."""
    ice = raw["icec"]
    oc = ocean & (ice < 0.15)
    F = {k: E.coarsen(raw[k])[E.R0:E.R1] for k in ("z925", "z850", "z700", "z600", "z500", "z400", "z300", "u300", "v300", "t500", "tcwv")}
    F["sst"] = E.coarsen(np.where(oc, raw["tsfc"], np.nan))[E.R0:E.R1]
    z = np.zeros_like(F["z500"])
    for k in ("u250", "v250", "u500", "v500", "u850", "v850"):
        F[k] = z
    F["t850"] = np.full_like(z, 280.0)
    if lead == 24:                                   # GFS fluxes are upward-positive; derived() negates, so hand it the downward-positive form
        F["lhf"] = -E.coarsen(raw["lhtfl"])[E.R0:E.R1]
        F["shf"] = -E.coarsen(raw["shtfl"])[E.R0:E.R1]
    else:
        F["lhf"] = F["shf"] = np.full_like(z, np.nan)
    return F


def do_cycle(args):
    try:
        return _cycle(args)
    except Exception as e:
        return args[0], f"ERROR {e!r}"


def _cycle(args):
    cyc, fixes = args
    out_p = os.path.join(WORK, "envfull", cyc + ".csv")
    if os.path.exists(out_p):
        return cyc, "skip"
    if cumulative() > HARD_STOP:
        return cyc, "HALT"
    day, cc = cyc[:8], cyc[8:10]
    ocean = land_ocean()
    rows_out, miss = [], []
    for lead in LEADS:
        fx = fixes[fixes.lead == lead]
        if not len(fx):
            continue
        recs = RECS + (FLUX_RECS if lead == 24 else [])
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
                    for var, level, key, ln, names in recs:
                        raw[key], b = read_rec(u, rows, var, level, ln, names)
                        nb += b
                    break
                except AF.Missing:
                    raise
                except Exception:
                    with open(LOG, "a") as f:            # bytes already spent on a failed attempt still count
                        f.write(f"{cyc} f{lead:03d} {nb}\n")
                    time.sleep(2 ** att)
            else:
                raise RuntimeError(f"unusable {u}")
        except AF.Missing:
            miss.append(lead)
            continue
        with open(LOG, "a") as f:
            f.write(f"{cyc} f{lead:03d} {nb}\n")
        F = fields(raw, ocean, lead)
        with np.errstate(all="ignore"):
            D = E.derived({k: v.astype(np.float64) for k, v in F.items()})
            for r in fx.itertuples():
                f = E.features(F, D, r.lat, r.lon, r.heading)
                rows_out.append((cyc, lead, r.track, r.lat, r.lon, r.heading) + tuple(f[c] for c in COLS))
    os.makedirs(os.path.join(WORK, "envfull_missing"), exist_ok=True)
    with open(os.path.join(WORK, "envfull_missing", cyc + ".txt"), "w") as f:
        f.write(" ".join(map(str, miss)) + "\n")
    tmp = out_p + ".tmp"
    pd.DataFrame(rows_out, columns=["cycle", "lead", "track", "lat", "lon", "heading"] + COLS).round(4).to_csv(tmp, index=False)
    os.replace(tmp, out_p)
    return cyc, f"missing {miss}" if miss else "ok"


def tracks():
    return pd.concat([pd.read_csv(os.path.join(GT, "results", f"gfs_tracks_{s}.csv.gz"), dtype={"cycle": str}).assign(season=s)
                      for s in range(2021, 2026)], ignore_index=True)


def cycles():
    T = tracks()
    P = T[["cycle", "track", "lead", "lat", "lon"]].assign(lead=lambda x: x.lead + 6).rename(columns={"lat": "lat_m6", "lon": "lon_m6"})
    T = T.merge(P, on=["cycle", "track", "lead"], how="left")
    dlat = T.lat - T.lat_m6
    dlon = ((T.lon - T.lon_m6 + 180) % 360 - 180) * np.cos(np.radians(T.lat))
    T["heading"] = np.degrees(np.arctan2(dlon, dlat)) % 360         # NaN where no fix 6 h earlier (intensity/fixes.py)
    T = T[T.lead.isin(LEADS)]
    return {c: g[["lead", "track", "lat", "lon", "heading"]] for c, g in T.groupby("cycle")}, T.cycle.nunique()


def assemble(outp):
    D = pd.concat([pd.read_csv(p, dtype={"cycle": str}) for p in sorted(glob.glob(os.path.join(WORK, "envfull", "*.csv")))], ignore_index=True)
    D.to_csv(outp, index=False, compression="gzip")
    print(len(D), "rows ->", outp)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--cycle")
    ap.add_argument("--assemble")
    a = ap.parse_args()
    if a.assemble:
        assemble(a.assemble)
        sys.exit()
    os.makedirs(os.path.join(WORK, "envfull"), exist_ok=True)
    land_ocean()
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
    print("HALTED: cumulative bytes exceeded 46 GB" if halted else "done", flush=True)
