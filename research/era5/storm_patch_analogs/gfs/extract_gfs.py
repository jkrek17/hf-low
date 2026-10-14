"""GFS f024 patches (HGT 500 mb, PRMSL, UGRD/VGRD 250 mb) at the matched lead 24 fixes of results/gfs_fixes_lead24.csv.gz.
usage: GFS_WORK=DIR extract_gfs.py [NPROC]     per cycle one npz in DIR/gfs_patches (resumable); bytes in DIR/bytes_gfs.log; stops at 9 GB.
Same patch recipe as the ERA5 patches (north-up, 13 x 13, 154.33 km) after env.coarsen to 1.5 degrees; anomalies vs clim_patch_0012.npz."""
import os, sys, glob, types
HERE = os.path.dirname(os.path.abspath(__file__)); RE = os.path.join(HERE, "..", "..")
WORK = os.environ.get("GFS_WORK") or os.path.join(HERE, "work"); os.environ["GFS_WORK"] = WORK; os.makedirs(WORK, exist_ok=True)
sys.path.insert(0, os.path.join(RE, "gfs_test")); sys.path.insert(0, os.path.join(RE, "intensity")); sys.path.insert(0, os.path.join(HERE, ".."))
import numpy as np, pandas as pd
from multiprocessing import Pool
import a_fetch as AF
import env as E
from common_patch import PR as _PR, PBANG, sample, destination, CH
from eccodes import codes_get, codes_get_values, codes_new_from_message, codes_release
RECS = [("HGT", "500 mb", 500, "z"), ("PRMSL", "mean sea level", 0, "m"), ("UGRD", "250 mb", 250, "u"), ("VGRD", "250 mb", 250, "v")]
HARD = 9e9
_C = None


def clim():
    global _C
    if _C is None:
        z = np.load(os.path.join(os.environ.get("STORM_PATCH_DIR", "/mnt/project-files/storm-patches"), "clim_patch_0012.npz")); _C = {k: z[k] for k in z.files}
    return _C


def decode_any(blob):
    g = codes_new_from_message(blob)
    try:
        assert (int(codes_get(g, "Ni")), int(codes_get(g, "Nj"))) == (1440, 721)
        v = np.array(codes_get_values(g), dtype=np.float32).reshape(721, 1440)
        assert not int(codes_get(g, "iScansNegatively")) and not int(codes_get(g, "jScansPositively"))
        return v, codes_get(g, "shortName")
    finally:
        codes_release(g)


def job(args):
    cyc, rows = args
    path = os.path.join(WORK, "gfs_patches", f"{cyc}.npz")
    if os.path.exists(path):
        return cyc, 0, "skip"
    spent = sum(float(l.split()[-1]) for l in open(os.path.join(WORK, "bytes_gfs.log"))) if os.path.exists(os.path.join(WORK, "bytes_gfs.log")) else 0
    if spent > HARD:
        raise SystemExit(f"stop: {spent / 1e9:.1f} GB")
    day, cc = cyc[:8], cyc[8:]
    u = AF.url(day, cc, 24)
    try:
        idx = AF.get(u + ".idx"); ir = AF.idx_rows(idx.decode()); nb = len(idx)
        F = {}
        for var, lev, _, tag in RECS:
            blob = AF.message(u, ir, var, lev); nb += len(blob)
            v, sn = decode_any(blob)
            F[tag] = E.coarsen(v)
    except AF.Missing:
        open(os.path.join(WORK, "missing.txt"), "a").write(cyc + "\n"); return cyc, 0, "missing"
    when = np.datetime64(f"{day[:4]}-{day[4:6]}-{day[6:8]}T{cc}") + np.timedelta64(24, "h")
    h = int((when.astype("datetime64[h]") - when.astype("datetime64[D]")).astype(int)) // 12
    doy = int((when.astype("datetime64[D]") - when.astype("datetime64[Y]").astype("datetime64[D]")).astype(int)) + 1
    C = clim(); assert int(C["hour"][h]) == h * 12 and C["doy"][doy - 1] == doy
    Fd = {"a_z500": F["z"] - C["z500"][h, doy - 1], "a_msl": F["m"] / 100.0 - C["msl"][h, doy - 1],
          "a_u250": F["u"] - C["u250"][h, doy - 1], "a_v250": F["v"] - C["v250"][h, doy - 1], "z500": F["z"], "msl": F["m"] / 100.0}
    out, keys = [], []
    for r in rows.itertuples():
        la, lo = destination(r.lat, r.lon, 90.0 + PBANG, _PR)
        out.append(np.stack([sample(Fd[c], la, lo) for c in CH]).astype(np.float16)); keys.append(f"{r.cycle}|{r.track}")
    np.savez_compressed(path + ".tmp.npz", keys=np.array(keys), x=np.stack(out))
    os.replace(path + ".tmp.npz", path)
    with open(os.path.join(WORK, "bytes_gfs.log"), "a") as f:
        f.write(f"{cyc} {nb}\n")
    return cyc, nb, "ok"


if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    os.makedirs(os.path.join(WORK, "gfs_patches"), exist_ok=True)
    S = pd.read_csv(os.path.join(HERE, "results", "gfs_fixes_lead24.csv.gz"), dtype={"cycle": str})
    J = [(c, g) for c, g in S.groupby("cycle")]
    print(len(J), "cycles,", len(S), "fixes; expected ~", round(len(J) * 2.97 / 1000, 1), "GB", flush=True)
    clim(); tot = 0
    with Pool(nproc) as p:
        for i, (c, nb, st) in enumerate(p.imap_unordered(job, J)):
            tot += nb
            if i % 100 == 0:
                print(f"{i + 1}/{len(J)} {c} {st} pulled {tot / 1e9:.2f} GB", flush=True)
    print("done, pulled this run", round(tot / 1e9, 2), "GB")
