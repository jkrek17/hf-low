"""Where pipeline A's gust index comes from, at Atlantic HF-strength fixes.

ERA5 proxy. Pipeline A (research/era5/hf_history) scores a low by the maximum
ERA5 instantaneous 10 m gust over ocean within 800 km, counting only grid
points nearer this low than any other. This script re-runs that exact
detection and ownership at every 6-hourly time when an Atlantic pipeline A
event track (seasons 2004-05 on) has an in-domain fix at or above the 71.7 kt
threshold, and records for that fix:

  - g800 recomputed (must equal the catalog value; checked in analyse.py)
  - where the maximum sits: lat, lon, distance and bearing from the centre,
    distance to Greenland, to Iceland and to any land, sea-ice concentration
    and the 10 m wind direction there
  - the index with terrain or ice masked out:
      g800_gl100  ocean points within 100 km of Greenland removed
      g800_gl300  ocean points within 300 km of Greenland removed
      g800_land50 ocean points within 50 km of any land removed
      g800_ice    points with sea-ice concentration > 0.15 removed
      g800_r400   only points within 400 km of the centre (a centre-driven
                  wind should survive this; a remote terrain jet will not)

Sea ice and 10 m wind are fetched only at times with a target fix at or north
of 55N, to keep the pull small. The ERA5 land-sea mask (< 0.5 = ocean) is the
same one pipeline A uses, so sea-ice-covered ocean counts as ocean in A.

Pull size (HEAD-checked 2026-10-08 on a 2015 chunk): gust 3.2 MB, MSLP 2.2 MB,
sea ice 0.4 MB, u10 and v10 3.3 MB each. About 3,400 times x 5.4 MB plus about
1,600 times x 7.0 MB = about 29 GB streamed, nothing stored but the output.

usage: gustloc.py [NPROC]   writes work/gustloc/<time>.csv, resumable;
                            then analyse.py combines them.
"""
import os, sys, time
import numpy as np, pandas as pd, numcodecs, gcsfs
from scipy import ndimage
from scipy.spatial import cKDTree
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
OUT = os.path.join(WORK, "gustloc")
RES = os.path.join(HERE, "..", "hf_history", "results")

ROOT = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
T0 = np.datetime64("1900-01-01T00")
R0, R1 = 60, 281                    # rows for 75N .. 20N, as extract.py
LAT = 90 - 0.25 * np.arange(R0, R1)
LON = 0.25 * np.arange(1440)
RE = 6371.0
KT = 1 / 0.514444
THR = 71.7
FIRST_SEASON = 2004
WIND_LAT = 55.0

_codec = numcodecs.Blosc(cname="lz4", clevel=5, shuffle=1)
_fs = None


def fs():
    global _fs
    if _fs is None:
        _fs = gcsfs.GCSFileSystem(token="anon")
    return _fs


def get(var, ti, full=False):
    for k in range(6):
        try:
            b = fs().cat(f"{ROOT}/{var}/{ti}.0.0")
            a = np.frombuffer(_codec.decode(b), "<f4").reshape(721, 1440)
            return a if full else a[R0:R1]
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {var} {ti}")


# ---- static fields: land-sea mask, land masses, distances -------------------
os.makedirs(WORK, exist_ok=True)
LSM_PATH = os.path.join(WORK, "lsm_full.npy")
if not os.path.exists(LSM_PATH):
    np.save(LSM_PATH, get("land_sea_mask", 876576, full=True))     # 2000-01-01 00Z, as extract.py
LSM_FULL = np.load(LSM_PATH)
LSM = LSM_FULL[R0:R1] < 0.5                                         # ocean, as extract.py
LA2, LO2 = np.meshgrid(np.radians(LAT), np.radians(LON), indexing="ij")
XYZ = np.stack([np.cos(LA2) * np.cos(LO2), np.cos(LA2) * np.sin(LO2), np.sin(LA2)], -1)


def _landmass(lat, lon):
    """Connected land component (LSM >= 0.5, full globe) containing a point."""
    land = LSM_FULL >= 0.5
    lab, _ = ndimage.label(land)
    return lab == lab[int(round((90 - lat) / 0.25)), int(round(lon / 0.25)) % 1440]


def _dist_to(mask_full):
    """km from every grid point in rows R0:R1 to the nearest True cell of a full-globe mask."""
    lat = np.radians(90 - 0.25 * np.arange(721))
    la, lo = np.meshgrid(lat, np.radians(LON), indexing="ij")
    xyz = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)
    tree = cKDTree(xyz[mask_full])
    d, _ = tree.query(XYZ.reshape(-1, 3), k=1)
    return (2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE).reshape(LSM.shape)


DIST_PATH = os.path.join(WORK, "dist_land.npz")
if not os.path.exists(DIST_PATH):
    gl = _landmass(72.0, 320.0)          # Greenland
    ic = _landmass(64.8, 342.0)          # Iceland
    assert gl.sum() > 1000 and ic.sum() > 50 and not (gl & ic).any()
    np.savez(DIST_PATH, gl=_dist_to(gl), ic=_dist_to(ic), land=_dist_to(LSM_FULL >= 0.5))
_D = np.load(DIST_PATH)
D_GL, D_IC, D_LAND = _D["gl"], _D["ic"], _D["land"]


# ---- pipeline A detection, copied from hf_history/extract.py ----------------
def lows(msl):
    sm = ndimage.gaussian_filter(msl, 2, mode=["nearest", "wrap"])
    mn = ndimage.minimum_filter(sm, size=(21, 21), mode=["nearest", "wrap"])
    r, c = np.where((sm == mn) & (sm < 1010.0))
    keep = (r > 2) & (r < len(LAT) - 3)
    r, c = r[keep], c[keep]
    pts = []
    for i, j in zip(r, c):
        rr = slice(max(i - 2, 0), i + 3)
        cols = np.arange(j - 2, j + 3) % 1440
        box = msl[rr][:, cols]
        a, b = np.unravel_index(np.argmin(box), box.shape)
        pts.append((box[a, b], max(i - 2, 0) + a, cols[b]))
    pts.sort()
    sel = []
    for p, i, j in pts:
        v = XYZ[i, j]
        if all(np.arccos(np.clip(v @ XYZ[a, b], -1, 1)) * RE > 400 for _, a, b in sel):
            sel.append((p, i, j))
    return sel


def bearing(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    y = np.sin(lo2 - lo1) * np.cos(la2)
    x = np.cos(la1) * np.sin(la2) - np.sin(la1) * np.cos(la2) * np.cos(lo2 - lo1)
    return (np.degrees(np.arctan2(y, x)) + 360) % 360


COLS = ["time", "track", "fix_lat", "fix_lon", "fix_g800", "cen_lat", "cen_lon", "cen_dkm", "msl",
        "g800", "max_lat", "max_lon", "max_dkm", "max_brg", "max_dgl", "max_dic", "max_dland",
        "max_ice", "max_wdir", "max_wspd",
        "g800_gl100", "g800_gl300", "g800_land50", "g800_ice", "g800_r400"]


def one(args):
    stamp, fixes = args
    path = f"{OUT}/{stamp}.csv"
    if os.path.exists(path):
        return stamp
    t = np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}")
    ti = int((t - T0) / np.timedelta64(1, "h"))
    msl = get("mean_sea_level_pressure", ti) / 100.0
    gust = get("instantaneous_10m_wind_gust", ti)
    need_wind = any(f[1] >= WIND_LAT for f in fixes)
    if need_wind:
        ice = get("sea_ice_cover", ti)
        u = get("10m_u_component_of_wind", ti)
        v = get("10m_v_component_of_wind", ti)
    sel = lows(msl)
    cen = np.array([XYZ[i, j] for _, i, j in sel])
    tree = cKDTree(cen)
    gi = np.where(LSM)                                   # every ocean point, no 17 m/s floor
    d, k = tree.query(XYZ[gi], k=1)
    dkm = 2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE
    g = gust[gi]
    rows = []
    for trk, flat, flon, fg in fixes:
        fv = XYZ[int(round((90 - flat) / 0.25)) - R0, int(round(flon / 0.25)) % 1440]
        dd, kk = tree.query(fv, k=1)
        cdist = 2 * np.arcsin(min(dd / 2, 1)) * RE
        p, ci, cj = sel[kk]
        own = (k == kk) & (dkm <= 800)
        if not own.any():
            continue
        gg = np.where(own, g, -1.0)
        a = int(np.argmax(gg))
        mi, mj = gi[0][a], gi[1][a]

        def mx(extra):
            w = own & extra
            return round(float(g[w].max()) * KT, 1) if w.any() else 0.0

        if need_wind:
            ice_pts = ice[gi]
            noice = ~(ice_pts > 0.15)                    # NaN (open ocean in ERA5) counts as ice-free
            mice = float(ice[mi, mj])
            uu, vv = float(u[mi, mj]), float(v[mi, mj])
            wdir = (np.degrees(np.arctan2(-uu, -vv)) + 360) % 360      # direction wind blows FROM
            wspd = np.hypot(uu, vv) * KT
            gice = mx(noice)
        else:
            mice, wdir, wspd, gice = np.nan, np.nan, np.nan, np.nan
        rows.append((stamp, trk, flat, flon, fg, LAT[ci], LON[cj], round(cdist, 1), round(float(p), 1),
                     round(float(g[a]) * KT, 1), LAT[mi], LON[mj], round(float(dkm[a]), 1),
                     round(float(bearing(LAT[ci], LON[cj], LAT[mi], LON[mj])), 1),
                     round(float(D_GL[mi, mj]), 1), round(float(D_IC[mi, mj]), 1), round(float(D_LAND[mi, mj]), 1),
                     round(mice, 3) if mice == mice else "", round(float(wdir), 1) if wdir == wdir else "",
                     round(float(wspd), 1) if wspd == wspd else "",
                     mx(D_GL[gi] > 100), mx(D_GL[gi] > 300), mx(D_LAND[gi] > 50),
                     gice if gice == gice else "", mx(dkm <= 400)))
    tmp = path + ".tmp"
    pd.DataFrame(rows, columns=COLS).to_csv(tmp, index=False)
    os.replace(tmp, path)
    return stamp


def targets():
    cat = pd.read_csv(os.path.join(RES, "era5_hf_catalog.csv"))
    trk = pd.read_csv(os.path.join(RES, "era5_hf_catalog_tracks.csv"))
    ev = cat[(cat.role == "event") & (cat.basin == "atl") & (cat.season >= FIRST_SEASON)]
    h = trk[trk.track.isin(ev.track) & (trk.basin == "atl") & (trk.g800 >= THR)]
    jobs = {}
    for r in h.itertuples():
        jobs.setdefault(str(r.time), []).append((r.track, r.lat, r.lon, r.g800))
    return sorted(jobs.items(), reverse=True)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    jobs = targets()
    nw = sum(any(f[1] >= WIND_LAT for f in fx) for _, fx in jobs)
    print(f"{len(jobs)} times, {nw} with wind and ice; est {(len(jobs) * 5.37 + nw * 7.0) / 1000:.1f} GB", flush=True)
    t = time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 12) as p:
        for i, s in enumerate(p.imap_unordered(one, jobs)):
            if i % 50 == 0:
                print(f"{s} {i + 1}/{len(jobs)} {time.time() - t:.0f}s", flush=True)
