"""Stage 1: regional gust outcome G_T and the two pressure ingredients, at every candidate time.

ERA5 proxy (ARCO-ERA5, 0.25 degree); lows are pipeline A's (research/era5/hf_history).
Definitions are in PREREGISTRATION.md (committed before this was run). For each
candidate time (00/12 UTC, Nov-Mar, 2004-05 on, a pipeline A track fix at 55-67N,
50W-15W) it fetches gust and MSLP and records

  G_T      largest gust (kt) over ocean points within 300 km of Greenland, 58-72N,
           52W-10W, at least 400 km from every pipeline-A-detected low; 0 if none
  max_lat, max_lon, max_dlow   where that maximum is and its distance to the nearest low
  GH       mean MSLP over Greenland land points 64-78N, 55W-30W (hPa)
  GRAD     mean MSLP over ocean 0-200 km from Greenland minus ocean 400-600 km from it,
           62-70N, 40W-10W (hPa)

and stores a 50-80N, 80W-0 window of MSLP for later use. Output (ignored, resumable):
$ERA5_WORK/greenland_jets/stage1/<stamp>.csv and mslp/<stamp>.npy.

usage: stage1.py [NPROC]     then build_cases.py
"""
import os, sys, time, urllib.request
import numpy as np, pandas as pd, numcodecs
from scipy import ndimage
from scipy.spatial import cKDTree
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
OUT = os.path.join(WORK, "stage1")
FIX = os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz")

ARCO = "https://storage.googleapis.com/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
T0 = np.datetime64("1900-01-01T00")
R0, R1 = 60, 281                       # 75N .. 20N, as pipeline A's detection band
LAT = 90 - 0.25 * np.arange(R0, R1)
LON = 0.25 * np.arange(1440)
RE = 6371.0
KT = 1 / 0.514444
THR = 71.7
_codec = numcodecs.Blosc()

# window stored: lat 80N..50N (rows 40..160), lon 280E..360E (cols 1120..1440)
WROW = slice(40, 160)
WCOL = slice(1120, 1440)


def fetch(var, ti):
    url = f"{ARCO}/{var}/{ti}.0.0"
    for k in range(7):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                raw = r.read()
            return np.frombuffer(_codec.decode(raw), "<f4").reshape(721, 1440), len(raw)
        except Exception:
            if k == 6:
                raise
            time.sleep(2 ** k)


def static():
    path = os.path.join(WORK, "static.npz")
    if os.path.exists(path):
        return dict(np.load(path))
    lsm, _ = fetch_const()
    land = lsm >= 0.5
    lab, _ = ndimage.label(land)
    gl = lab == lab[int(round((90 - 72.0) / 0.25)), int(round(320.0 / 0.25))]
    assert gl.sum() > 1000
    la = np.radians(90 - 0.25 * np.arange(721))
    lo = np.radians(0.25 * np.arange(1440))
    LA2, LO2 = np.meshgrid(la, lo, indexing="ij")
    xyz = np.stack([np.cos(LA2) * np.cos(LO2), np.cos(LA2) * np.sin(LO2), np.sin(LA2)], -1)
    tree = cKDTree(xyz[gl])
    rows = slice(0, 721)
    d, _ = tree.query(xyz.reshape(-1, 3), k=1)
    dgl = (2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE).reshape(721, 1440)
    np.savez(path, lsm=lsm, gl=gl, dgl=dgl)
    return dict(lsm=lsm, gl=gl, dgl=dgl)


def fetch_const():
    # land-sea mask at 2000-01-01 00Z, as pipeline A
    return fetch("land_sea_mask", 876576)


S = None
LA_FULL = 90 - 0.25 * np.arange(721)
XYZ = None


def init():
    global S, XYZ, ELIG, GHM, NEAR, FAR
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(WORK, "mslp"), exist_ok=True)
    S = static()
    la, lo = np.meshgrid(np.radians(LA_FULL), np.radians(LON), indexing="ij")
    XYZ = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)
    LA2, LO2 = np.meshgrid(LA_FULL, LON, indexing="ij")
    ocean = S["lsm"] < 0.5
    inlon = (LO2 >= 308.0) & (LO2 <= 350.0)
    ELIG = ocean & (S["dgl"] <= 300.0) & (LA2 >= 58.0) & (LA2 <= 72.0) & inlon
    GHM = S["gl"] & (LA2 >= 64.0) & (LA2 <= 78.0) & (LO2 >= 305.0) & (LO2 <= 330.0)
    sect = ocean & (LA2 >= 62.0) & (LA2 <= 70.0) & (LO2 >= 320.0) & (LO2 <= 350.0)
    NEAR = sect & (S["dgl"] <= 200.0)
    FAR = sect & (S["dgl"] >= 400.0) & (S["dgl"] <= 600.0)


def lows(msl):
    """Pipeline A's detection (hf_history/extract.py): rows are the 75N..20N band."""
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
        v = XYZ[i + R0, j]
        if all(np.arccos(np.clip(v @ XYZ[a + R0, b], -1, 1)) * RE > 400 for _, a, b in sel):
            sel.append((p, i, j))
    return sel


def one(stamp):
    path = f"{OUT}/{stamp}.csv"
    if os.path.exists(path):
        return stamp
    t = np.datetime64(f"{stamp[:4]}-{stamp[4:6]}-{stamp[6:8]}T{stamp[8:10]}")
    ti = int((t - T0) / np.timedelta64(1, "h"))
    msl_full, n1 = fetch("mean_sea_level_pressure", ti)
    gust_full, n2 = fetch("instantaneous_10m_wind_gust", ti)
    msl = msl_full / 100.0
    sel = lows(msl[R0:R1])
    cen = np.array([XYZ[i + R0, j] for _, i, j in sel])
    tree = cKDTree(cen)
    ei, ej = np.where(ELIG)
    d, _ = tree.query(XYZ[ei, ej], k=1)
    dlow = 2 * np.arcsin(np.clip(d / 2, 0, 1)) * RE
    ok = dlow >= 400.0
    gt, mlat, mlon, mdl, mdgl = 0.0, np.nan, np.nan, np.nan, np.nan
    n_ok = int(ok.sum())
    if n_ok:
        g = gust_full[ei[ok], ej[ok]] * KT
        a = int(np.argmax(g))
        gt = float(g[a])
        mlat, mlon = LA_FULL[ei[ok][a]], LON[ej[ok][a]]
        mdl = float(dlow[ok][a])
    gh = float(msl[GHM].mean())
    grad = float(msl[NEAR].mean() - msl[FAR].mean())
    np.save(f"{WORK}/mslp/{stamp}.npy", msl_full[WROW, WCOL].astype("float32"))
    row = dict(time=stamp, G_T=round(gt, 1), max_lat=mlat, max_lon=mlon, max_dlow=round(mdl, 1) if mdl == mdl else np.nan,
               n_elig=int(ELIG.sum()), n_ok=n_ok, n_lows=len(sel), GH=round(gh, 2), GRAD=round(grad, 2),
               gmax_all_elig=round(float(gust_full[ELIG].max()) * KT, 1), bytes=n1 + n2)
    tmp = path + ".tmp"
    pd.DataFrame([row]).to_csv(tmp, index=False)
    os.replace(tmp, path)
    return stamp


def targets():
    f = pd.read_csv(FIX)
    a = f[f.basin == "atl"].copy()
    a["lon"] = np.where(a.lon > 180, a.lon - 360, a.lon)
    a["mon"] = (a.time // 10000) % 100
    r = a[a.lat.between(55, 67) & a.lon.between(-50, -15) & a.mon.isin([11, 12, 1, 2, 3])]
    return sorted(str(t) for t in r.time.unique())


if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True)
    init()
    tg = targets()
    print(f"{len(tg)} candidate times, est {len(tg) * 5.4 / 1000:.1f} GB", flush=True)
    if len(sys.argv) > 2 and sys.argv[2] == "test":
        for s in tg[:2]:
            print(one(s)); print(open(f"{OUT}/{s}.csv").read())
        sys.exit()
    t = time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8, initializer=init) as p:
        for i, s in enumerate(p.imap_unordered(one, tg)):
            if i % 100 == 0:
                print(f"{s} {i + 1}/{len(tg)} {time.time() - t:.0f}s", flush=True)
