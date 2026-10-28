"""Match ERA5 MSLP-derived gradient to ERA5 gust at ocean cells, western Atlantic, October-April 2004-05 on.
Usage: python3 build_pairs.py OUT.npz     (reads the shared folder caches; no download)
Plan: PREREGISTRATION.md. ERA5 proxy throughout."""
import glob, os, sys, datetime as dt
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ascat_check"))
from blosc_nolib import decode

SH = "/mnt/project-files"
KT = 1.943844
NM = 1852.0
RHO = 1.25
ROW0, COL0 = 55.0, 270.0                      # MSLP crop: 55N down to 20N, 270E to 310E, 0.25 degree
LATS = np.arange(30.0, 45.01, 2.0)            # every 8th cell
LONS = np.arange(282.0, 304.01, 2.0)


def key_to_h1900(k):
    d = dt.datetime(k // 1000000, k // 10000 % 100, k // 100 % 100, k % 100)
    return int((d - dt.datetime(1900, 1, 1)).total_seconds() // 3600)


def box(a, n):
    """mean over an n x n box (n odd), edges left as NaN."""
    c = np.cumsum(np.cumsum(np.pad(a, ((1, 0), (1, 0))), 0), 1)
    s = c[n:, n:] - c[:-n, n:] - c[n:, :-n] + c[:-n, :-n]
    out = np.full(a.shape, np.nan, np.float64); h = n // 2
    out[h:-h, h:-h] = s / (n * n)
    return out


def derivs(p, la):
    """px, py (Pa/m), pxx, pyy, pxy (Pa/m2) by central differences over 2 cells either side; NaN within 2 cells of the edge."""
    dy = 111_000.0; dxr = (dy * np.cos(np.radians(la)))[:, None]
    h = 2; nan = lambda: np.full(p.shape, np.nan)
    px, py, pxx, pyy, pxy = nan(), nan(), nan(), nan(), nan()
    py[h:-h] = (p[:-2 * h] - p[2 * h:]) / (2 * h * 0.25 * dy)
    px[:, h:-h] = (p[:, 2 * h:] - p[:, :-2 * h]) / (2 * h * 0.25 * dxr)
    pyy[h:-h] = (p[:-2 * h] - 2 * p[h:-h] + p[2 * h:]) / (h * 0.25 * dy) ** 2
    pxx[:, h:-h] = (p[:, 2 * h:] - 2 * p[:, h:-h] + p[:, :-2 * h]) / (h * 0.25 * dxr) ** 2
    pxy[h:-h] = (px[:-2 * h] - px[2 * h:]) / (2 * h * 0.25 * dy)
    return px, py, pxx, pyy, pxy


def slide_min(a, w):
    from numpy.lib.stride_tricks import sliding_window_view as sw
    h = w // 2; b = np.pad(a, h, constant_values=np.inf)
    b = sw(b, w, axis=0).min(-1); return sw(b, w, axis=1).min(-1)


def fields(p):
    """p: MSLP Pa (141,161). Returns gradient (hPa per 100 nm), G (kt), gradient wind V (kt, NaN if no real root), core mask."""
    la = ROW0 - 0.25 * np.arange(p.shape[0]); lo = COL0 + 0.25 * np.arange(p.shape[1])
    f = (2 * 7.292e-5 * np.sin(np.radians(la)))[:, None]
    px, py, *_ = derivs(p, la); g = np.hypot(px, py)
    G = g / (RHO * f) * KT
    ps = box(p, 5)                                                 # 1 degree smoothing for the isobar curvature
    sx, sy, sxx, syy, sxy = derivs(ps, la); gs = np.hypot(sx, sy)
    kap = (sxx * sy ** 2 - 2 * sxy * sx * sy + syy * sx ** 2) / np.maximum(gs, 1e-9) ** 3      # cyclonic positive, 1/m
    Gms = g / (RHO * f)
    disc = f ** 2 + 4 * kap * f * Gms
    V = np.where(disc >= 0, 2 * f * Gms / (f + np.sqrt(np.maximum(disc, 0))), np.nan) * KT
    # core: within 216 nm of a local minimum (+-3 degree box) below 1000 hPa, at least 3 degrees inside the crop edge
    mn = slide_min(p, 25); cand = (p == mn) & (p < 100000.0)
    cand[:12] = False; cand[-12:] = False; cand[:, :12] = False; cand[:, -12:] = False
    core = np.zeros(p.shape, bool)
    for r, c in zip(*np.nonzero(cand)):
        dy = (la - la[r]) * 111.0; dx = (lo - lo[c])[None, :] * 111.0 * np.cos(np.radians(la[r]))
        core |= np.hypot(dy[:, None], dx) <= 216 * 1.852
    return g * NM, G, V, core                                   # Pa/m x 1852 m/nm = Pa per nm; x 100 nm / 100 Pa = hPa per 100 nm


def main(out):
    lsm = np.load(f"{SH}/era5-cache/lsm_nh_0p25.npy")
    ri = np.rint((90 - LATS) / 0.25).astype(int); ci = np.rint(LONS / 0.25).astype(int)
    ok = np.ones((len(LATS), len(LONS)), bool)
    for a, r in enumerate(ri):
        for b, c in enumerate(ci):
            ok[a, b] = bool((lsm[r - 4:r + 5, c - 4:c + 5] < 0.5).all())      # ocean, no land within 1 degree
    mr = np.rint((ROW0 - LATS) / 0.25).astype(int); mc = np.rint((LONS - COL0) / 0.25).astype(int)
    gfiles = sorted(glob.glob(f"{SH}/era5-cache/gust_nh_12h/gust_*.npz"))
    hmap = {}
    for fi, fn in enumerate(gfiles):
        for j, h in enumerate(np.load(fn)["hours_since_1900"]): hmap[int(h)] = (fi, j)
    t2 = {}; sst = {}
    for fn in sorted(glob.glob(f"{SH}/gulf-stream/era5temp/era5_temp_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["t2_time"], z["t2"]): t2[int(t)] = a
        for t, a in zip(z["sst_time"], z["sst"]): sst[int(t) // 100] = a
    tr = np.rint((45 - LATS) / 0.25).astype(int); tc = np.rint((LONS - 275) / 0.25).astype(int)
    cols = {k: [] for k in "winter key lat lon grad G V gust core dT".split()}
    gcache = {}
    for fn in sorted(glob.glob(f"{SH}/gulf-stream/msl/msl_*.npz")):
        z = np.load(fn); win = int(os.path.basename(fn)[4:8]); print(win, flush=True)
        for t, a in zip(z["t"], z["msl"]):
            h = key_to_h1900(int(t))
            if h not in hmap: continue
            fi, j = hmap[h]
            if fi not in gcache:
                gcache.clear(); w = np.load(gfiles[fi])
                gcache[fi] = decode(w["data"].tobytes()).reshape(w["shape"])
            gu = gcache[fi][j, 140:281, 1080:1241].astype(np.float32) / 2.0       # lat 55..20, lon 270..310 (rows (90-lat)/.25)
            p = a.astype(np.float64) + 85000.0
            grad, G, V, core = fields(p)
            R, C = np.meshgrid(mr, mc, indexing="ij")
            sel = ok & np.isfinite(grad[R, C]) & np.isfinite(G[R, C])
            tt = t2.get(int(t)); ss = sst.get(int(t) // 100)
            if tt is not None and ss is not None:
                TA = tt[tr][:, tc].astype(float) / 100 + 200; SS = ss[tr][:, tc].astype(float) / 100 + 200
                dT = np.where(SS > 250, SS - TA, np.nan)
            else: dT = np.full(sel.shape, np.nan)
            LA, LO = np.meshgrid(LATS, LONS, indexing="ij")
            for k, v in dict(winter=np.full(sel.shape, win), key=np.full(sel.shape, int(t)), lat=LA, lon=LO, grad=grad[R, C], G=G[R, C],
                             V=V[R, C], gust=gu[R, C], core=core[R, C], dT=dT).items():
                cols[k].append(v[sel])
    np.savez_compressed(out, **{k: np.concatenate(v) for k, v in cols.items()})
    print("rows", sum(len(x) for x in cols["gust"]))


if __name__ == "__main__":
    main(sys.argv[1])
