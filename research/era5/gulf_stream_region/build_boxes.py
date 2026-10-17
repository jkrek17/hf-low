"""Collocate ASCAT 1-degree box passes with ERA5 controls (table A, winters 2013-14..2024-25) and GFS 24 h forecasts (table B, 2021-22..2025-26).
Usage: python3 build_boxes.py SHARED_DIR OUT_DIR   (SHARED_DIR holds ascat/, msl/, era5temp/, gfs/ ; /mnt/project-files/gulf-stream)
Vg: ERA5 MSLP proxy gradient wind; dT: ERA5 SST - 2 m T (A) or GFS skin - 2 m T (B)."""
import sys, os, glob, numpy as np, pandas as pd
from scipy.ndimage import uniform_filter
SH, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
KT = 1.94384449
LSM = np.load("/mnt/project-files/era5-cache/lsm_nh_0p25.npy")          # 361 x 1440, 90N..0N, 0..359.75E
LATB, LONB = np.arange(29, 41), np.arange(280, 305)                       # box SW corners (deg N, deg E): 29-41N, 80-55W
def land_box(la, lo):
    r0, r1 = int((90 - (la + 1 + 0.5)) / 0.25), int((90 - (la - 0.5)) / 0.25) + 1
    c0, c1 = int((lo - 0.5) / 0.25), int((lo + 1 + 0.5) / 0.25) + 1
    return (LSM[r0:r1, c0:c1] >= 0.5).any()
LAND = np.array([[land_box(la, lo) for lo in LONB] for la in LATB])      # (12, 25)
def dec_t(x): return x.astype(np.float64)
def boxes_for_file(fn):
    z = np.load(fn); data = z["data"]; dates = z["dates"]; lat, lon = z["lat"], z["lon"]
    ri = np.nonzero((lat >= 29) & (lat < 41))[0]; ci = np.nonzero((lon >= 280) & (lon < 305))[0]
    d = data[:, :, :, ri][:, :, :, :, ci]                                   # (day, seg, var, 48, 100)
    nd = d.shape[0]
    t, w, di, rn = [d[:, :, k].astype(np.float64) for k in range(4)]
    ok = (t <= 250) & (w <= 250) & ((d[:, :, 3] & 0b11111101) == 0)
    sh = (nd, 2, 12, 4, 25, 4)
    def bx(a): return a.reshape(sh)
    okb = bx(ok); n = okb.sum(axis=(3, 5))
    wb = np.where(okb, bx(w), 0).sum(axis=(3, 5)) / np.maximum(n, 1) * 0.2 * KT
    tb = np.where(okb, bx(t), 0).sum(axis=(3, 5)) / np.maximum(n, 1)
    ang = np.radians(bx(di) * 1.5 + 180.0)          # RSS direction byte is oceanographic (going-to); +180 gives the meteorological from-direction (checked on buoys)
    sx = np.where(okb, np.sin(ang), 0).sum(axis=(3, 5)); sy = np.where(okb, np.cos(ang), 0).sum(axis=(3, 5))
    dirb = (np.degrees(np.arctan2(sx, sy)) % 360)
    keep = (n >= 8) & ~LAND[None, None]
    ii = np.nonzero(keep)
    day0 = pd.to_datetime(dates.astype(str), format="%Y%m%d").values.astype("datetime64[m]").astype(np.int64)
    tmin = day0[ii[0]] + np.rint(tb[ii] * 6).astype(np.int64)                  # minutes since 1970
    return pd.DataFrame(dict(day=dates[ii[0]].astype(int), seg=ii[1], lat=LATB[ii[2]] + 0.5, lon=LONB[ii[3]] + 0.5,
                             tmin=tmin, v_ascat=wb[ii], dir_ascat=dirb[ii], ncell=n[ii]))
def load_ascat():
    fr = []
    for fn in sorted(glob.glob(f"{SH}/ascat/ascat_*.npz")):
        sen, w = os.path.basename(fn)[6:-4].split("_"); df = boxes_for_file(fn); df["sat"] = sen; df["winter"] = int(w); fr.append(df)
    return pd.concat(fr, ignore_index=True)
def nearest_12h(tmin):
    h = tmin / 60.0; k = np.rint(h / 12.0) * 12.0
    return k, np.abs(h - k)
def tkey(hours):                                                           # hours since 1970 -> YYYYMMDDHH int
    return pd.to_datetime(np.asarray(hours).astype("int64") * 3600, unit="s").strftime("%Y%m%d%H").astype(np.int64)
def vg_field(p):
    """gradient wind (kt) on the MSLP crop grid, central difference +-0.5 deg (2 cells); NaN at the edges."""
    la = 55 - 0.25 * np.arange(p.shape[0])
    dy = 1.0 * 111_000.0
    dpdy = np.full(p.shape, np.nan); dpdx = np.full(p.shape, np.nan)
    dpdy[2:-2] = (p[:-4] - p[4:]) / dy
    dx = (dy * np.cos(np.radians(la)))[:, None]
    dpdx[:, 2:-2] = (p[:, 4:] - p[:, :-4]) / dx
    f = (2 * 7.292e-5 * np.sin(np.radians(la)))[:, None]
    return np.hypot(dpdx, dpdy) / (1.25 * f) * KT
def front_lat(sst):
    """latitude of the north wall per box-centre longitude: strongest southward SST increase (warm to the south) between 32N and 42N,
    1 degree box-smoothed SST, at least 0.01 K/km; sst int16 crop (81,141), 45N..25N, 275E..310E."""
    s = sst.astype(np.float64) / 100.0 + 200.0; m = (s > 250).astype(float)
    num = uniform_filter(np.where(m > 0, s, 0), 4, mode="nearest"); den = uniform_filter(m, 4, mode="nearest")
    S = np.where(den > 0.5, num / np.maximum(den, 1e-9), np.nan)
    g = np.full(S.shape, np.nan); g[1:-1] = (S[2:] - S[:-2]) / 55.6              # rows run north to south: south minus north, K/km
    lat = 45 - 0.25 * np.arange(81); out = np.full(25, np.nan)
    for k, lo in enumerate(LONB):
        col = g[:, int((lo + 0.5 - 275) / 0.25)]
        sel = (lat >= 32) & (lat <= 42) & np.isfinite(col)
        if sel.any():
            j = np.nonzero(sel)[0][np.nanargmax(col[sel])]
            if col[j] >= 0.01 and 32.0 < lat[j] < 42.0: out[k] = lat[j]      # a maximum on the search edge means the front lies outside the window: undefined
    return out
def temp_index(rowlat, collon, lat, lon):                                       # ERA5 temp / GFS crop: 45N..25N, 275E..310E
    return np.rint((45 - lat) / 0.25).astype(int), np.rint((lon - 275) / 0.25).astype(int)
def build_A(asc):
    asc = asc[asc.winter <= 2024].copy()
    hh, dt_h = nearest_12h(asc.tmin.values); asc["ctl_h"] = hh; asc["dt_h"] = dt_h; asc = asc[asc.dt_h <= 4].copy()
    asc["key"] = tkey(asc.ctl_h.values); asc["datekey"] = (asc.key // 100).astype(int)
    vg = {}; t2 = {}; sst = {}; front = {}
    for fn in sorted(glob.glob(f"{SH}/msl/msl_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["t"], z["msl"]): vg[int(t)] = vg_field(a.astype(np.float32) + 85000)
    for fn in sorted(glob.glob(f"{SH}/era5temp/era5_temp_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["t2_time"], z["t2"]): t2[int(t)] = a
        for t, a in zip(z["sst_time"], z["sst"]): sst[int(t) // 100] = a; front[int(t) // 100] = front_lat(a)
    ri = np.rint((55 - asc.lat.values) / 0.25).astype(int); ci = np.rint((asc.lon.values - 270) / 0.25).astype(int)
    ti, tj = temp_index(None, None, asc.lat.values, asc.lon.values)
    out = []
    for key, g in asc.groupby("key"):
        if key not in vg or key not in t2 or (key // 100) not in sst: continue
        ix = g.index.values; r = ri[asc.index.get_indexer(ix)]; c = ci[asc.index.get_indexer(ix)]; a_ = ti[asc.index.get_indexer(ix)]; b_ = tj[asc.index.get_indexer(ix)]
        v = vg[key][r, c]; tt = t2[key][a_, b_].astype(float) / 100 + 200
        s_ = sst[key // 100][a_, b_].astype(float) / 100 + 200; s_[s_ < 250] = np.nan
        fr = front[key // 100][((g.lon.values - 280.5)).astype(int)]
        out.append(pd.DataFrame(dict(idx=ix, vg=v, t2m=tt, sst=s_, front=fr)))
    c = pd.concat(out).set_index("idx"); asc = asc.join(c, how="inner")
    asc["dT"] = (asc.sst - asc.t2m).clip(-3, 12); asc["dist"] = asc.lat - asc.front
    asc["band"] = np.where(asc.dist.abs() <= 1, "in", np.where(asc.dist.abs() >= 2, "out", "mid")); asc.loc[asc.front.isna(), "band"] = "undef"
    asc["month"] = pd.to_datetime(asc.day.astype(str), format="%Y%m%d").dt.month
    return asc.drop(columns=["key", "datekey", "ctl_h"])
def build_B(asc):
    asc = asc[asc.winter >= 2021].copy()
    hh, dt_h = nearest_12h(asc.tmin.values); asc["dt_h"] = dt_h; asc = asc[asc.dt_h <= 4].copy(); asc["key"] = tkey(hh[asc.index.values] if False else nearest_12h(asc.tmin.values)[0])
    gf = {}
    for fn in sorted(glob.glob(f"{SH}/gfs/gfs_f024_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["valid"], z["data"]): gf[int(t)] = a
    sst = {}; front = {}
    for fn in sorted(glob.glob(f"{SH}/era5temp/era5_temp_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["sst_time"], z["sst"]): sst[int(t) // 100] = a; front[int(t) // 100] = front_lat(a)
    ti, tj = temp_index(None, None, asc.lat.values, asc.lon.values); pos = {ix: k for k, ix in enumerate(asc.index.values)}
    out = []
    for key, g in asc.groupby("key"):
        if key not in gf or (key // 100) not in front: continue
        k = np.array([pos[i] for i in g.index.values]); a = gf[key]
        u = a[0][ti[k], tj[k]] / 100.0; vv = a[1][ti[k], tj[k]] / 100.0
        t2 = a[2][ti[k], tj[k]] / 100.0 + 200; sk = a[3][ti[k], tj[k]] / 100.0 + 200
        fr = front[key // 100][((g.lon.values - 280.5)).astype(int)]
        out.append(pd.DataFrame(dict(idx=g.index.values, v_gfs=np.hypot(u, vv) * KT, dir_gfs=(np.degrees(np.arctan2(-u, -vv)) % 360), t2m=t2, skt=sk, front=fr)))
    c = pd.concat(out).set_index("idx"); asc = asc.join(c, how="inner")
    asc["dT"] = (asc.skt - asc.t2m).clip(-3, 12); asc["dist"] = asc.lat - asc.front
    asc["band"] = np.where(asc.dist.abs() <= 1, "in", np.where(asc.dist.abs() >= 2, "out", "mid")); asc.loc[asc.front.isna(), "band"] = "undef"
    asc["month"] = pd.to_datetime(asc.day.astype(str), format="%Y%m%d").dt.month
    return asc.drop(columns=["key"])
if __name__ == "__main__":
    asc = load_ascat(); print("ASCAT box-passes (land-masked, >=8 cells):", len(asc), asc.groupby("sat").size().to_dict())
    A = build_A(asc); A.to_csv(f"{OUT}/boxes_A.csv.gz", index=False); print("A", len(A), A.winter.nunique(), "winters")
    B = build_B(asc); B.to_csv(f"{OUT}/boxes_B.csv.gz", index=False); print("B", len(B), B.winter.nunique(), "winters")
