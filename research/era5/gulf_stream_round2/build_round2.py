"""Add round 2 columns to the round 1 boxes tables.
A: d2m (ERA5 2 m dewpoint, K), tadv (thermal advection by the pressure-gradient wind, K/day), curv (geostrophic relative vorticity, 1e-5 /s), sst_mur (box-mean MUR SST, K),
   front_mur, dist_mur, band_mur (MUR front, 1 degree longitude strip), cur_u, cur_v (NOAA geostrophic current at the box centre, m/s; winters 2017-18 on).
B: dpt, pres (hPa), uflx, vflx (GFS 24 h, 18-24 h average stress), sst_mur, front_mur, dist_mur, band_mur.
Usage: python3 build_round2.py SHARED_DIR TABLES_IN OUT_DIR   (SHARED_DIR /mnt/project-files/gulf-stream; TABLES_IN holds boxes_A.csv.gz, boxes_B.csv.gz)"""
import sys, os, glob, numpy as np, pandas as pd
from scipy.ndimage import uniform_filter
SH, TIN, OUT = sys.argv[1], sys.argv[2], sys.argv[3]; os.makedirs(OUT, exist_ok=True)
RHO, OM = 1.25, 7.292e-5
def nearest_12h(tmin):
    h = tmin / 60.0; k = np.rint(h / 12.0) * 12.0; return k
def tkey(hours): return pd.to_datetime(np.asarray(hours).astype("int64") * 3600, unit="s").strftime("%Y%m%d%H").astype(np.int64)
# ---- MUR
def mur_day(sst):
    """sst int16 (401,701), lat 25..45 step .05, lon -85..-50 step .05 -> (front per box longitude [25], box-mean SST [12,25])"""
    s = np.where(sst > -32768, sst / 100.0 + 200.0, np.nan); ok = np.isfinite(s).astype(float)
    num = uniform_filter(np.nan_to_num(s), 11, mode="nearest"); den = uniform_filter(ok, 11, mode="nearest")
    S = np.where(den > 0.5, num / np.maximum(den, 1e-9), np.nan)
    g = np.full(S.shape, np.nan); g[5:-5] = (S[:-10] - S[10:]) / 55.6            # south minus north over +-0.25 deg of latitude... rows run south to north: (S[j-5]-S[j+5]) > 0 when warm to the south
    lat = 25 + 0.05 * np.arange(401); front = np.full(25, np.nan); box = np.full((12, 25), np.nan)
    for k in range(25):
        lo = 280.5 + k - 360; ci = int(round((lo + 85) / 0.05))
        col = np.nanmean(g[:, ci - 10:ci + 11], axis=1) if np.isfinite(g[:, ci - 10:ci + 11]).any() else np.full(401, np.nan)
        sel = (lat >= 32) & (lat <= 42) & np.isfinite(col)
        if sel.any():
            j = np.nonzero(sel)[0][np.nanargmax(col[sel])]
            if col[j] >= 0.01 and 32.0 < lat[j] < 42.0: front[k] = lat[j]
        for m in range(12):
            cj = int(round((29.5 + m - 25) / 0.05)); box[m, k] = np.nanmean(s[cj - 10:cj + 11, ci - 10:ci + 11]) if np.isfinite(s[cj - 10:cj + 11, ci - 10:ci + 11]).any() else np.nan
    return front, box
def load_mur(days):
    fr, bx = {}, {}
    for fn in sorted(glob.glob(f"{SH}/mur/mur_*.npz")):
        z = np.load(fn)
        for d, a in zip(z["day"], z["sst"]):
            if int(d) in days: fr[int(d)], bx[int(d)] = mur_day(a)
    return fr, bx
def add_mur(df, fr, bx):
    lat_i = (df.lat.values - 29.5).astype(int); lon_i = (df.lon.values - 280.5).astype(int)
    f = np.full(len(df), np.nan); b = np.full(len(df), np.nan)
    for d, idx in df.groupby("day").indices.items():
        if int(d) in fr: f[idx] = fr[int(d)][lon_i[idx]]; b[idx] = bx[int(d)][lat_i[idx], lon_i[idx]]
    df["sst_mur"] = b; df["front_mur"] = f; df["dist_mur"] = df.lat - f
    df["band_mur"] = np.where(df.dist_mur.abs() <= 1, "in", np.where(df.dist_mur.abs() >= 2, "out", "mid")); df.loc[df.front_mur.isna(), "band_mur"] = "undef"
    return df
# ---- A
def add_A(A):
    A = A.reset_index(drop=True); key = tkey(nearest_12h(A.tmin.values)); A["_key"] = key
    msl, t2, d2 = {}, {}, {}
    for fn in sorted(glob.glob(f"{SH}/msl/msl_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["t"], z["msl"]): msl[int(t)] = a
    for fn in sorted(glob.glob(f"{SH}/era5temp/era5_temp_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["t2_time"], z["t2"]): t2[int(t)] = a
    for fn in sorted(glob.glob(f"{SH}/era5d2m/era5_d2m_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["d2m_time"], z["d2m"]): d2[int(t)] = a
    ri = np.rint((55 - A.lat.values) / 0.25).astype(int); ci = np.rint((A.lon.values - 270) / 0.25).astype(int)
    ti = np.rint((45 - A.lat.values) / 0.25).astype(int); tj = np.rint((A.lon.values - 275) / 0.25).astype(int)
    for c in ("d2m", "tadv", "curv", "msl_hpa"): A[c] = np.nan
    dy = 55_500.0; f = 2 * OM * np.sin(np.radians(A.lat.values)); dx = dy * np.cos(np.radians(A.lat.values))
    for k, idx in A.groupby("_key").indices.items():
        if k not in msl or k not in t2 or k not in d2: continue
        p = msl[k].astype(np.float64) + 85000; T = t2[k].astype(np.float64) / 100 + 200; r, c, a, b = ri[idx], ci[idx], ti[idx], tj[idx]
        dpdy = (p[r - 2, c] - p[r + 2, c]) / dy; dpdx = (p[r, c + 2] - p[r, c - 2]) / dx[idx]
        dTdy = (T[a - 2, b] - T[a + 2, b]) / dy; dTdx = (T[a, b + 2] - T[a, b - 2]) / dx[idx]
        ug = -dpdy / (RHO * f[idx]); vg = dpdx / (RHO * f[idx])
        A.loc[idx, "tadv"] = -(ug * dTdx + vg * dTdy) * 86400.0
        lap = (p[r, c + 2] + p[r, c - 2] - 2 * p[r, c]) / dx[idx] ** 2 + (p[r - 2, c] + p[r + 2, c] - 2 * p[r, c]) / dy ** 2
        A.loc[idx, "curv"] = lap / (RHO * f[idx]) * 1e5
        A.loc[idx, "d2m"] = d2[k][a, b].astype(np.float64) / 100 + 200; A.loc[idx, "msl_hpa"] = p[r, c] / 100.0
    fr, bx = load_mur(set(A.day.unique())); A = add_mur(A, fr, bx)
    A["cur_u"] = np.nan; A["cur_v"] = np.nan; lat0, lon0 = 25.125, -84.875
    cur = {}
    for fn in sorted(glob.glob(f"{SH}/currents/cur_*.npz")):
        z = np.load(fn)
        for d, a in zip(z["day"], z["uv"]): cur[int(d)] = np.where(a > -32768, a / 1000.0, np.nan)
    fi = (A.lat.values - lat0) / 0.25; fj = (A.lon.values - 360 - lon0) / 0.25; i0 = np.floor(fi).astype(int); j0 = np.floor(fj).astype(int); wi = fi - i0; wj = fj - j0
    for d, idx in A.groupby("day").indices.items():
        if int(d) not in cur: continue
        a = cur[int(d)]
        for n, nm in enumerate(("cur_u", "cur_v")):
            g = a[n]; I, J, WI, WJ = i0[idx], j0[idx], wi[idx], wj[idx]
            A.loc[idx, nm] = (g[I, J] * (1 - WI) * (1 - WJ) + g[I + 1, J] * WI * (1 - WJ) + g[I, J + 1] * (1 - WI) * WJ + g[I + 1, J + 1] * WI * WJ)
    return A.drop(columns=["_key"])
# ---- B
def add_B(B):
    B = B.reset_index(drop=True); key = tkey(nearest_12h(B.tmin.values)); B["_key"] = key; gf = {}
    for fn in sorted(glob.glob(f"{SH}/gfs2/gfs2_f024_*.npz")):
        z = np.load(fn)
        for t, a in zip(z["valid"], z["data"]): gf[int(t)] = a
    ti = np.rint((45 - B.lat.values) / 0.25).astype(int); tj = np.rint((B.lon.values - 275) / 0.25).astype(int)
    for c in ("dpt", "pres", "uflx", "vflx"): B[c] = np.nan
    for k, idx in B.groupby("_key").indices.items():
        if k not in gf: continue
        a = gf[k]; r, c = ti[idx], tj[idx]
        B.loc[idx, "dpt"] = a[0][r, c] / 100.0 + 200; B.loc[idx, "pres"] = (a[1][r, c] / 0.1 + 80000) / 100.0
        B.loc[idx, "uflx"] = a[2][r, c] / 1000.0; B.loc[idx, "vflx"] = a[3][r, c] / 1000.0
    fr, bx = load_mur(set(B.day.unique())); B = add_mur(B, fr, bx)
    return B.drop(columns=["_key"])
if __name__ == "__main__":
    A = add_A(pd.read_csv(f"{TIN}/boxes_A.csv.gz")); A.to_csv(f"{OUT}/boxes2_A.csv.gz", index=False); print("A", len(A), A[["d2m", "tadv", "curv", "sst_mur", "front_mur", "cur_u"]].notna().mean().round(3).to_dict(), flush=True)
    B = add_B(pd.read_csv(f"{TIN}/boxes_B.csv.gz")); B.to_csv(f"{OUT}/boxes2_B.csv.gz", index=False); print("B", len(B), B[["dpt", "pres", "uflx", "sst_mur", "front_mur"]].notna().mean().round(3).to_dict(), flush=True)
