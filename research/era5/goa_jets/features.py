"""Stage 2: ingredients per pulled time and zone (ERA5 PROXY fields from stage1_pull.py). Definitions are in PREREGISTRATION.md.

Writes results/features.csv, one row per (time, zone). usage: features.py
"""
import os, sys
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "goa_jets", "fields")
CACHE = "/mnt/project-files/era5-cache"
KT = 1.943844
R = 6371.0
ROW0, NR, COL0, NC = 80, 105, 740, 201
lat = 90 - 0.25 * (ROW0 + np.arange(NR))                 # 70N .. 43.75N
lon = 0.25 * (COL0 + np.arange(NC)) - 360                # -175 .. -125
LA, LO = np.meshgrid(lat, lon, indexing="ij")
lsm = np.load(os.path.join(CACHE, "lsm_nh_0p25.npy"))[ROW0:ROW0 + NR, COL0:COL0 + NC]
ocean = lsm < 0.5
ZONES = {"Z1": (56, 61, -150, -138), "Z2": (53, 58, -140, -130), "Z3": (55, 60, -158, -148)}


def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]) * R


def hav(la0, lo0, la, lo):
    p1, p2 = np.radians(la0), np.radians(la); dl = np.radians(lo - lo0)
    h = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(h))


# distance to nearest land cell (km), chord approximation is fine at these scales
land = ~ocean
dland = np.full(ocean.shape, np.nan)
tree = cKDTree(xyz(LA[land], LO[land]))
dd, _ = tree.query(xyz(LA[ocean], LO[ocean]))
dland[ocean] = dd

box = lambda s, n, w, e: (LA >= s) & (LA <= n) & (LO >= w) & (LO <= e)
M_LOW = box(44, 62, -175, -125) & ocean
M_INT = box(62, 66, -155, -145)
M_EXIT = box(56, 59, -156, -150) & ocean
M_STAB = box(50, 58, -160, -145) & ocean
ZM = {z: box(*b) & ocean for z, b in ZONES.items()}
# coast pile-up: ocean cells within the zone's longitudes, 0-100 km and 300-400 km from land
ZC = {z: (box(44, 66, b[2], b[3]) & ocean) for z, b in ZONES.items()}
dy = 0.25 * 111.19


def grad_mag(m):
    gx = np.gradient(m, axis=1) / (dy * np.cos(np.radians(LA)))
    gy = np.gradient(m, axis=0) / dy
    return np.hypot(gx, gy) * 100.0       # hPa per 100 km


def run():
    pt = pd.read_csv(os.path.join(HERE, "results", "pull_times.csv"))
    s0 = pd.read_csv(os.path.join(HERE, "results", "stage0_zone_gust.csv")).set_index("time")
    rows = []
    for t in pt.time:
        st = t[:4] + t[5:7] + t[8:10] + t[11:13]
        z = np.load(os.path.join(WORK, st + ".npz"))
        u, v, p = z["u10"], z["v10"], z["msl"] / 100.0
        t2, sst = z["t2m"], z["sst"]
        i = np.where(M_LOW, p, np.inf); k = np.unravel_index(np.argmin(i), i.shape)
        lowp, lowla, lowlo = p[k], LA[k], LO[k]
        gm = grad_mag(p)
        interior_p = p[M_INT].mean(); exit_p = p[M_EXIT].mean()
        interior_t = t2[M_INT].mean() - 273.15
        stab = np.nanmean(sst[M_STAB] - t2[M_STAB])
        for zn in ZONES:
            m = ZM[zn]
            um, vm = u[m].mean(), v[m].mean()
            wdir = (np.degrees(np.arctan2(-um, -vm))) % 360          # direction wind blows FROM
            spd = float(np.hypot(um, vm) * KT)
            g, glat, glon = s0.loc[t, [f"{zn}_gust_kt", f"{zn}_lat", f"{zn}_lon"]]
            dist = hav(glat, glon, LA, LO)
            near = ocean & (dist <= 200)
            near_all = near.any()
            pile = np.nan
            if zn in ("Z1", "Z2"):
                zc = ZC[zn]
                a = p[zc & (dland <= 100)]; b = p[zc & (dland >= 300) & (dland <= 400)]
                pile = a.mean() - b.mean() if len(a) and len(b) else np.nan
            rows.append(dict(time=t, zone=zn, gust_kt=g, glat=glat, glon=glon, wdir=wdir, wspd_kt=spd,
                             DEPTH=-lowp, low_lat=lowla, low_lon=lowlo,
                             LOW_DIST=hav(lowla, lowlo, glat, glon),
                             GRAD=gm[near].mean() if near_all else np.nan,
                             COASTPILE=pile,
                             INTERIOR_DP=(interior_p - exit_p) if zn == "Z3" else np.nan,
                             INTERIOR_T2M=interior_t if zn == "Z3" else np.nan,
                             STAB=stab))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "results", "features.csv"), index=False)
    print(len(df), "rows"); print(df.describe().round(2).T.to_string())


if __name__ == "__main__":
    run()
