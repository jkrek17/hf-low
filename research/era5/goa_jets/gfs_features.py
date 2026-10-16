"""Follow-up stage 3 (PREREGISTRATION_followup.md): ingredients from the GFS FORECAST fields of gfs_pull.py, defined exactly as features.py defines them
for ERA5 (same boxes, land mask, regime rule). Forecast zone gust = largest GFS gust over the zone's ocean cells; its location gives LOW_DIST.
GFS gust is a different index from ERA5 gust and is not corrected. STAB is not computed (no SST pulled).
Writes results/gfs_features.csv, one row per (lead, valid time, zone). usage: gfs_features.py
"""
import os, sys, glob
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import features as F

GW = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "goa_jets", "gfs")


def run():
    rows = []; miss = []
    for fn in sorted(glob.glob(os.path.join(GW, "*_*.npz"))):
        if fn.endswith(".tmp.npz"): continue
        lead, st = os.path.basename(fn)[:-4].split("_")
        z = np.load(fn); p, u, v, g = z["msl"], z["u"], z["v"], z["gust"]
        t = pd.Timestamp(st[:4] + "-" + st[4:6] + "-" + st[6:8] + " " + st[8:10] + ":00")
        i = np.where(F.M_LOW, p, np.inf); k = np.unravel_index(np.argmin(i), i.shape)
        lowp, lowla, lowlo = p[k], F.LA[k], F.LO[k]
        interior_p = p[F.M_INT].mean(); exit_p = p[F.M_EXIT].mean()
        for zn in F.ZONES:
            m = F.ZM[zn]
            um, vm = u[m].mean(), v[m].mean()
            wdir = (np.degrees(np.arctan2(-um, -vm))) % 360
            gz = np.where(m, g, -1.0); kk = np.unravel_index(np.argmax(gz), gz.shape)
            glat, glon, gmax = F.LA[kk], F.LO[kk], float(gz[kk])
            pile = np.nan
            if zn in ("Z1", "Z2"):
                zc = F.ZC[zn]
                a = p[zc & (F.dland <= 100)]; b = p[zc & (F.dland >= 300) & (F.dland <= 400)]
                pile = a.mean() - b.mean() if len(a) and len(b) else np.nan
            rows.append(dict(lead=int(lead), time=t, zone=zn, wdir=wdir, wspd_kt=float(np.hypot(um, vm)), GFSGUST=gmax, glat=glat, glon=glon,
                             DEPTH=-float(lowp), low_lat=float(lowla), low_lon=float(lowlo), LOW_DIST=float(F.hav(lowla, lowlo, glat, glon)),
                             COASTPILE=pile, INTERIOR_DP=(interior_p - exit_p) if zn == "Z3" else np.nan))
    df = pd.DataFrame(rows).sort_values(["lead", "time", "zone"])
    df.to_csv(os.path.join(HERE, "results", "gfs_features.csv"), index=False)
    print(len(df), "rows;", df.groupby(["lead", "zone"]).size().to_dict())


if __name__ == "__main__":
    run()
