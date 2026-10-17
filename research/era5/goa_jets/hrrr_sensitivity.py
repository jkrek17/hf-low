"""POST HOC sensitivity (PREREGISTRATION_hrrr.md, Deviations): HRRR zone-maximum gust restricted to ocean cells at least 10 km and at least 25 km from land, because the registered
zone maximum sat within 10 km of land for most HRRR events (results/hrrr_h2.csv). Compares with the ERA5 zone maximum (0.25 degree cells, about 25 km) at the same times.
Writes results/hrrr_sens.csv (per time and zone) and results/hrrr_sens_summary.csv, hrrr_sens.txt. usage: hrrr_sensitivity.py
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hrrr_analyse as A
from hrrr_pull import vtimes
RES = A.RES


def main():
    la, lo, land = A.load_grid(); dl, _ = A.dland_grid(la, lo, land); ocean = ~land
    box = lambda s, n, w, e: (la >= s) & (la <= n) & (lo >= w) & (lo <= e)
    ZM = {z: box(*b) & ocean for z, b in A.ZONES.items()}
    rows = []
    for t in vtimes():
        fn = os.path.join(A.WORK, f"{t:%Y%m%d%H}.npz")
        if not os.path.exists(fn): continue
        z = np.load(fn); g = z["gust"].astype(float); u = z["u"].astype(float); v = z["v"].astype(float)
        for zn, m in ZM.items():
            r = dict(time=t, zone=zn)
            for k in (0, 10, 25):
                mm = m & (dl >= k); r[f"gust_ge{k}km"] = float(g[mm].max()) if mm.any() else np.nan
            rows.append(r)
    S = pd.DataFrame(rows); S.to_csv(os.path.join(RES, "hrrr_sens.csv"), index=False)
    f = pd.read_csv(os.path.join(RES, "hrrr_features.csv"), parse_dates=["time"])
    s0 = pd.read_csv(os.path.join(RES, "stage0_zone_gust.csv"), parse_dates=["time"])
    e = pd.concat([s0[["time", f"{z}_gust_kt"]].set_axis(["time", "era5_gust"], axis=1).assign(zone=z) for z in A.ZONES])
    d = f[["time", "zone", "wdir"]].merge(S, on=["time", "zone"]).merge(e, on=["time", "zone"])
    d["regime"] = [A.regime(z, w) for z, w in zip(d.zone, d.wdir)]
    out = []
    for zn, rg in (("Z1", "BARRIER"), ("Z2", "BARRIER"), ("Z3", "GAP")):
        x = d[(d.zone == zn) & (d.regime == rg)]
        for k in (0, 10, 25):
            c = f"gust_ge{k}km"
            out.append(dict(zone=zn, regime=rg, min_km_from_land=k, n=len(x), median_ratio_to_era5=float((x[c] / x.era5_gust).median()), spearman=float(x[[c, "era5_gust"]].corr(method="spearman").iloc[0, 1]),
                            hrrr54=int((x[c] >= 54).sum()), era5_54=int((x.era5_gust >= 54).sum()), both54=int(((x[c] >= 54) & (x.era5_gust >= 54)).sum()),
                            pct_hrrr54=100 * float((x[c] >= 54).mean()), pct_era5_54=100 * float((x.era5_gust >= 54).mean()), pct_hrrr64=100 * float((x[c] >= 64).mean()), pct_era5_64=100 * float((x.era5_gust >= 64).mean())))
    O = pd.DataFrame(out); O.to_csv(os.path.join(RES, "hrrr_sens_summary.csv"), index=False)
    txt = "POST HOC sensitivity: HRRR zone-max gust restricted to ocean cells at least k km from land (barrier / gap regime from HRRR wind)\n" + O.round(3).to_string(index=False)
    open(os.path.join(RES, "hrrr_sens.txt"), "w").write(txt); print(txt)


if __name__ == "__main__":
    main()
