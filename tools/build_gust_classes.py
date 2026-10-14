#!/usr/bin/env python3
"""Build docs/data/gust-classes.js: monthly ERA5 gale / storm / hurricane-force gust hours for the map.

Derived data only. ERA5 PROXY (instantaneous 10 m gust at 00 and 12 UTC, 1 June 2004 to 31 May 2026, seasons 2004-05 to 2025-26,
16,070 times), not direct observation. ERA5 reads low in extreme storms and the samples are 12-hourly, so every number is
a lower bound on true gusts. The classes are NOMINAL PROXY classes on the gust index, not calibrated on observed events:
gale at or above 38 kt, storm at or above 54 kt, hurricane force (HF) at or above 71.7 kt (the HF cutoff used everywhere else).

    python3 tools/build_gust_classes.py --counts <gust_monthly_counts.npz> --lsm <lsm_nh_0p25.npy> --source <commit>

Inputs: research/era5/gust_class_maps (extract_cache.py, monthly_counts.py on the research branch; the cache is in the
project's shared folder, era5-cache/). The 0.25 degree grid is reduced to 1 degree blocks (mean of the ocean cells in a block)
over lat 20-80N and the two basin windows the rest of the map uses (Atlantic 100W-30E, Pacific 120E-110W).

Wire format (append-only):
  HF_GUST_CLASSES.win[basin] = {lat0, lon0, dlat, dlon, nlat, nlon, f: {g38|g54|g717: [12 arrays, January first]}}
  each array is row-major from (lat0, lon0), south to north then west to east; value = hours per season in that calendar
  month, in units of 0.1 h (sampled times at or above the cutoff x 12 h / 22 seasons); -1 = no ocean.
  Hours per season add across months.
"""
import argparse, json, os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WINDOWS = {"atl": (-100, 30), "pac": (120, 250)}
LAT0, LAT1 = 20, 80
CLASSES = (("g38", 38.0), ("g54", 54.0), ("g717", 71.7))
SEASONS = 22


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--lsm", required=True)
    ap.add_argument("--source", required=True)
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "data", "gust-classes.js"))
    ap.add_argument("--seasons", type=int, default=SEASONS, help="seasons in the counts (22 for 2004-05 on, 25 for 1979-80 to 2003-04)")
    ap.add_argument("--var", default="HF_GUST_CLASSES", help="window variable name (HF_GUST_CLASSES_1979 for the 1979-2003 set)")
    ap.add_argument("--period", default="1 June 2004 to 31 May 2026, seasons 2004-05 to 2025-26", help="period text for meta.what")
    a = ap.parse_args()
    seasons = a.seasons
    z = np.load(a.counts)
    cnt, cut, nt = z["counts"], list(map(float, z["cutoffs_kt"])), z["n_times"]
    ocean = np.load(a.lsm).astype(float) < 0.5                       # 361 x 1440, rows 90N..0N
    lat = 90 - 0.25 * np.arange(361)
    win, summary = {}, {}
    for b, (w0, w1) in WINDOWS.items():
        nlat, nlon = LAT1 - LAT0, w1 - w0
        # row of the southern edge of block i is lat LAT0+i; the block holds the 4 rows going north
        rows = np.array([[int(round((90 - (LAT0 + i)) / 0.25)) - 1 - k for k in range(4)] for i in range(nlat)])
        cols = np.array([[(int(round(((w0 + j) % 360) / 0.25)) + k) % 1440 for k in range(4)] for j in range(nlon)])
        f = {}
        for name, th in CLASSES:
            ki = int(np.argmin(np.abs(np.array(cut) - th)))
            assert abs(cut[ki] - th) < 1e-6
            months = []
            for m in range(12):
                if nt[m] == 0:                                         # month not sampled in this set (1979-2003 set is Oct-Apr)
                    months.append([])
                    continue
                h = cnt[m, ki].astype(float) * 12.0 / seasons         # hours per season in month m
                out = np.full((nlat, nlon), -1, int)
                for i in range(nlat):
                    r = rows[i]
                    for j in range(nlon):
                        c = cols[j]
                        om = ocean[np.ix_(r, c)]
                        if om.any():
                            out[i, j] = int(round(h[np.ix_(r, c)][om].mean() * 10))
                months.append(out.ravel().tolist())
                summary[(b, name, m + 1)] = max(out.max(), 0) / 10
            f[name] = months
        win[b] = {"lat0": LAT0, "lon0": w0, "dlat": 1, "dlon": 1, "nlat": nlat, "nlon": nlon, "f": f}
    meta = {"what": "ERA5 PROXY gust hours per season by calendar month (research/era5/gust_class_maps): instantaneous 10 m gust at "
                    "00 and 12 UTC, " + a.period + ", ocean only, 1 degree blocks. Nominal proxy "
                    "classes on the gust index: gale >= 38 kt, storm >= 54 kt, hurricane force >= 71.7 kt. ERA5 reads low in extreme "
                    "storms and samples are 12-hourly, so every value is a lower bound.",
            "source": a.source, "n_times": int(nt.sum()), "n_times_by_month": nt.tolist(), "seasons": seasons,
            "months": [m + 1 for m in range(12) if nt[m]],
            "cutoffs_kt": {n: t for n, t in CLASSES},
            "fields": ["g38/g54/g717: [12 months] hours per season in that month, 0.1 h units, -1 = no ocean"]}
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("/* generated by tools/build_gust_classes.py - do not edit. ERA5 PROXY, not direct observation. */\n")
        fh.write("window." + a.var + " = " + json.dumps({"meta": meta, "win": win}, separators=(",", ":")) + ";\n")
    print("%.3f MB -> %s" % (os.path.getsize(a.out) / 1e6, a.out))
    for n, _ in CLASSES:
        print(n, "max block h/season by month:", [round(summary.get(("atl", n, m), 0), 1) for m in range(1, 13)], "(atl)")
        print(n, "                              ", [round(summary.get(("pac", n, m), 0), 1) for m in range(1, 13)], "(pac)")


if __name__ == "__main__":
    main()
