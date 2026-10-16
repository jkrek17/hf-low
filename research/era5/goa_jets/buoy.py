"""Q3 (descriptive): does the ERA5 PROXY under-read at moored buoys in the Gulf of Alaska at event times?

NDBC standard meteorological history (hourly files, UTC; wind at buoy anemometer height, about 5 m; no height adjustment), nearest
record within 40 min of 00/12 UTC. ERA5 10 m wind from the pulled fields at the nearest cell; ERA5 gust from the shared gust cache.
Station positions are those on the NDBC station pages (check_positions re-reads them); a station whose ERA5 cell is mostly land is reported with its land fraction.
usage: buoy.py check | buoy.py cells  (decode gust at the station cells, slow, once)  |  buoy.py run
"""
import os, sys, glob, gzip, io, urllib.request
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "goa_jets")
CACHE = "/mnt/project-files/era5-cache"
# Positions as published on the NDBC station pages (read 2026-10-16). check_positions() re-reads them and fails on a difference above 0.05 degree.
# The first version of this file had 46077 at 156.94W from memory, 160 km from its true position, which put it in a land cell (failures/2026-10-16-buoy-position-from-memory.md).
ST = {"46001": (56.296, -148.027), "46061": (60.231, -146.840), "46080": (57.910, -150.129), "46076": (59.508, -148.005),
      "46082": (59.670, -143.353), "46083": (58.276, -138.024), "46077": (57.869, -154.211)}
KT = 1.943844   # m/s to kt
T0 = np.datetime64("1900-01-01T00")


def cell(la, lo):
    return int(round((90 - la) / 0.25)), int(round((lo % 360) / 0.25))


def check_positions():
    import re
    for s_, (la, lo) in ST.items():
        h = urllib.request.urlopen(f"https://www.ndbc.noaa.gov/station_page.php?station={s_}", timeout=60).read().decode("latin-1")
        m = re.search(r"([0-9]+\.[0-9]+) ?N ([0-9]+\.[0-9]+) ?W", h)
        pla, plo = float(m.group(1)), -float(m.group(2))
        assert abs(pla - la) < 0.05 and abs(plo - lo) < 0.05, (s_, (la, lo), (pla, plo))
        lsm = np.load(os.path.join(CACHE, "lsm_nh_0p25.npy")); r, c = cell(la, lo)
        print(s_, "NDBC", pla, plo, "ERA5 cell land fraction", round(float(lsm[r, c]), 2))


def download():
    d = os.path.join(WORK, "buoys"); os.makedirs(d, exist_ok=True)
    for s in ST:
        for y in range(2004, 2027):
            f = os.path.join(d, f"{s}_{y}.txt")
            if os.path.exists(f): continue
            try:
                raw = urllib.request.urlopen(f"https://www.ndbc.noaa.gov/data/historical/stdmet/{s}h{y}.txt.gz", timeout=60).read()
                open(f, "wb").write(gzip.decompress(raw))
            except Exception as e:
                open(f, "w").write("")   # missing year


def read_buoy(s):
    parts = []
    for f in sorted(glob.glob(os.path.join(WORK, "buoys", f"{s}_*.txt"))):
        if os.path.getsize(f) == 0: continue
        rows = []; has_mm = False
        for ln in open(f):
            if ln.startswith("#"):
                has_mm = has_mm or (" mm" in ln.lower() and "YY" in ln.upper() or "mm" in ln.split()); continue
            p = ln.split()
            try:
                yr, mo, dy, hh = int(p[0]), int(p[1]), int(p[2]), int(p[3]); mm = int(p[4]) if has_mm else 0
                if yr < 100: yr += 1900 if yr > 50 else 2000
                if has_mm: wd, ws, gs = float(p[5]), float(p[6]), float(p[7])
                else: wd, ws, gs = float(p[4]), float(p[5]), float(p[6])
            except Exception:
                continue
            rows.append((pd.Timestamp(yr, mo, dy, hh, mm), wd, ws, gs))
        parts.append(pd.DataFrame(rows, columns=["t", "wd", "ws", "gs"]))
    if not parts: return pd.DataFrame(columns=["t", "wd", "ws", "gs"])
    d = pd.concat(parts)
    d.loc[d.ws >= 98, "ws"] = np.nan; d.loc[d.gs >= 98, "gs"] = np.nan
    return d.sort_values("t").drop_duplicates("t")


def cells():
    import numcodecs
    from multiprocessing import Pool
    files = sorted(glob.glob(os.path.join(CACHE, "gust_nh_12h", "gust_*.npz")))
    with Pool(4) as p:
        rows = [r for rr in p.imap(_cells_one, files) for r in rr]
    pd.DataFrame(rows, columns=["time", "station", "era5_gust_kt"]).to_csv(os.path.join(RES, "buoy_era5_gust.csv"), index=False)


def _cells_one(f):
    import numcodecs
    codec = numcodecs.Blosc()
    z = np.load(f); t = T0 + z["hours_since_1900"].astype("timedelta64[h]")
    a = np.frombuffer(codec.decode(z["data"].tobytes()), np.uint8).reshape(z["shape"])
    m = t.astype("datetime64[M]").astype(int) % 12
    out = []
    for i in np.where(np.isin(m, [9, 10, 11, 0, 1, 2, 3]))[0]:
        for s, (la, lo) in ST.items():
            r, c = cell(la, lo); out.append((str(t[i]), s, a[i, r, c] / 2.0))
    return out


def run():
    from features import ocean, LA, LO, ROW0, COL0   # same window as the pulled fields
    ff = pd.read_csv(os.path.join(RES, "features.csv")); ff["time"] = pd.to_datetime(ff.time)
    pt = pd.read_csv(os.path.join(RES, "pull_times.csv")); pt["time"] = pd.to_datetime(pt.time)
    g = pd.read_csv(os.path.join(RES, "buoy_era5_gust.csv")); g["time"] = pd.to_datetime(g.time); g["station"] = g.station.astype(str)
    lsm = np.load(os.path.join(CACHE, "lsm_nh_0p25.npy"))
    rows = []
    for s, (la, lo) in ST.items():
        b = read_buoy(s)
        r, c = cell(la, lo)
        land_cell = bool(lsm[r, c] >= 0.5); land_frac = float(lsm[r, c])
        for t in pt.time:
            st = t.strftime("%Y%m%d%H"); z = np.load(os.path.join(WORK, "fields", st + ".npz"))
            rr, cc = r - ROW0, c - COL0
            if not (0 <= rr < z["u10"].shape[0] and 0 <= cc < z["u10"].shape[1]): continue
            e10 = float(np.hypot(z["u10"][rr, cc], z["v10"][rr, cc]) * KT)
            w = b[(b.t - t).abs() <= pd.Timedelta(minutes=40)]
            if w.empty: continue
            w = w.iloc[(w.t - t).abs().argmin()]
            eg = g[(g.time == t) & (g.station == s)].era5_gust_kt
            rows.append(dict(time=t, station=s, land_cell=land_cell, land_frac=land_frac, buoy_ws_kt=w.ws * KT, buoy_gust_kt=w.gs * KT, era5_wind_kt=e10,
                             era5_gust_kt=float(eg.iloc[0]) if len(eg) else np.nan))
    B = pd.DataFrame(rows)
    z = ff[ff.zone.isin(["Z1", "Z2", "Z3"])][["time", "zone", "regime", "event"]] if "regime" in ff else None
    B.to_csv(os.path.join(RES, "buoy_pairs.csv"), index=False)
    out = []
    for s, d in B.groupby("station"):
        for lab, m in (("all pulled times", d.buoy_ws_kt.notna()), ("buoy wind >= 20 kt", d.buoy_ws_kt >= 20), ("buoy wind >= 30 kt", d.buoy_ws_kt >= 30), ("buoy gust >= 45 kt", d.buoy_gust_kt >= 45)):
            x = d[m]
            if len(x) < 5: continue
            out.append(dict(station=s, subset=lab, n=len(x), land_cell=bool(d.land_cell.iloc[0]),
                            median_ratio_wind=float((x.era5_wind_kt / x.buoy_ws_kt).median()),
                            median_ratio_gust=float((x.era5_gust_kt / x.buoy_gust_kt).median()),
                            era5_gust_ge54=int((x.era5_gust_kt >= 54).sum()), buoy_gust_ge54=int((x.buoy_gust_kt >= 54).sum())))
    O = pd.DataFrame(out); O.to_csv(os.path.join(RES, "buoy_summary.csv"), index=False)
    print(O.round(2).to_string(index=False))
    # when ERA5's own cell shows a storm-force gust, what did the buoy report?
    rr = []
    for s_, d in B.groupby("station"):
        x = d[d.era5_gust_kt >= 54]
        if len(x) >= 5 and not bool(d.land_cell.iloc[0]):
            rr.append(dict(station=s_, n=len(x), median_era5_gust=float(x.era5_gust_kt.median()), median_buoy_gust=float(x.buoy_gust_kt.median()),
                           share_buoy_gust_ge48=float((x.buoy_gust_kt >= 48).mean()), share_buoy_gust_ge54=float((x.buoy_gust_kt >= 54).mean())))
    R2 = pd.DataFrame(rr); R2.to_csv(os.path.join(RES, "buoy_when_era5_54.csv"), index=False)
    print(R2.round(2).to_string(index=False))


if __name__ == "__main__":
    sys.path.insert(0, HERE)
    if sys.argv[1] == "download": download()
    elif sys.argv[1] == "check": check_positions()
    elif sys.argv[1] == "cells": cells()
    else: run()
