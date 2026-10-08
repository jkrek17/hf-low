"""Test 1 stage B extraction: 0.25 degree storm-relative environment at t0, t0-12, t0-24 h of converters and non-converters. ERA5 proxy, pipeline A low detection (helpers from hf_vs_storm).
Fields per time: MSLP, 2 m dewpoint, 2 m temperature, SST (surface); 250 hPa u,v; 850 hPa T, q (single levels by blosc range requests). Gust is not pulled (the 800 km owned gust is not recomputed; the position check stands in).
Position = the track's 00/12 UTC fix snapped to the nearest re-detected low within 150 km, else missing (lag 0 additionally needs <= 25 km at analysis). usage: HFVS_SUB=hf_conversion_pull ERA5_WORK=... extract_conv.py NPROC [--dry|N]"""
import os, sys, time, datetime as dt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "hf_vs_storm"))
os.environ.setdefault("HFVS_SUB", "hf_conversion_pull")
import numpy as np, pandas as pd
from multiprocessing import Pool
import extract_storm as E
import fields_large as FL
WORK = E.WORK; SNAP = 150.0
DY = np.radians(0.25) * E.RE                                  # km per 0.25 degree of latitude


def thetae(T, q, p=850.0):
    """Bolton (1980) equivalent potential temperature (K) from T (K), specific humidity (kg/kg) at p (hPa)."""
    r = q / (1 - q) * 1000.0                                   # g/kg
    e = np.maximum(q * p / (0.622 + 0.378 * q), 1e-3)
    ln = np.log(e / 6.112); Td = 243.5 * ln / (17.67 - ln) + 273.15
    TL = 1 / (1 / (Td - 56) + np.log(T / Td) / 800) + 56
    return T * (1000 / p) ** (0.2854 * (1 - 0.28e-3 * r)) * np.exp((3.376 / TL - 0.00254) * r * (1 + 0.81e-3 * r))


def one_time(job):
    stamp, fixes = job
    spath = os.path.join(WORK, "stats", f"{stamp}.csv")
    if os.path.exists(spath):
        return stamp, 0
    s = stamp; ti = int((dt.datetime(int(s[:4]), int(s[4:6]), int(s[6:8]), int(s[8:10])) - dt.datetime(1900, 1, 1)).total_seconds() // 3600); nbytes = 0; F = {}
    for name, var in (("msl", "mean_sea_level_pressure"), ("d2m", "2m_dewpoint_temperature"), ("t2m", "2m_temperature"), ("sst", "sea_surface_temperature")):
        F[name], nb = E.fetch(f"{E.ARCO}/{var}/{ti}.0.0"); nbytes += nb
    for name, var, lev in (("u250", "u_component_of_wind", 250), ("v250", "v_component_of_wind", 250), ("t850", "temperature", 850), ("q850", "specific_humidity", 850)):
        F[name], nb = FL.arco_level(var, ti, lev); nbytes += nb
    msl = F["msl"] / 100.0; d2m = F["d2m"] - 273.15; t2m = F["t2m"] - 273.15; sst = F["sst"] - 273.15
    ws250 = np.hypot(F["u250"], F["v250"]); te = thetae(F["t850"], F["q850"])
    stab = sst - t2m                                           # NaN over land (SST undefined)
    gy, gx = np.gradient(F["t850"], DY, axis=0), np.gradient(F["t850"], axis=1) / (DY * np.cos(np.radians(E.LAT))[:, None])
    gt = np.hypot(gy, gx) * 100.0                              # K per 100 km
    lows = E.detect(msl); cen = np.array([E.XYZ[i, j] for _, i, j in lows])
    rows, boxes = [], {}
    for fx in fixes:
        v = E.XYZ[int(round((90 - fx["lat"]) / 0.25)), int(round(fx["lon"] / 0.25)) % 1440]
        dk = np.arccos(np.clip(cen @ v, -1, 1)) * E.RE; k = int(np.argmin(dk))
        base = dict(time=stamp, track=fx["track"], anchor=fx["anchor"], basin=fx["basin"], lat=fx["lat"], lon=fx["lon"], heading=fx["heading"])
        if dk[k] > SNAP:
            rows.append(dict(base, snapped=False)); continue
        _, ci, cj = lows[k]; clat, clon = E.LAT[ci], E.LON[cj]; hd = fx["heading"]
        bn = {n: E.box(f, clat, clon, 90.0) for n, f in (("msl", msl), ("ws250", ws250), ("te", te), ("stab", stab), ("d2m", d2m), ("gt", gt))}
        ok = np.isfinite(bn["ws250"]) & (E.BR <= 1500)
        j = np.argmax(np.where(ok, bn["ws250"], -1)); jy, jx = np.unravel_index(j, ok.shape)
        jet_d = float(E.BR[jy, jx]); jet_v = float(bn["ws250"][jy, jx])
        if np.isfinite(hd) and jet_d > 0:
            brg = (np.degrees(np.arctan2(E.BX[jy, jx], E.BY[jy, jx]))) % 360                 # north-up box: +x east, +y north
            jet_cos = float(np.cos(np.radians(brg - hd)))
        else:
            jet_cos = np.nan
        m5 = E.BR <= 500; oc = m5 & np.isfinite(bn["stab"])
        rows.append(dict(base, snapped=True, match_km=round(float(dk[k]), 1), jet_dist=round(jet_d, 0), jet_speed=round(jet_v, 1), jet_cos=round(jet_cos, 3) if np.isfinite(jet_cos) else np.nan,
                         thetae850=round(float(np.nanmean(bn["te"][m5])), 2), stab=round(float(np.nanmean(bn["stab"][oc])), 3) if oc.sum() >= 50 else np.nan, n_ocean=int(oc.sum()),
                         baroc=round(float(np.nanmean(bn["gt"][m5])), 3), pc=round(float(np.nanmin(np.where(E.BR <= 100, bn["msl"], np.nan))), 1), d2m_500=round(float(np.nanmean(bn["d2m"][m5])), 2)))
        if np.isfinite(hd):
            for n in ("msl", "ws250", "te", "stab", "d2m", "gt"):
                boxes[f"{fx['track']}_{fx['anchor']}_{n}"] = E.box({"msl": msl, "ws250": ws250, "te": te, "stab": stab, "d2m": d2m, "gt": gt}[n], clat, clon, hd).astype(np.float16)
    np.savez_compressed(os.path.join(WORK, "boxes", f"{stamp}.npz"), **boxes)
    pd.DataFrame(rows).to_csv(spath + ".tmp", index=False); os.replace(spath + ".tmp", spath)
    with open(os.path.join(WORK, "bytes.log"), "a") as f:
        f.write(f"{stamp},{nbytes}\n")
    return stamp, nbytes


def jobs():
    T = pd.read_csv(os.environ.get("HFVS_TIMES", os.path.join(HERE, "results", "times_conv.csv")))
    return [(str(t), g[["track", "anchor", "basin", "lat", "lon", "heading"]].to_dict("records")) for t, g in T.groupby("time")]


if __name__ == "__main__":
    J = jobs()
    if "--dry" in sys.argv:
        print(len(J), "times", sum(len(j[1]) for j in J), "anchors"); sys.exit()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
    if lim: J = J[:lim]
    t0, tot = time.time(), 0
    with Pool(n) as p:
        for i, (s, nb) in enumerate(p.imap_unordered(one_time, J, chunksize=2)):
            tot += nb
            if i % 100 == 0: print(f"{i + 1}/{len(J)} {tot / 1e9:.2f} GB {time.time() - t0:.0f}s", flush=True)
    print(f"done {len(J)} times, {tot / 1e9:.2f} GB this run")
