"""Candidate predictors for the P(HF) model, from WeatherBench2 1.5 degree ERA5.

Tier 1 of PREREG.md. For every fix in the fixes file with time < 2023-01-10,
compute the columns of groups A, M, D, J, L, P, K (group H comes from the track
history and needs no pull, see features.py). Same grid, discs and cos(lat)
weighting as ../intensity/env.py, whose helpers this file imports.

usage: extract.py FIXES_CSV OUT_DIR [NPROC] [--tier2]
One CSV per WB2 time chunk in OUT_DIR, resumable. `--tier2` adds N, R and W
(stability, trough, omega700); it must not be used before Jason's go-ahead.
"""
import os, sys, time
import numpy as np, pandas as pd, numcodecs
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "intensity"))
import env
from env import LAT2, LON2, R0, R1, WB2, WB2_T0, WB2_END, WB2_LEVELS, KT, cat, wb2meta, wmean, ddx, ddy
import hart

SINGLE = {"integrated_vapor_transport": "ivt", "mean_sea_level_pressure": "msl",
          "2m_temperature": "t2m", "total_precipitation_6hr": "pr6",
          "sea_surface_temperature": "sst"}
LEVELED = {"wind_speed": (250, 850)}                      # Tier 1
LEVELED2 = {"temperature": (925, 700), "geopotential": (500,), "vertical_velocity": (700,)}   # Tier 2
_lsm = None


def lsm():
    """Land-sea mask on the [lat, lon] sub-grid, 1 = land."""
    global _lsm
    if _lsm is None:
        m = wb2meta()["land_sea_mask/.zarray"]
        a = np.frombuffer(numcodecs.get_codec(m["compressor"]).decode(cat(f"{WB2}/land_sea_mask/0.0")),
                          m["dtype"]).reshape(m["chunks"])
        _lsm = a.T[R0:R1].astype(np.float64)
    return _lsm


def read_chunk(chunk, tier2):
    m = wb2meta()
    out = {}
    todo = [(v, None) for v in SINGLE] + list(LEVELED.items())
    if tier2:
        todo += list(LEVELED2.items())
    for var, levs in todo:
        za = m[f"{var}/.zarray"]
        codec = numcodecs.get_codec(za["compressor"])
        key = f"{chunk}.0.0.0" if levs else f"{chunk}.0.0"
        a = np.frombuffer(codec.decode(cat(f"{WB2}/{var}/{key}")), za["dtype"]).reshape(za["chunks"])
        a = a[:min(za["chunks"][0], za["shape"][0] - chunk * za["chunks"][0])]
        if levs:
            for L in levs:
                out[f"{var}{L}"] = a[:, WB2_LEVELS.index(L)].transpose(0, 2, 1)[:, R0:R1]
        else:
            out[SINGLE[var]] = a.transpose(0, 2, 1)[:, R0:R1]
    return out


def displacement(la, lo, la2, lo2):
    """East and north displacement (km) from (la, lo) to (la2, lo2)."""
    n = (la2 - la) * 111.19
    e = ((lo2 - lo + 180) % 360 - 180) * 111.19 * np.cos(np.radians(0.5 * (la + la2)))
    return e, n


def features(F, la, lo, heading, tier2):
    d = hart.great_circle_km(LAT2, LON2, la, lo)
    m3, m5, m8, m10, m15, m20 = (d <= 300), (d <= 500), (d <= 800), (d <= 1000), (d <= 1500), (d <= 2000)
    out = {}
    land = lsm()
    ocean = (land < 0.5) & np.isfinite(F["sst"])
    mo5 = m5 & ocean
    out["airsea"] = wmean(F["sst"] - F["t2m"], mo5)
    out["ivt500"] = wmean(F["ivt"], m5)
    out["ivtmax"] = float(np.nanmax(np.where(m10, F["ivt"], np.nan)))
    out["precip6"] = wmean(F["pr6"], m5)
    ws250 = F["wind_speed250"]
    k = np.nanargmax(np.where(m15, ws250, np.nan))
    i, j = np.unravel_index(k, ws250.shape)
    e, n = displacement(la, lo, LAT2[i, j], LON2[i, j])
    th = np.radians(heading if np.isfinite(heading) else 90.0)
    out["jet_along"] = (e * np.sin(th) + n * np.cos(th)) / 100.0
    out["jet_cross"] = (e * np.cos(th) - n * np.sin(th)) / 100.0
    out["ws850"] = float(np.nanmax(np.where(m5, F["wind_speed850"], np.nan))) * KT
    msl = F["msl"] / 100.0                                   # Pa -> hPa
    if not m3.any():
        m3 = d <= d.min() + 1e-6
    out["dphigh"] = float(np.nanmax(np.where(m20, msl, np.nan)) - np.nanmin(np.where(m3, msl, np.nan)))
    gr = np.hypot(ddx(msl), ddy(msl)) * 1e5                  # hPa per m -> hPa per 100 km
    ring = (d > 300) & (d <= 800)
    out["ringgrad"] = wmean(gr, ring)
    out["landfrac"] = wmean((land > 0.5).astype(float), m5)
    if tier2:
        th925 = F["temperature925"] * (1000 / 925) ** 0.2857
        th700 = F["temperature700"] * (1000 / 700) ** 0.2857
        out["stab"] = wmean((th700 - th925) / 225.0, mo5)
        ann = (d > 1500) & (d <= 2500)
        out["trough"] = wmean(F["geopotential500"] / env.G, ann) - float(np.nanmin(np.where(m10, F["geopotential500"] / env.G, np.nan)))
        out["omega700"] = float(np.nanmin(np.where(m5, F["vertical_velocity700"], np.nan)))
    return out


def run_block(job):
    key, path, fixes, tier2 = job
    if os.path.exists(path):
        return path
    Fall = read_chunk(key, tier2)
    t0 = WB2_T0 + np.timedelta64(int(key) * 8 * 6, "h")
    nt = Fall["msl"].shape[0]
    rows = []
    for i in range(nt):
        t = t0 + np.timedelta64(6 * i, "h")
        stamp = str(t).replace("-", "").replace("T", "")[:10]
        sub = fixes[fixes.time == stamp]
        if not len(sub):
            continue
        F = {k: v[i].astype(np.float64) for k, v in Fall.items()}
        for r in sub.itertuples():
            f = features(F, r.lat, r.lon, r.heading, tier2)
            f.update(track=r.track, time=stamp)
            rows.append(f)
    tmp = path + ".tmp"
    pd.DataFrame(rows).to_csv(tmp, index=False)
    os.replace(tmp, path)
    return path


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    tier2 = "--tier2" in sys.argv
    fixes = pd.read_csv(args[0], dtype={"time": str})
    outdir = args[1]
    nproc = int(args[2]) if len(args) > 2 else 6
    t = pd.to_datetime(fixes.time, format="%Y%m%d%H").values.astype("datetime64[h]")
    fixes = fixes[(t < WB2_END) & (fixes.season >= 2004) & (fixes.season <= 2021)].copy()
    t = pd.to_datetime(fixes.time, format="%Y%m%d%H").values.astype("datetime64[h]")
    fixes["key"] = ((t - WB2_T0) / np.timedelta64(6, "h")).astype(int) // 8
    os.makedirs(outdir, exist_ok=True)
    J = [(int(k), f"{outdir}/wb2_{k:05d}.csv", g.drop(columns="key"), tier2) for k, g in fixes.groupby("key")]
    print(len(fixes), "fixes", len(J), "chunks", flush=True)
    t0 = time.time()
    with Pool(nproc) as p:
        for n, path in enumerate(p.imap_unordered(run_block, J)):
            if n % 100 == 0:
                print(f"{n + 1}/{len(J)} {time.time() - t0:.0f}s", flush=True)
