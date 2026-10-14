#!/usr/bin/env python3
"""Storm-scale environment diagnostics along one storm's track, for the storm-review pages.

For each pipeline A fix of a storm, compute from WeatherBench2 ERA5 at 1.5 degrees (6-hourly, 13 pressure levels):
  B        Hart thermal asymmetry, 900-600 hPa thickness right minus left of motion, R 500 km (m)
  VTL,VTU  Hart -V_T, lower (900-600 hPa) and upper (600-300 hPa) troposphere, R 500 km (m); positive = warm core
  jet250   highest 250 hPa wind within 1000 km (kt)
  div300   mean 300 hPa divergence within 500 km (1e-5 /s)
  vadv500  mean 500 hPa absolute-vorticity advection within 500 km (1e-9 /s2); positive = cyclonic vorticity advected in
  eady     mean 850-500 hPa Eady growth rate within 500 km (per day)
  sst, sstgrad, sst_t500, flux, tcwv   sea-surface temperature (C), its largest gradient (K per 100 km), SST minus 500 hPa
           temperature (K), upward surface heat flux (W m-2), total column water vapour (kg m-2), all within 500 km
Same definitions as research/era5/intensity/env.py on claude/exciting-fermat-8vcvgq (8fa87cc), reduced to the WeatherBench2
branch (storms up to 2023-01-09) and made importable. ERA5 proxy; pipeline A tracks. Not observation.

usage: storm_review_env.py FIXES_CSV OUT_CSV      FIXES_CSV needs track,time(YYYYMMDDHH),lat,lon[,heading]
"""
import os, sys, time, json, urllib.request
import numpy as np, pandas as pd, numcodecs

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hart_cps as hart

WB2 = "weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
WB2_T0 = np.datetime64("1959-01-01T00")
G = 9.80665
RE = 6.371e6
OMEGA = 7.292e-5
KT = 1 / 0.514444
LAT = -90 + 1.5 * np.arange(121)
LON = 1.5 * np.arange(240)
R0, R1 = 70, 116                              # 15N .. 82.5N
LATS = LAT[R0:R1]
LAT2, LON2 = np.meshgrid(LATS, LON, indexing="ij")
PL = {"geopotential": (925, 850, 700, 600, 500, 400, 300),
      "temperature": (850, 500),
      "u_component_of_wind": (850, 500, 300, 250),
      "v_component_of_wind": (850, 500, 300, 250)}
SL = ("sea_surface_temperature", "mean_surface_latent_heat_flux", "mean_surface_sensible_heat_flux", "total_column_water_vapour")
WB2_LEVELS = (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000)
SHORT = {"geopotential": "z", "temperature": "t", "u_component_of_wind": "u", "v_component_of_wind": "v",
         "sea_surface_temperature": "sst", "mean_surface_latent_heat_flux": "lhf",
         "mean_surface_sensible_heat_flux": "shf", "total_column_water_vapour": "tcwv"}
_meta = None


def cat(key):
    for k in range(6):
        try:
            with urllib.request.urlopen("https://storage.googleapis.com/" + key, timeout=180) as r:
                return r.read()
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {key}")


def wb2meta():
    global _meta
    if _meta is None:
        _meta = json.loads(cat(f"{WB2}/.zmetadata"))["metadata"]
    return _meta


def read_wb2(chunk):
    """All 8 times of one WB2 time chunk -> {name: array[ntime, len(LATS), 240]}; SI units except z in m."""
    m = wb2meta()
    out = {}
    for var, levs in list(PL.items()) + [(v, None) for v in SL]:
        za = m[f"{var}/.zarray"]
        codec = numcodecs.get_codec(za["compressor"])
        key = f"{chunk}.0.0.0" if levs else f"{chunk}.0.0"
        shape = list(za["chunks"])
        a = np.frombuffer(codec.decode(cat(f"{WB2}/{var}/{key}")), za["dtype"]).reshape(shape)
        a = a[:min(shape[0], za["shape"][0] - chunk * shape[0])]
        if levs:
            for L in levs:
                f = a[:, WB2_LEVELS.index(L)].transpose(0, 2, 1)[:, R0:R1]
                out[f"{SHORT[var]}{L}"] = f / G if var == "geopotential" else f
        else:
            out[SHORT[var]] = a.transpose(0, 2, 1)[:, R0:R1]
    return out


def ddx(f):
    dl = np.radians(1.5)
    return (np.roll(f, -1, -1) - np.roll(f, 1, -1)) / (2 * dl * RE * np.cos(np.radians(LAT2)))


def ddy(f):
    d = np.full_like(f, np.nan)
    d[1:-1] = (f[2:] - f[:-2]) / (2 * np.radians(1.5) * RE)
    return d


def derived(F):
    cosl = np.cos(np.radians(LAT2))
    fcor = 2 * OMEGA * np.sin(np.radians(LAT2))
    D = {}
    D["z900"] = F["z925"] + (F["z850"] - F["z925"]) * np.log(925 / 900) / np.log(925 / 850)
    D["ws250"] = np.hypot(F["u250"], F["v250"]) * KT
    D["div300"] = (ddx(F["u300"]) + ddy(F["v300"] * cosl) / cosl) * 1e5
    zeta = ddx(F["v500"]) - ddy(F["u500"] * cosl) / cosl + fcor
    D["vadv500"] = -(F["u500"] * ddx(zeta) + F["v500"] * ddy(zeta)) * 1e9
    th850 = F["t850"] * (1000 / 850) ** 0.2857
    th500 = F["t500"] * (1000 / 500) ** 0.2857
    dz = F["z500"] - F["z850"]
    n2 = G / (0.5 * (th850 + th500)) * (th500 - th850) / dz
    shear = np.hypot(F["u500"] - F["u850"], F["v500"] - F["v850"]) / dz
    with np.errstate(invalid="ignore"):
        D["eady"] = 0.3098 * np.abs(fcor) * shear / np.sqrt(np.where(n2 > 1e-6, n2, 1e-6)) * 86400
    D["sst"] = F["sst"] - 273.15
    D["sstgrad"] = np.hypot(ddx(F["sst"]), ddy(F["sst"])) * 1e5
    D["t500"] = F["t500"] - 273.15
    D["flux"] = -(F["lhf"] + F["shf"])          # ECMWF fluxes are positive downward
    return D


def wmean(x, m):
    m = m & np.isfinite(x)
    if not m.any():
        return np.nan
    w = np.cos(np.radians(LAT2[m]))
    return float(np.sum(x[m] * w) / np.sum(w))


def features(F, D, la, lo, heading):
    d = hart.great_circle_km(LAT2, LON2, la, lo)
    m5, m10 = d <= 500, d <= 1000
    zs = np.stack([D["z900"], F["z850"], F["z700"], F["z600"], F["z500"], F["z400"], F["z300"]])
    tw = hart.thermal_wind((900, 850, 700, 600, 500, 400, 300), zs, LAT2, LON2, la, lo)
    hd = heading if np.isfinite(heading) else 90.0
    sg = np.where(m5, D["sstgrad"], np.nan)
    out = dict(B=hart.parameter_b(D["z900"], F["z600"], LAT2, LON2, la, lo, hd), VTL=tw["VTL"], VTU=tw["VTU"],
               jet250=float(np.nanmax(np.where(m10, D["ws250"], np.nan))),
               div300=wmean(D["div300"], m5), vadv500=wmean(D["vadv500"], m5), eady=wmean(D["eady"], m5), sst=wmean(D["sst"], m5),
               sstgrad=float(np.nanmax(sg)) if np.isfinite(sg).any() else np.nan, flux=wmean(D["flux"], m5), tcwv=wmean(F["tcwv"], m5))
    out["sst_t500"] = out["sst"] - wmean(D["t500"], m5)
    return out


def env_for_fixes(fixes):
    """fixes: DataFrame(track, time 'YYYYMMDDHH', lat, lon, heading) at 6-hourly times before 2023-01-10 -> DataFrame of features."""
    t = pd.to_datetime(fixes.time, format="%Y%m%d%H").values.astype("datetime64[h]")
    chunk = ((t - WB2_T0) / np.timedelta64(6, "h")).astype(int) // 8
    rows = []
    for ck, g in fixes.groupby(chunk):
        Fall = read_wb2(int(ck))
        t0 = WB2_T0 + np.timedelta64(int(ck) * 8 * 6, "h")
        for i in range(Fall["z500"].shape[0]):
            stamp = str(t0 + np.timedelta64(6 * i, "h")).replace("-", "").replace("T", "")[:10]
            sub = g[g.time == stamp]
            if not len(sub):
                continue
            F = {k: v[i].astype(np.float64) for k, v in Fall.items()}
            D = derived(F)
            for r in sub.itertuples():
                f = features(F, D, r.lat, r.lon, r.heading)
                f.update(track=r.track, time=stamp)
                rows.append(f)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    fx = pd.read_csv(sys.argv[1], dtype={"time": str})
    env_for_fixes(fx).to_csv(sys.argv[2], index=False)
