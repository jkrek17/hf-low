"""Near-storm environment along ERA5 cyclone tracks, on a 1.5 degree grid.

For every fix in a fixes CSV (time, lat, lon, heading), compute Hart's cyclone
phase space terms and a set of environmental predictors from ERA5 at 1.5
degrees, the resolution class Hart (2003) used.

Sources
- 1979 .. 2023-01-09: WeatherBench2's conservative 1.5 degree regrid of ERA5
  (13 pressure levels, 6-hourly, 8 times per chunk).
- 2023-01-10 on: the hourly 0.25 degree ARCO-ERA5 store, coarsened here to the
  same 1.5 degree grid with the same cell edges (see coarsen()).
Both are the same reanalysis; check_overlap() compares them where they meet.

Predictors (all centred on the fix; R = radius)
  B        Hart thermal asymmetry, 900-600 hPa thickness right minus left of
           motion, R 500 km (m). 900 hPa is interpolated in ln p from 925/850.
  VTL,VTU  Hart -V_T, slope of max-min height against ln p over 900-600 and
           600-300 hPa, R 500 km (m). Levels 900,850,700,600 / 600,500,400,300,
           not Hart's 50 hPa spacing (hart.py is a copy of awips-tools' module).
  jet250   max 250 hPa wind speed, R 1000 km (kt)
  div300   mean 300 hPa divergence, R 500 km (1e-5 s-1)
  vadv500  mean 500 hPa absolute vorticity advection, R 500 km (1e-9 s-2)
  eady     mean 850-500 hPa Eady growth rate, R 500 km (day-1)
  sst      mean SST over ocean, R 500 km (C)
  sstgrad  max |grad SST|, R 500 km (K per 100 km)
  sst_t500 sst minus mean 500 hPa temperature, R 500 km (K)
  flux     mean upward surface sensible + latent heat flux, R 500 km (W m-2)
  tcwv     mean total column water vapour, R 500 km (kg m-2)

usage: env.py FIXES_CSV OUT_DIR [NPROC]
Writes one CSV per source chunk into OUT_DIR, resumable; cat them afterwards.
"""
import os, sys, time
import numpy as np, pandas as pd, numcodecs, gcsfs, json
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hart

WB2 = "weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
ARCO = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
WB2_T0 = np.datetime64("1959-01-01T00")
ARCO_T0 = np.datetime64("1900-01-01T00")
WB2_END = np.datetime64("2023-01-10T00")      # first time not in WB2 chunks we trust
G = 9.80665
RE = 6.371e6
OMEGA = 7.292e-5
KT = 1 / 0.514444

LAT = -90 + 1.5 * np.arange(121)
LON = 1.5 * np.arange(240)
R0, R1 = 70, 116                              # 15N .. 82.5N, enough for 20-75N fixes + 1000 km
LATS = LAT[R0:R1]
LAT2, LON2 = np.meshgrid(LATS, LON, indexing="ij")

PL = {"geopotential": (925, 850, 700, 600, 500, 400, 300),
      "temperature": (850, 500),
      "u_component_of_wind": (850, 500, 300, 250),
      "v_component_of_wind": (850, 500, 300, 250)}
SL = ("sea_surface_temperature", "mean_surface_latent_heat_flux",
      "mean_surface_sensible_heat_flux", "total_column_water_vapour")
WB2_LEVELS = (50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000)
ARCO_LEVELS = (1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400,
               450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000)

_fs = None
_wb2meta = None


def fs():
    global _fs
    if _fs is None:
        _fs = gcsfs.GCSFileSystem(token="anon")
    return _fs


def cat(key):
    for k in range(6):
        try:
            return fs().cat(key)
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError(f"fetch failed {key}")


def wb2meta():
    global _wb2meta
    if _wb2meta is None:
        _wb2meta = json.loads(cat(f"{WB2}/.zmetadata"))["metadata"]
    return _wb2meta


# ---------------------------------------------------------------- field readers
# Every reader returns {name: array[ntime, len(LATS), 240]} with names like
# "z500", "t850", "u250", "sst", "lhf", "shf", "tcwv"; SI units except z in m.

SHORT = {"geopotential": "z", "temperature": "t", "u_component_of_wind": "u", "v_component_of_wind": "v",
         "sea_surface_temperature": "sst", "mean_surface_latent_heat_flux": "lhf",
         "mean_surface_sensible_heat_flux": "shf", "total_column_water_vapour": "tcwv"}


def read_wb2(chunk):
    """All 8 times of one WB2 time chunk. WB2 arrays are (time, [level,] lon, lat)."""
    m = wb2meta()
    out = {}
    for var, levs in list(PL.items()) + [(v, None) for v in SL]:
        za = m[f"{var}/.zarray"]
        codec = numcodecs.get_codec(za["compressor"])
        key = f"{chunk}.0.0.0" if levs else f"{chunk}.0.0"
        shape = list(za["chunks"])
        a = np.frombuffer(codec.decode(cat(f"{WB2}/{var}/{key}")), za["dtype"]).reshape(shape)
        nt = min(shape[0], za["shape"][0] - chunk * shape[0])
        a = a[:nt]
        if levs:
            for L in levs:
                f = a[:, WB2_LEVELS.index(L)].transpose(0, 2, 1)[:, R0:R1]
                out[f"{SHORT[var]}{L}"] = f / G if var == "geopotential" else f
        else:
            out[SHORT[var]] = a.transpose(0, 2, 1)[:, R0:R1]
    return out


# Conservative 0.25 -> 1.5 degree: a 1.5 degree cell centred on a 0.25 degree
# point spans 7 points, the two end ones shared half with the neighbour.
_W7 = np.array([0.5, 1, 1, 1, 1, 1, 0.5]) / 6.0


def coarsen(a):
    """0.25 degree [721, 1440] (lat 90..-90) -> 1.5 degree [121, 240] (lat -90..90)."""
    a = a.astype(np.float64)
    nan = np.isnan(a)
    w = (~nan).astype(np.float64)
    a = np.where(nan, 0.0, a)
    def smooth(x):
        # longitude, periodic: cells centred on 0, 1.5, ... -> 0.25 index 0, 6, ...
        xl = np.concatenate([x[:, -3:], x, x[:, :3]], axis=1)
        y = sum(_W7[k] * xl[:, k:k + 1440] for k in range(7))[:, ::6]
        # latitude: centres at 90, 88.5, ... ; pad by edge reflection at the poles
        yl = np.concatenate([y[3:0:-1], y, y[-2:-5:-1]], axis=0)
        z = sum(_W7[k] * yl[k:k + 721] for k in range(7))[::6]
        return z
    num, den = smooth(a), smooth(w)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0.5 * smooth(np.ones_like(w)), num / den, np.nan)
    return out[::-1]


def read_arco(ti):
    """One hourly ARCO time index, coarsened to 1.5 degrees."""
    codec = numcodecs.Blosc()
    out = {}
    for var, levs in PL.items():
        a = np.frombuffer(codec.decode(cat(f"{ARCO}/{var}/{ti}.0.0.0")), "<f4").reshape(37, 721, 1440)
        for L in levs:
            f = coarsen(a[ARCO_LEVELS.index(L)])[R0:R1][None]
            out[f"{SHORT[var]}{L}"] = f / G if var == "geopotential" else f
    for var in SL:
        a = np.frombuffer(codec.decode(cat(f"{ARCO}/{var}/{ti}.0.0")), "<f4").reshape(721, 1440)
        out[SHORT[var]] = coarsen(a)[R0:R1][None]
    return out


# ---------------------------------------------------------------- diagnostics

def ddx(f):
    """d/dx on the sphere, periodic in longitude (per m)."""
    dl = np.radians(1.5)
    return (np.roll(f, -1, -1) - np.roll(f, 1, -1)) / (2 * dl * RE * np.cos(np.radians(LAT2)))


def ddy(f):
    d = np.full_like(f, np.nan)
    d[1:-1] = (f[2:] - f[:-2]) / (2 * np.radians(1.5) * RE)
    return d


def derived(F):
    """2-D derived fields for one time from the {name: [lat, lon]} dict."""
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
    sst = F["sst"] - 273.15
    D["sst"] = sst
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
    out = dict(B=hart.parameter_b(D["z900"], F["z600"], LAT2, LON2, la, lo, hd),
               VTL=tw["VTL"], VTU=tw["VTU"],
               jet250=float(np.nanmax(np.where(m10, D["ws250"], np.nan))),
               div300=wmean(D["div300"], m5), vadv500=wmean(D["vadv500"], m5),
               eady=wmean(D["eady"], m5), sst=wmean(D["sst"], m5),
               sstgrad=float(np.nanmax(np.where(m5, D["sstgrad"], np.nan))) if np.isfinite(np.where(m5, D["sstgrad"], np.nan)).any() else np.nan,
               flux=wmean(D["flux"], m5), tcwv=wmean(F["tcwv"], m5))
    out["sst_t500"] = out["sst"] - wmean(D["t500"], m5)
    return out


def run_block(job):
    src, key, path, fixes = job
    if os.path.exists(path):
        return path
    if src == "wb2":
        Fall = read_wb2(key)
        t0 = WB2_T0 + np.timedelta64(int(key) * 8 * 6, "h")
        times = [t0 + np.timedelta64(6 * i, "h") for i in range(Fall["z500"].shape[0])]
    else:
        Fall = read_arco(key)
        times = [ARCO_T0 + np.timedelta64(int(key), "h")]
    rows = []
    for i, t in enumerate(times):
        stamp = str(t).replace("-", "").replace("T", "")[:10]
        sub = fixes[fixes.time == stamp]
        if not len(sub):
            continue
        F = {k: v[i].astype(np.float64) for k, v in Fall.items()}
        D = derived(F)
        for r in sub.itertuples():
            f = features(F, D, r.lat, r.lon, r.heading)
            f.update(track=r.track, time=stamp)
            rows.append(f)
    tmp = path + ".tmp"
    pd.DataFrame(rows).to_csv(tmp, index=False)
    os.replace(tmp, path)
    return path


def jobs(fixes, outdir):
    t = pd.to_datetime(fixes.time, format="%Y%m%d%H").values.astype("datetime64[h]")
    wb = t < WB2_END
    key = np.where(wb, ((t - WB2_T0) / np.timedelta64(6, "h")).astype(int) // 8,
                   ((t - ARCO_T0) / np.timedelta64(1, "h")).astype(int))
    out = []
    for (w, k), g in fixes.groupby([wb, key]):
        src = "wb2" if w else "arco"
        name = f"wb2_{k:05d}" if w else f"arco_{k}"
        out.append((src, int(k), f"{outdir}/{name}.csv", g))
    out.sort(key=lambda j: j[0] != "arco")      # ARCO (slowest) first
    return out


if __name__ == "__main__":
    fixes = pd.read_csv(sys.argv[1], dtype={"time": str})
    outdir = sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    os.makedirs(outdir, exist_ok=True)
    J = jobs(fixes, outdir)
    t0 = time.time()
    with Pool(nproc) as p:
        for i, path in enumerate(p.imap_unordered(run_block, J)):
            if i % 50 == 0:
                print(f"{i + 1}/{len(J)} {time.time() - t0:.0f}s {os.path.basename(path)}", flush=True)
