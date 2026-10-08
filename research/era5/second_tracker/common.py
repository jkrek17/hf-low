"""Shared pieces for the second-tracker study (RA-15): WeatherBench2 1.5 degree ERA5 chunk access,
2D extremum detection with biquadratic refinement, and track linking for trackers M and V.

All of it is a proxy built on ERA5 at 1.5 degrees; nothing here uses pipeline A.
Grid: arrays are [time, lon(240), lat(121)], lat ascending from -90, lon 0..358.5, 6-hourly from 1959-01-01 00 UTC,
chunks of 8 steps (2 days).
"""
import time, urllib.request, datetime as dt
import numpy as np, numcodecs
from scipy import ndimage
from scipy.optimize import linear_sum_assignment

WB2 = "https://storage.googleapis.com/weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
T0 = dt.datetime(1959, 1, 1)
CODEC = numcodecs.Blosc()
BYTES = [0]
H = 1.5
LAT = -90 + H * np.arange(121)
LON = H * np.arange(240)
R0, R1 = 73, 111              # rows 19.5N..75N (20-75N domain plus one row of margin)
LATD = LAT[R0:R1]
RE = 6371.0
LEV850 = 10                   # level index of 850 hPa in the 13-level axis


def chunk_of(t):
    """chunk index and in-chunk step of datetime t (6-hourly grid)."""
    h = int((t - T0).total_seconds() // 3600)
    assert h % 6 == 0
    return h // 48, (h % 48) // 6


def stamp(c, k):
    t = T0 + dt.timedelta(hours=48 * c + 6 * k)
    return int(t.strftime("%Y%m%d%H"))


def fetch(var, c):
    """one 2-day chunk; MSLP -> (8, lat, lon) in hPa; vorticity 850 -> (8, lat, lon) in 1/s. Retries."""
    key = f"{c}.0.0" if var == "mean_sea_level_pressure" else f"{c}.0.0.0"
    for k in range(6):
        try:
            with urllib.request.urlopen(f"{WB2}/{var}/{key}", timeout=120) as r:
                raw = r.read()
            BYTES[0] += len(raw)
            a = np.frombuffer(CODEC.decode(raw), "<f4")
            if var == "mean_sea_level_pressure":
                a = a.reshape(8, 240, 121) / 100.0
            else:
                a = a.reshape(8, 13, 240, 121)[:, LEV850]
            return np.ascontiguousarray(a.transpose(0, 2, 1))[:, R0:R1, :]
        except Exception:
            if k == 5:
                raise
            time.sleep(2 ** k)


def refine(f, i, j, sign):
    """biquadratic fit through the 3x3 around (i, j) (periodic in lon); returns (lat, lon, value)
    of the extremum; sign=+1 for a minimum, -1 for a maximum (f is negated first)."""
    g = sign * f
    ii = np.clip([i - 1, i, i + 1], 0, g.shape[0] - 1)
    jj = np.array([j - 1, j, j + 1]) % g.shape[1]
    v = g[np.ix_(ii, jj)]
    b = (v[1, 2] - v[1, 0]) / 2
    c = (v[2, 1] - v[0, 1]) / 2
    d = (v[1, 2] - 2 * v[1, 1] + v[1, 0]) / 2
    ff = (v[2, 1] - 2 * v[1, 1] + v[0, 1]) / 2
    e = (v[2, 2] - v[2, 0] - v[0, 2] + v[0, 0]) / 4
    M = np.array([[2 * d, e], [e, 2 * ff]])
    try:
        x, y = np.clip(np.linalg.solve(M, [-b, -c]), -1, 1) if abs(np.linalg.det(M)) > 1e-12 else (0.0, 0.0)
    except np.linalg.LinAlgError:
        x, y = 0.0, 0.0
    val = v[1, 1] + b * x + c * y + d * x * x + e * x * y + ff * y * y
    return LATD[i] + y * H, (LON[j] + x * H) % 360, sign * val


def detect_m(msl):
    """tracker M detector on one [lat, lon] MSLP field in hPa: lowest within 5 grid points (11x11),
    Laplacian 4*(mean of 8 neighbours - centre)/H^2 >= 0.2 hPa/deg^2, pressure < 1020 hPa."""
    mn = ndimage.minimum_filter(msl, size=11, mode=["nearest", "wrap"])
    nb = ndimage.uniform_filter(msl, size=3, mode=["nearest", "wrap"]) * 9
    lap = 4 * ((nb - msl) / 8 - msl) / H ** 2
    r, c = np.where((msl == mn) & (lap >= 0.2) & (msl < 1020.0))
    keep = (r > 0) & (r < msl.shape[0] - 1)
    out = []
    for i, j in zip(r[keep], c[keep]):
        la, lo, p = refine(msl, i, j, +1)
        if 20 <= la <= 75:
            out.append((round(la, 3), round(lo, 3), round(float(p), 2)))
    return out


def detect_v(zeta):
    """tracker V detector on one [lat, lon] 850 hPa relative vorticity field (1/s): Gaussian sigma 1 grid point,
    maxima >= 1e-5 /s that are the largest within 5 grid points."""
    sm = ndimage.gaussian_filter(zeta, 1, mode=["nearest", "wrap"])
    mx = ndimage.maximum_filter(sm, size=11, mode=["nearest", "wrap"])
    r, c = np.where((sm == mx) & (sm >= 1e-5))
    keep = (r > 0) & (r < sm.shape[0] - 1)
    out = []
    for i, j in zip(r[keep], c[keep]):
        la, lo, z = refine(sm, i, j, -1)
        if 20 <= la <= 75:
            out.append((round(la, 3), round(lo, 3), float(z)))
    return out


def dist(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def _xy(la0, lo0, la, lo):
    """east, north displacement in km of (la, lo) from (la0, lo0) (tangent plane)."""
    dlo = ((np.asarray(lo) - lo0 + 180) % 360) - 180
    return RE * np.radians(dlo) * np.cos(np.radians(la0)), RE * np.radians(np.asarray(la) - la0)


def link(det, mode, min_len=4):
    """link detections [(time_int, lat, lon, val)] sorted by time into tracks.
    mode 'M': Hungarian on distance from the position extrapolated by the last full displacement
              (distance only), gates 1000 km from the extrapolation and 1500 km from the last fix.
    mode 'V': Hodges-style: Hungarian on the triplet cost 0.5*(1-cos turn)/2 + 0.5*(1-2*sqrt(s*s')/(s+s')),
              (direction change and speed change, equal weights; 1-fix tracks cost 0.5*d/1000 km), gate 1000 km/6 h.
    No gap bridging. Returns array of track ids (-1 for tracks shorter than min_len) aligned with det."""
    import pandas as pd
    df = pd.DataFrame(det, columns=["time", "lat", "lon", "val"])
    tid = np.full(len(df), -1)
    active = {}               # id -> list of last two row indices
    nxt = 0
    for t, g in df.groupby("time", sort=True):
        idx = g.index.values
        la, lo = g.lat.values, g.lon.values
        ks = list(active)
        cost = np.full((len(ks), len(idx)), 1e9)
        for a, k in enumerate(ks):
            rows = active[k]
            i = rows[-1]
            pla, plo = df.lat[i], df.lon[i]
            d0 = dist(pla, plo, la, lo)
            if len(rows) >= 2 and mode == "M":
                ip = rows[-2]
                dla = pla - df.lat[ip]
                dlo = ((plo - df.lon[ip] + 180) % 360) - 180
                ela, elo = pla + dla, (plo + dlo) % 360
            else:
                ela, elo = pla, plo
            de = dist(ela, elo, la, lo)
            ok = (de < 1000) & (d0 < 1500) if mode == "M" else (d0 < 1000)
            if mode == "M":
                cost[a, ok] = de[ok]
            else:
                if len(rows) >= 2:
                    ip = rows[-2]
                    ux, uy = _xy(df.lat[ip], df.lon[ip], pla, plo)       # previous step
                    vx, vy = _xy(pla, plo, la, lo)                       # candidate step
                    su, sv = np.hypot(ux, uy), np.hypot(vx, vy)
                    cosang = (ux * vx + uy * vy) / np.maximum(su * sv, 1e-6)
                    phi = 0.5 * (1 - cosang) / 2 + 0.5 * (1 - 2 * np.sqrt(su * sv) / np.maximum(su + sv, 1e-6))
                    cost[a, ok] = phi[ok]
                else:
                    cost[a, ok] = (0.5 * d0 / 1000.0)[ok]
        used_j = set()
        new_active = {}
        if len(ks):
            ra, cj = linear_sum_assignment(cost)
            for a, j in zip(ra, cj):
                if cost[a, j] < 1e8:
                    k = ks[a]
                    tid[idx[j]] = k
                    new_active[k] = (active[k] + [idx[j]])[-2:]
                    used_j.add(j)
        for j in range(len(idx)):
            if j not in used_j:
                tid[idx[j]] = nxt
                new_active[nxt] = [idx[j]]
                nxt += 1
        active = new_active
    df["track"] = tid
    n = df.groupby("track").size()
    df.loc[df.track.map(n) < min_len, "track"] = -1
    return df
