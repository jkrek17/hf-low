#!/usr/bin/env python3
"""Context overlays for the Outlook pages, from one GFS 00 UTC run: lows, 500 hPa trough axes, 250 hPa jet cores, and the Coupled flag per low.

None of these carries a probability. They are what a forecaster would look for next to the shading, drawn from the same GFS run.

  lows     local minima of GFS mean sea level pressure (1 degree box mean, a minimum over 3 degrees each way, at or below 992 hPa, 20-75N),
           merged within 5 degrees; blocks with terrain of 500 m or more (Greenland ice, Rockies, Himalaya) are skipped because pressure there is reduced to sea level, the hPa shown is the deepest 0.25 degree value within 1 degree. Display only: the shading never uses them.
  coupled  per low, the two environment terms of research/era5/intensity/env.py (ERA5 proxy, pipeline A) computed from the GFS fields on the
           same 1.5 degree grid: mean 500 hPa absolute vorticity advection within 500 km >= 0.372 (1e-9 s-2) AND mean 850-500 hPa Eady growth
           rate within 500 km >= 0.760 (per day). Both cuts are the upper thirds of the pipeline A ERA5 fixes of 2004-05 to 2025-26 (research/era5/bomb_env/results/
           diagrams.txt), picked post hoc; the same cuts are applied to GFS here, which was not tested for this display.
  jet      1.5 degree cells where the 250 hPa wind speed is at least 110 kt, stored as runs [lat, first lon, cells east] (20-75N, basin views only).
  trough   500 hPa trough axes: along each latitude row, a minimum of the zonal-anomaly height (smoothed) at least 60 m below the mean of the
           heights 12 degrees east and west; minima in neighbouring rows (within 4.5 degrees of longitude) are linked into lines of at least
           6 rows, 25-70N. The research trough finder often picked the western side of a storm's own upper low rather than a separate trough
           (research/era5/hf_jet_trough/README.md), and this one has the same weakness. Context only.

Valid at forecast hour 24, 48, 72, 96 (the end of forecast days 1 to 4). The environment algebra is copied from research/era5/intensity/env.py
(coarsen, derived) so this runs without the research tree. Standard library plus numpy."""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import grib2_t53 as G
import outlook_live as O

LEADS_OV = (24, 48, 72, 96)
RECS = [("PRMSL", "mean sea level", "prmsl"), ("UGRD", "250 mb", "u250"), ("VGRD", "250 mb", "v250"),
        ("UGRD", "500 mb", "u500"), ("VGRD", "500 mb", "v500"), ("UGRD", "850 mb", "u850"), ("VGRD", "850 mb", "v850"),
        ("TMP", "500 mb", "t500"), ("TMP", "850 mb", "t850"), ("HGT", "500 mb", "z500"), ("HGT", "850 mb", "z850")]
GRAV = 9.80665; RE = 6.371e6; OMEGA = 7.292e-5; KT = 1 / 0.514444
VADV_CUT, EADY_CUT = 0.372, 0.760
LOW_HPA, JET_KT, TROUGH_DEPTH_M = 992.0, 110.0, 60.0

TERRAIN = np.unpackbits(np.load(os.path.join(os.path.dirname(HERE), 'data', 'outlook_live', 'high_terrain_1deg.npz'))['m'], axis=1)[:, :360].astype(bool)   # 1 degree blocks (lat 89.5..-89.5) with GFS terrain of 500 m or more

LAT = -90 + 1.5 * np.arange(121); LON = 1.5 * np.arange(240)
R0, R1 = 70, 116                                          # 15N .. 82.5N on the 1.5 degree grid
LATS = LAT[R0:R1]
LAT2, LON2 = np.meshgrid(LATS, LON, indexing="ij")

_W7 = np.array([0.5, 1, 1, 1, 1, 1, 0.5]) / 6.0


def coarsen(a):
    """0.25 degree [721, 1440] (lat 90..-90) -> 1.5 degree [121, 240] (lat -90..90); the cell average of intensity/env.py."""
    a = a.astype(np.float64)
    nan = np.isnan(a); w = (~nan).astype(np.float64); a = np.where(nan, 0.0, a)

    def smooth(x):
        xl = np.concatenate([x[:, -3:], x, x[:, :3]], axis=1)
        y = sum(_W7[k] * xl[:, k:k + 1440] for k in range(7))[:, ::6]
        yl = np.concatenate([y[3:0:-1], y, y[-2:-5:-1]], axis=0)
        return sum(_W7[k] * yl[k:k + 721] for k in range(7))[::6]
    num, den = smooth(a), smooth(w)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0.5 * smooth(np.ones_like(w)), num / den, np.nan)
    return out[::-1]


def ddx(f):
    dl = np.radians(1.5)
    return (np.roll(f, -1, -1) - np.roll(f, 1, -1)) / (2 * dl * RE * np.cos(np.radians(LAT2)))


def ddy(f):
    d = np.full_like(f, np.nan)
    d[1:-1] = (f[2:] - f[:-2]) / (2 * np.radians(1.5) * RE)
    return d


def derived(F):
    """vadv500 (1e-9 s-2), eady (day-1) and ws250 (kt); the same algebra as env.derived()."""
    cosl = np.cos(np.radians(LAT2)); fcor = 2 * OMEGA * np.sin(np.radians(LAT2))
    zeta = ddx(F["v500"]) - ddy(F["u500"] * cosl) / cosl + fcor
    vadv = -(F["u500"] * ddx(zeta) + F["v500"] * ddy(zeta)) * 1e9
    th850 = F["t850"] * (1000 / 850) ** 0.2857; th500 = F["t500"] * (1000 / 500) ** 0.2857
    dz = F["z500"] - F["z850"]
    n2 = GRAV / (0.5 * (th850 + th500)) * (th500 - th850) / dz
    shear = np.hypot(F["u500"] - F["u850"], F["v500"] - F["v850"]) / dz
    with np.errstate(invalid="ignore"):
        eady = 0.3098 * np.abs(fcor) * shear / np.sqrt(np.where(n2 > 1e-6, n2, 1e-6)) * 86400
    return {"vadv500": vadv, "eady": eady, "ws250": np.hypot(F["u250"], F["v250"]) * KT}


def wmean(x, m):
    m = m & np.isfinite(x)
    if not m.any():
        return np.nan
    w = np.cos(np.radians(LAT2[m]))
    return float(np.sum(x[m] * w) / np.sum(w))


def great_circle_km(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lo2 - lo1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def pull(day, lead):
    """The eleven raw 721 x 1440 fields of one lead, as {name: array}; raises O.Missing if any is absent."""
    u = f"{O.BUCKET}/gfs.{day}/00/atmos/gfs.t00z.pgrb2.0p25.f{lead:03d}"
    idx = O.get(u + ".idx"); rows = []
    for line in idx.decode().splitlines():
        p = line.split(":")
        if len(p) >= 6 and p[1].isdigit():
            rows.append((int(p[1]), p[3], p[4]))
    out = {}
    for var, level, name in RECS:
        for i, (off, v, l) in enumerate(rows):
            if v == var and l == level:
                blob = O.get(u, f"bytes={off}-{rows[i + 1][0] - 1}")
                try:
                    out[name] = G.decode_message(blob)
                except (AssertionError, O.struct.error, IndexError, KeyError):
                    out[name] = G.decode_message(O.walk(u, i))
                break
        else:
            raise O.Missing(f"{u} no {var} {level}")
    return out


def find_lows(prmsl):
    """[(lat, lon, hPa)] sorted deepest first; lon in 0..360."""
    p = prmsl[:720].reshape(180, 4, 360, 4).mean(axis=(1, 3)) / 100.0         # 1 degree box means, lat 89.5 .. -89.5
    lat = 89.5 - np.arange(180); lon = 0.5 + np.arange(360)
    r0, r1 = int(np.argmax(lat <= 75)), int(np.argmax(lat < 20))
    pad = np.pad(p, ((3, 3), (0, 0)), constant_values=np.inf); pad = np.concatenate([pad[:, -3:], pad, pad[:, :3]], axis=1)
    is_min = np.ones(p.shape, bool)
    for dy in range(-3, 4):
        for dx in range(-3, 4):
            if dy or dx:
                is_min &= p <= pad[3 + dy:3 + dy + 180, 3 + dx:3 + dx + 360]
    cand = [(p[i, j], i, j) for i in range(r0, r1) for j in range(360) if is_min[i, j] and p[i, j] <= LOW_HPA and not TERRAIN[i, j]]
    cand.sort()
    keep = []
    for v, i, j in cand:
        if any(abs(i - k[1]) <= 5 and min(abs(j - k[2]), 360 - abs(j - k[2])) <= 5 / max(0.2, np.cos(np.radians(lat[i]))) for k in keep):
            continue
        keep.append((v, i, j))
    out = []
    for v, i, j in keep:
        la0 = lat[i]; lo0 = lon[j]
        ri = np.arange(max(0, int(round((90 - la0 - 1) / 0.25))), min(721, int(round((90 - la0 + 1) / 0.25)) + 1))
        ci = np.arange(int(round((lo0 - 1) / 0.25)), int(round((lo0 + 1) / 0.25)) + 1) % 1440
        sub = prmsl[np.ix_(ri, ci)]
        out.append((float(la0), float(lo0), float(sub.min() / 100.0)))
    return out


def trough_lines(z500):
    """Trough axes from the 1.5 degree 500 hPa height (gpm) [lat 15..82.5N, lon]: list of lines of [lat, lon] points."""
    zp = z500 - np.nanmean(z500, axis=1, keepdims=True)
    k = np.array([1, 2, 1]) / 4.0
    pts = {}
    for r in range(zp.shape[0]):
        lat = LATS[r]
        if lat < 25 or lat > 70:
            continue
        row = zp[r]
        for _ in range(2):
            row = np.convolve(np.concatenate([row[-1:], row, row[:1]]), k, mode="valid")
        ms = []
        for j in range(240):
            if row[j] < row[(j - 1) % 240] and row[j] <= row[(j + 1) % 240]:
                edge = 0.5 * (row[(j - 8) % 240] + row[(j + 8) % 240])
                if edge - row[j] >= TROUGH_DEPTH_M:
                    ms.append(j)
        pts[r] = ms
    lines, used = [], set()
    for r in sorted(pts):
        for j in pts[r]:
            if (r, j) in used:
                continue
            line = [(r, j)]; used.add((r, j))
            while True:
                rr, jj = line[-1]
                nxt = [j2 for j2 in pts.get(rr + 1, []) if (rr + 1, j2) not in used and min(abs(j2 - jj), 240 - abs(j2 - jj)) <= 3]
                if not nxt:
                    break
                j2 = min(nxt, key=lambda x: min(abs(x - jj), 240 - abs(x - jj)))
                line.append((rr + 1, j2)); used.add((rr + 1, j2))
            if len(line) >= 6:
                lines.append([[float(LATS[a]), float(LON[b])] for a, b in line])
    return lines


def overlays_for(raw):
    """raw: {name: 721 x 1440 array} for one lead -> dict(lows=[[lat, lon, hPa, coupled 0/1]], jet=[[lat, first lon, n cells]], trough=[[[lat, lon], ...]])."""
    F = {k: coarsen(raw[k])[R0:R1] for k in ("u250", "v250", "u500", "v500", "u850", "v850", "t500", "t850", "z500", "z850")}
    D = derived(F)
    lows = []
    for la, lo, hpa in find_lows(raw["prmsl"]):
        if not in_views(lo):
            continue
        d = great_circle_km(LAT2, LON2, la, lo); m5 = d <= 500
        va, ea = wmean(D["vadv500"], m5), wmean(D["eady"], m5)
        lows.append([round(la, 1), round(lo, 1), round(hpa), int(bool(va >= VADV_CUT and ea >= EADY_CUT))])
    jet = []
    for r in range(len(LATS)):
        if not 20 <= LATS[r] <= 75:
            continue
        hot = (D["ws250"][r] >= JET_KT) & np.array([in_views(x) for x in LON])
        j = 0
        while j < 240:
            if hot[j]:
                k = j
                while k + 1 < 240 and hot[k + 1]:
                    k += 1
                jet.append([float(LATS[r]), float(LON[j]), k - j + 1]); j = k + 1
            else:
                j += 1
    tr = [[[round(a, 1), round(b, 1)] for a, b in ln] for ln in trough_lines(F["z500"]) if in_views(ln[len(ln) // 2][1])]
    return {"lows": lows, "jet": jet, "trough": tr}


def in_views(lon):
    """Atlantic view 100W-10E (260..360 and 0..10) or Pacific view 120E-120W (120..240)."""
    return lon >= 260 or lon <= 10.5 or 119.5 <= lon <= 240.5


def for_cycle(day, workers=4):
    """[overlays for forecast day 1..4 or None]; never raises on a missing or undecodable lead (the page then shows no overlay for that day)."""
    from concurrent.futures import ThreadPoolExecutor

    def one(lead):
        try:
            return overlays_for(pull(day, lead))
        except O.Missing:
            return None
        except Exception as e:                                   # overlays are context; a bad field must not stop the shading
            print(f"::warning::overlay lead {lead} failed: {e!r}")
            return None
    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(one, LEADS_OV))
