"""Jet-streak (J2) and upstream-trough (T2) features at one fix (see PREREGISTRATION.md).

fix_features(F, lat, lon, heading) -> dict. F is the per-time dict from time_fields(): grid fields [121, 240]
(lat -90..90 south to north, lon 0..358.5): u250 v250 u300 v300 and zp (Z500 minus its zonal mean, m), plus the
zonal-mean 250 hPa speed per latitude row, vbar (kt).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_env_composites"))
import numpy as np
from scipy import ndimage
from common import LAT, LON, RE, BOX_N, BOX_D, BX, BY, BR, destination, sample, G

KT = 1.943844
N0 = BOX_N // 2
THETA = np.arctan2(BY, BX)
SIN, COS = np.sin(THETA), np.cos(THETA)
BRG = np.degrees(np.arctan2(BX, BY)) % 360                    # bearing from north, clockwise
RING = np.rint(BR / BOX_D).astype(int)
KERN = np.outer([1, 2, 1], [1, 2, 1]) / 16.0
J2 = ["nojet", "vmaxp", "Lhalf", "RE", "LE", "RX", "LX", "s", "n", "dvds", "div300max", "div300dist"]
T2 = ["notrough", "tdepth", "tdist", "tbear_cos", "tbear_sin", "tamp", "ttilt", "tphase", "nohead"]
EXTRA = ["vmax2500", "div300mean1000"]
COLS = J2 + T2 + EXTRA


def time_fields(U250, V250, U300, V300, Z500):
    """Grid fields [121, 240] from chunk arrays already sliced at one time ([lon, lat] order)."""
    u250, v250, u300, v300 = (a.T.astype(np.float64) for a in (U250, V250, U300, V300))
    z = Z500.T.astype(np.float64) / G
    return dict(u250=u250, v250=v250, u300=u300, v300=v300, zp=z - z.mean(1, keepdims=True),
                vbar=np.hypot(u250, v250).mean(1) * KT)


def smooth(a):
    return ndimage.convolve(a, KERN, mode="nearest")


def remove_vortex(u, v, rm=1500):
    """Subtract the azimuthal-mean tangential wind (100 km rings; weight 1 inside rm km, 0 beyond rm + 500 km). rm = 0: none."""
    if rm <= 0:
        return u, v
    taper = np.clip((rm + 500.0 - BR) / 500.0, 0.0, 1.0)
    vt = u * -SIN + v * COS
    vt0 = np.where(np.isfinite(vt), vt, 0.0)
    cnt = np.bincount(RING.ravel(), weights=np.isfinite(vt).ravel().astype(float), minlength=RING.max() + 1)
    sm = np.bincount(RING.ravel(), weights=vt0.ravel(), minlength=RING.max() + 1)
    mean = np.where(cnt > 0, sm / np.maximum(cnt, 1), 0.0)
    mean[0] = 0.0
    corr = mean[RING] * taper
    return u - corr * -SIN, v - corr * COS


def box(F, name, lat, lon, la, lo):
    return sample(F[name], la, lo)


def div1e5(u, v):
    return (np.gradient(u, BOX_D * 1e3, axis=1) + np.gradient(v, BOX_D * 1e3, axis=0)) * 1e5


def idx(x_km, y_km):
    return [y_km / BOX_D + N0, x_km / BOX_D + N0]


def fix_features(F, lat, lon, heading, dbg=None, rm=1500, excl=750):
    """dbg: optional dict that receives the box arrays used (for the Q3 drawings); does not change the result."""
    out = {c: np.nan for c in COLS}
    la, lo = destination(lat, lon, 90.0 + np.degrees(np.arctan2(-BY, BX)), BR)
    ok = np.isfinite(sample(F["u250"], la, lo))
    S = {k: np.where(ok, sample(F[k], la, lo), 0.0) for k in ("u250", "v250", "u300", "v300", "zp")}
    vbar_box = np.interp(la, LAT, F["vbar"])
    if dbg is not None:
        dbg.update(ok=ok)
    # ---- jet streak ----
    u, v = remove_vortex(S["u250"], S["v250"], rm)
    u, v = smooth(u), smooth(v)
    V = np.hypot(u, v) * KT
    Vp = V - vbar_box
    out["vmax2500"] = float(V[BR <= 2500].max())
    mask = (V >= 60) & (Vp >= 20) & ok
    near = mask & (BR <= 2500)
    if not near.any():
        out["nojet"] = 1.0
        for c in ("RE", "LE", "RX", "LX"):
            out[c] = 0.0
    else:
        out["nojet"] = 0.0
        lab, _ = ndimage.label(mask, structure=np.ones((3, 3)))
        i0 = np.unravel_index(np.argmax(np.where(near, V, -1)), V.shape)
        comp = lab == lab[i0]
        j = np.unravel_index(np.argmax(np.where(comp, V, -1)), V.shape)
        xm, ym = BX[j], BY[j]
        if dbg is not None:
            dbg.update(V=V, comp=comp, jmax=(xm, ym), a=a)
        out["vmaxp"] = float(Vp[comp].max())
        w = V[comp] - 60.0
        x, y = BX[comp], BY[comp]
        mu, mv = (u[comp] * w).sum(), (v[comp] * w).sum()
        if comp.sum() >= 3 and w.sum() > 0:
            cx, cy = (w * x).sum() / w.sum(), (w * y).sum() / w.sum()
            C = np.cov(np.vstack([x - cx, y - cy]), aweights=w + 1e-9)
            ev, evec = np.linalg.eigh(C)
            a = evec[:, -1]
        else:
            a = np.array([mu, mv]) / (np.hypot(mu, mv) + 1e-9)
        if a[0] * mu + a[1] * mv < 0:
            a = -a
        vmax = V[j]
        half = comp & (V >= 0.5 * (60 + vmax))
        sp = x[:0]
        proj = BX[half] * a[0] + BY[half] * a[1]
        out["Lhalf"] = float((proj.max() - proj.min()) / 2 / 1000.0)
        s = -(xm * a[0] + ym * a[1])
        n = xm * a[1] - ym * a[0]
        out["s"], out["n"] = float(s / 1000), float(n / 1000)
        inq = abs(s) <= 2000 and abs(n) <= 1500
        for c, cond in (("RE", s < 0 and n < 0), ("LE", s < 0 and n > 0), ("RX", s > 0 and n < 0), ("LX", s > 0 and n > 0)):
            out[c] = float(inq and cond)
        px = lambda sp_: ndimage.map_coordinates(V, idx(np.array([xm + sp_ * a[0]]), np.array([ym + sp_ * a[1]])), order=1, mode="nearest")[0]
        out["dvds"] = float((px(s + 300) - px(s - 300)) / 0.6)
    # ---- 300 hPa divergence of the vortex-removed wind ----
    u3, v3 = remove_vortex(S["u300"], S["v300"], rm)
    d = div1e5(smooth(u3), smooth(v3))
    d = np.where(ok, d, np.nan)
    m15 = BR <= 1500
    if np.isfinite(d[m15]).any():
        k = np.nanargmax(np.where(m15, d, -np.inf))
        out["div300max"] = float(d.ravel()[k])
        out["div300dist"] = float(BR.ravel()[k] / 1000)
    m10 = BR <= 1000
    out["div300mean1000"] = float(np.nanmean(d[m10]))
    # ---- upstream trough ----
    if not np.isfinite(heading):
        out["nohead"], out["notrough"] = 1.0, 0.0
        return out
    out["nohead"] = 0.0
    zs = smooth(S["zp"])
    opp = (heading + 180.0) % 360
    dif = np.abs((BRG - opp + 180.0) % 360 - 180.0)
    sect = (BR >= excl) & (BR <= 3000) & (dif <= 75) & ok
    lmin = zs <= ndimage.minimum_filter(zs, size=3, mode="nearest")
    cand = sect & lmin
    if not cand.any():
        out["notrough"] = 1.0
        return out
    out["notrough"] = 0.0
    k = np.argmin(np.where(cand, zs, np.inf))
    ij = np.unravel_index(k, zs.shape)
    xt, yt, zt = BX[ij], BY[ij], zs[ij]
    if dbg is not None:
        dbg.update(zs=zs, trough=(xt, yt), sect=sect)
    out["tdepth"], out["tdist"] = float(-zt), float(BR[ij] / 1000)
    rel = np.radians(BRG[ij] - heading)
    out["tbear_cos"], out["tbear_sin"] = float(np.cos(rel)), float(np.sin(rel))
    east = (BX > xt) & (BR <= 3000) & ok
    out["tamp"] = float(zs[east].max() - zt) if east.any() else np.nan
    # tilt from the grid rows around the trough
    lat_t, lon_t = destination(lat, lon, BRG[ij], BR[ij])
    lon_t = float(lon_t)
    r0 = int(round((lat_t + 90) / 1.5))
    rows, lons = [], []
    zp = F["zp"]
    for r in range(max(r0 - 5, 1), min(r0 + 5, 119) + 1):
        c0 = int(round(lon_t / 1.5))
        cols = (np.arange(c0 - 20, c0 + 21)) % 240
        seg = zp[r, cols]
        m = int(np.argmin(seg))
        if m in (0, len(seg) - 1):
            continue
        rows.append(LAT[r])
        lons.append((m - 20) * 1.5 + (c0 * 1.5 - lon_t))
    if len(rows) >= 5:
        b = np.polyfit(rows, lons, 1)[0]
        cphi = np.cos(np.radians(lat_t))
        out["ttilt"] = float(np.degrees(np.arctan(b * cphi)))
        dvec = np.array([b * cphi, 1.0])
        dvec /= np.hypot(*dvec)
        e = np.array([dvec[1], -dvec[0]])
        out["tphase"] = float((-xt * e[0] - yt * e[1]) / 1000)
    return out
