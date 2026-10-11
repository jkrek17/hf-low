"""Ascent (W) and ascent-coupling (C) features at one fix (see PREREGISTRATION.md).

fix_features_w(F, lat, lon, heading) -> dict. F = time_fields_w(): A500, A700 grids [121, 240] (lat -90..90, lon 0..358.5),
A = -omega (Pa/s, positive = rising). Box +-5000 km about the low, 100 km spacing, bilinear, great-circle placement.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_env_composites"))
import numpy as np
from scipy import ndimage
from common import destination, sample

NB, D = 101, 100.0
N0 = NB // 2
_ax = (np.arange(NB) - N0) * D
BX, BY = np.meshgrid(_ax, _ax)                      # x east, y north (km); row = y
BR = np.hypot(BX, BY)
BRG = np.degrees(np.arctan2(BX, BY)) % 360          # bearing from north, clockwise
KERN = np.outer([1, 2, 1], [1, 2, 1]) / 16.0
W = ["a500_1000", "a700_1000", "a500_500", "a700_500", "aarea500", "adx", "ady", "afwd500", "a500_x750"]
C = ["nasc", "asc2", "asc2dist", "asc2cos", "asc2sin", "couple_w"]
CV = [c + v for v in ("_o", "_m") for c in C]
COLS_W = W + C + CV + ["nohead_w"]
COUPLE_KM = 2500.0


def time_fields_w(W500, W700):
    """Chunk arrays sliced at one time and level ([lon, lat] order, Pa/s) -> A = -omega on [121, 240]."""
    return dict(a500=-W500.T.astype(np.float64), a700=-W700.T.astype(np.float64))


def smooth(a):
    return ndimage.convolve(a, KERN, mode="nearest")


def fix_features_w(F, lat, lon, heading):
    out = {c: np.nan for c in COLS_W}
    la, lo = destination(lat, lon, BRG, BR)
    A5, A7 = sample(F["a500"], la, lo), sample(F["a700"], la, lo)
    ok = np.isfinite(A5) & np.isfinite(A7)
    A5, A7 = smooth(np.where(ok, A5, 0.0)), smooth(np.where(ok, A7, 0.0))
    A5, A7 = np.where(ok, A5, -np.inf), np.where(ok, A7, -np.inf)
    mx = lambda A, m: float(np.max(np.where(m, A, -np.inf))) if (m & ok).any() else np.nan
    mn = lambda A, m: float(A[m & ok].mean()) if (m & ok).any() else np.nan
    out["a500_1000"], out["a700_1000"] = mx(A5, BR <= 1000), mx(A7, BR <= 1000)
    out["a500_500"], out["a700_500"] = mn(A5, BR <= 500), mn(A7, BR <= 500)
    out["aarea500"] = float(((A5 >= 0.3) & (BR <= 3000) & ok).sum() * D * D / 1e6)
    out["a500_x750"] = mx(A5, (BR > 750) & (BR <= 2500))
    if np.isfinite(heading):
        out["nohead_w"] = 0.0
        h = np.radians(heading)
        e = np.array([np.sin(h), np.cos(h)])               # forward (x east, y north)
        lf = np.array([-np.cos(h), np.sin(h)])             # left of motion
        near = (BR <= 2500) & ok
        if near.any():
            k = np.argmax(np.where(near, A5, -np.inf))
            x, y = BX.ravel()[k], BY.ravel()[k]
            out["adx"], out["ady"] = float((x * e[0] + y * e[1]) / 1000), float((x * lf[0] + y * lf[1]) / 1000)
        al = BX * e[0] + BY * e[1]
        m2 = (BR <= 2000) & ok & (BR > 0)
        fw, rr = m2 & (al > 0), m2 & (al < 0)
        if fw.any() and rr.any():
            out["afwd500"] = float(A5[fw].mean() - A5[rr].mean())
    else:
        out["nohead_w"] = 1.0
    # ---- ascent centres (greedy, >= sep km apart); primary 0.8 Pa/s, 1,500 km, within 3,000 km; variants _o (0.4, 1,000, 4,000: as first written) and _m (0.6, 1,500, 3,000)
    for suf, thr, sep, rmax in (("", 0.8, 1500.0, 3000.0), ("_o", 0.4, 1000.0, 4000.0), ("_m", 0.6, 1500.0, 3000.0)):
        cen = centres(A5, ok, thr, sep, rmax)
        out["nasc" + suf] = float(len(cen))
        out["asc2" + suf] = float(cen[1][2]) if len(cen) > 1 else 0.0
        out["couple_w" + suf] = 0.0
        if len(cen) > 1:
            x, y, _ = cen[1]
            out["asc2dist" + suf] = float(np.hypot(x, y) / 1000)
            out["couple_w" + suf] = float(np.hypot(x, y) <= COUPLE_KM)
            if np.isfinite(heading):
                rel = np.radians(np.degrees(np.arctan2(x, y)) - heading)
                out["asc2cos" + suf], out["asc2sin" + suf] = float(np.cos(rel)), float(np.sin(rel))
    return out


def centres(A5, ok, thr, sep, rmax):
    cand = (A5 >= thr) & (A5 >= ndimage.maximum_filter(A5, size=3, mode="nearest")) & (BR <= rmax) & ok
    ys, xs = np.nonzero(cand)
    cen = []
    for i in np.argsort(-A5[ys, xs]):
        x, y = BX[ys[i], xs[i]], BY[ys[i], xs[i]]
        if all(np.hypot(x - cx, y - cy) >= sep for cx, cy, _ in cen):
            cen.append((x, y, A5[ys[i], xs[i]]))
    return cen
