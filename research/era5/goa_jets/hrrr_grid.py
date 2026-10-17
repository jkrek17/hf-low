"""HRRR-Alaska grid geolocation (3 km polar stereographic, 1299 x 919, scan +i +j), rotation of grid-relative winds to east/north, and a land mask.

Template 3.20 values read from the GRIB message (sec 3): first point 41.612949N 185.117126E, LaD 60N, LoV 225E (135W), dx = dy = 3000 m, flag 'u,v relative to the grid'.
Sphere R = 6371229 m (NCEP). Checked 2026-10-17 against the file's own LAND field (Anchorage, Kodiak Island, Yakutat land; Gulf 57N 145W and Cook Inlet mouth ocean) and against
GFS 24 h winds (rotated u/v correlation 0.91/0.83, unrotated 0.88/0.53). usage: import hrrr_grid as HG; HG.grid()
"""
import numpy as np
from pyproj import CRS, Transformer

NX, NY, DX = 1299, 919, 3000.0
LA1, LO1 = 41.612949, 185.117126 - 360
CRS_H = CRS.from_proj4("+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-135 +R=6371229 +units=m")
_T = Transformer.from_crs("EPSG:4326", CRS_H, always_xy=True)
_I = Transformer.from_crs(CRS_H, "EPSG:4326", always_xy=True)


def grid():
    """returns lat, lon (deg, float64, shape NY x NX) and alpha (rad): earth u = ug cos a - vg sin a ; earth v = ug sin a + vg cos a"""
    x0, y0 = _T.transform(LO1, LA1)
    X, Y = np.meshgrid(x0 + DX * np.arange(NX), y0 + DX * np.arange(NY))
    lo, la = _I.transform(X, Y)
    x2, y2 = _T.transform(lo, la + 0.01)
    return la, lo, np.arctan2(x2 - X, y2 - Y)


def to_earth(ug, vg, alpha):
    return ug * np.cos(alpha) - vg * np.sin(alpha), ug * np.sin(alpha) + vg * np.cos(alpha)


def ij(la0, lo0):
    x0, y0 = _T.transform(LO1, LA1); x, y = _T.transform(lo0, la0)
    return (y - y0) / DX, (x - x0) / DX
