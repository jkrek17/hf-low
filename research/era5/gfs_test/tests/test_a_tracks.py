"""Synthetic-field test of the GFS tracker: python tests/test_a_tracks.py (no network)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import a_tracks as A


def field(lat0, lon0, depth=30.0):
    la = (90 - 0.25 * np.arange(721))[:, None]
    lo = (0.25 * np.arange(1440))[None, :]
    dlo = (lo - lon0 + 180) % 360 - 180
    d2 = ((la - lat0) ** 2 + (dlo * np.cos(np.radians(lat0))) ** 2)     # deg^2
    return 1015.0 - depth * np.exp(-d2 / (2 * 4.0 ** 2)), la, lo


def run(lat0, lon0):
    msl, la, lo = field(lat0, lon0)
    gust = np.full((721, 1440), 5.0)
    patch = (np.abs(la - (lat0 - 3)) < 1) & (np.abs((lo - lon0 + 180) % 360 - 180) < 1)   # ~400 km from the low
    gust[patch] = 40.0
    far = (np.abs(la - lat0) < 1) & (np.abs((lo - (lon0 + 30) + 180) % 360 - 180) < 1)      # far outside 800 km
    gust[far] = 60.0
    ocean = np.ones((721, 1440), bool)
    out = A.lows_in(msl, gust, ocean)
    assert len(out) == 1, out
    lat, lon, p, g500, g800, g1200 = out[0]
    assert abs(lat - lat0) <= 0.5 and abs(lon - lon0) <= 0.5, (lat, lon)
    assert abs(p - 985.0) < 0.2, p
    assert abs(g800 - 40.0 / 0.514444) < 0.2 and abs(g1200 - g800) < 1e-6, (g800, g1200)
    assert g500 > 0
    # same patch over land is ignored
    out2 = A.lows_in(msl, gust, ~patch & ~far)
    assert out2[0][4] == 0.0, out2


run(50.0, 330.0)
run(45.0, 359.75)      # a low on the longitude seam must be found once, at the seam
run(55.0, 0.0)
print("ok")
