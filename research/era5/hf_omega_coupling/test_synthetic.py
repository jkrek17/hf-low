"""Synthetic-field test of omega_features.py: two Gaussian ascent blobs at known places; run before any real field is read.
usage: python3 test_synthetic.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "hf_env_composites"))
import numpy as np
from common import LAT, LON, destination
from omega_features import *

def blob(lat0, lon0, amp, rad_km):
    LA, LO = np.meshgrid(LAT, LON, indexing="ij")
    p, l, p0, l0 = map(np.radians, (LA, LO, lat0, lon0))
    d = 6371 * np.arccos(np.clip(np.sin(p) * np.sin(p0) + np.cos(p) * np.cos(p0) * np.cos(l - l0), -1, 1))
    return amp * np.exp(-(d / rad_km) ** 2)

lat, lon = 50.0, 330.0
la2, lo2 = destination(lat, lon, 90.0, 1800.0)                      # 1,800 km due east
A = blob(lat, lon, 0.9, 400) + blob(float(la2), float(lo2), 0.9, 400)
F = dict(a500=A, a700=0.5 * A)
f = fix_features_w(F, lat, lon, 90.0)
print({k: round(v, 3) for k, v in f.items()})
assert f["nasc"] == 2 and f["nasc_o"] == 2 and abs(f["asc2dist"] - 1.8) < 0.1 and f["couple_w"] == 1
assert abs(f["asc2cos"] - 1) < 0.05 and abs(f["asc2sin"]) < 0.1             # ahead of an eastward-moving low
assert 0.7 < f["a500_1000"] <= 0.95 and abs(f["a700_1000"] / f["a500_1000"] - 0.5) < 0.02
f2 = fix_features_w(F, lat, lon, 270.0)                                     # low moving west: the second centre is behind
assert abs(f2["asc2cos"] + 1) < 0.05
f3 = fix_features_w(dict(a500=blob(lat, lon, 0.9, 400), a700=blob(lat, lon, 0.9, 400)), lat, lon, 90.0)
assert f3["nasc"] == 1 and f3["asc2"] == 0 and f3["couple_w"] == 0 and np.isnan(f3["asc2dist"])
f4 = fix_features_w(F, lat, lon, np.nan)
assert f4["nohead_w"] == 1 and np.isnan(f4["asc2cos"]) and f4["nasc"] == 2
print("synthetic tests passed")
