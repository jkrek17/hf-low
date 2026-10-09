"""Pull the WeatherBench2 ERA5 climatology needed for the anomalies (about 1.5 GB) and save it compactly.

usage: clim.py   (writes $ERA5_WORK/clim.npz with u250 v250 u500 v500 z500, each [hour 4, dayofyear 366, lat 121, lon 240])
The climatology chunks are [4 hours, 30 days, 3 levels, lon, lat]; levels 200/250/300 are level chunk 1 and
400/500/600 are level chunk 2.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import *

NEED = {"u250": ("u_component_of_wind", 250), "v250": ("v_component_of_wind", 250),
        "u500": ("u_component_of_wind", 500), "v500": ("v_component_of_wind", 500),
        "z500": ("geopotential", 500)}

if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True)
    doy, hr = read_small(CLIM, "dayofyear"), read_small(CLIM, "hour")
    print("dayofyear", doy[:3], doy[-3:], len(doy), "hour", hr)
    out, nbytes = {}, 0
    for name, (var, lev) in NEED.items():
        li = LEVELS.index(lev)
        lc, lp = li // 3, li % 3
        arr = np.zeros((4, 366, 121, 240), np.float32)
        for dc in range(13):
            a, nb = decode(CLIM, var, f"0.{dc}.{lc}.0.0")
            nbytes += nb
            n = min(30, 366 - dc * 30)
            f = a[:, :n, lp]                                   # [hour, day, lon, lat]
            arr[:, dc * 30:dc * 30 + n] = f.transpose(0, 1, 3, 2) / (G if var == "geopotential" else 1)
        out[name] = arr
        print(name, "read", round(nbytes / 1e9, 2), "GB so far", float(np.nanmean(arr)), flush=True)
    np.savez(os.path.join(WORK, "clim.npz"), doy=doy, hour=hr, **out)
    open(os.path.join(WORK, "bytes_clim.txt"), "w").write(str(nbytes))
    print("total bytes", nbytes)
