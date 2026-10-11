"""WeatherBench2 1990-2019 mean sea level pressure climatology (about 0.1 GB) -> $ERA5_WORK/clim_msl.npz [hour 4, dayofyear 366, lat 121, lon 240] hPa.
usage: ERA5_WORK=DIR clim_msl.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "hf_env_composites"))
import numpy as np
from common import *
if __name__ == "__main__":
    arr, nbytes = np.zeros((4, 366, 121, 240), np.float32), 0
    for dc in range(13):
        a, nb = decode(CLIM, "mean_sea_level_pressure", f"0.{dc}.0.0")
        nbytes += nb
        n = min(30, 366 - dc * 30)
        arr[:, dc * 30:dc * 30 + n] = a[:, :n].transpose(0, 1, 3, 2) / 100.0
    np.savez(os.path.join(WORK, "clim_msl.npz"), msl=arr)
    open(os.path.join(WORK, "bytes_climmsl.txt"), "w").write(str(nbytes))
    print("MB", nbytes / 1e6, "mean hPa", float(arr.mean()))
