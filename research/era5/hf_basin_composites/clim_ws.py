"""WeatherBench2 1990-2019 climatological MEAN SPEED at 250 hPa (the mean of |V|, not |mean V|) -> $ERA5_WORK/clim_ws250.npz [hour 4, dayofyear 366, lat 121, lon 240] m/s.
Needed because |V| - |V_clim| (V_clim = climatological mean vector) is biased high everywhere (speed of the mean wind is below the mean speed).
usage: ERA5_WORK=DIR clim_ws.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "hf_env_composites"))
import numpy as np
from common import *
if __name__ == "__main__":
    li = LEVELS.index(250); lc, lp = li // 3, li % 3
    arr, nbytes = np.zeros((4, 366, 121, 240), np.float32), 0
    for dc in range(13):
        a, nb = decode(CLIM, "wind_speed", f"0.{dc}.{lc}.0.0")
        nbytes += nb
        n = min(30, 366 - dc * 30)
        arr[:, dc * 30:dc * 30 + n] = a[:, :n, lp].transpose(0, 1, 3, 2)
    np.savez(os.path.join(WORK, "clim_ws250.npz"), ws250=arr)
    open(os.path.join(WORK, "bytes_climws.txt"), "w").write(str(nbytes))
    print("MB", nbytes / 1e6, "mean m/s", float(arr.mean()), "at 45N 30W Jan 1 12Z:", float(arr[2, 0, int((45 + 90) / 1.5), int(330 / 1.5)]))
