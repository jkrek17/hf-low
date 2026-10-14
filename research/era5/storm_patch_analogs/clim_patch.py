"""WeatherBench2 1990-2019 climatology for the patch channels -> $ERA5_WORK/clim_patch.npz (z500 m, msl hPa, u250, v250;
each [hour 4, dayofyear 366, lat 121, lon 240]). Same read logic as hf_env_composites/clim.py and hf_basin_composites/clim_msl.py."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common_patch import *

NEED = {"u250": ("u_component_of_wind", 250), "v250": ("v_component_of_wind", 250), "z500": ("geopotential", 500)}
if __name__ == "__main__":
    os.makedirs(WORK, exist_ok=True)
    doy, hr = read_small(CLIM, "dayofyear"), read_small(CLIM, "hour")
    out, nbytes = {}, 0
    for name, (var, lev) in NEED.items():
        li = LEVELS.index(lev); lc, lp = li // 3, li % 3
        arr = np.zeros((4, 366, 121, 240), np.float32)
        for dc in range(13):
            a, nb = decode(CLIM, var, f"0.{dc}.{lc}.0.0"); nbytes += nb
            n = min(30, 366 - dc * 30)
            arr[:, dc * 30:dc * 30 + n] = a[:, :n, lp].transpose(0, 1, 3, 2) / (G if var == "geopotential" else 1)
        out[name] = arr; print(name, round(nbytes / 1e9, 2), "GB", float(arr.mean()), flush=True)
    arr = np.zeros((4, 366, 121, 240), np.float32)
    for dc in range(13):
        a, nb = decode(CLIM, "mean_sea_level_pressure", f"0.{dc}.0.0"); nbytes += nb
        n = min(30, 366 - dc * 30)
        arr[:, dc * 30:dc * 30 + n] = a[:, :n].transpose(0, 1, 3, 2) / 100.0
    out["msl"] = arr
    np.savez(os.path.join(WORK, "clim_patch.npz"), doy=doy, hour=hr, **out)
    open(os.path.join(WORK, "bytes_clim.txt"), "w").write(str(nbytes))
    print("total GB", nbytes / 1e9, "msl mean hPa", float(arr.mean()))
