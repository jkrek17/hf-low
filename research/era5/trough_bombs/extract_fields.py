"""Z500, U250 and V250 at 00 and 12 UTC, 5.625 degree grid, Northern Hemisphere 14-87 N, seasons 1979-80..2014-15.

ERA5 proxy, outcome-free. Reuses `../jet_trough/extract_fields.py` (WeatherBench2 64 x 32 reader) with V250
dropped; the WB2 store ends 2023-01-09 so no ARCO splice is needed for these seasons.

usage: extract_fields.py OUT_DIR Y0 Y1            (resumable: one npz per season)
Writes fields_<season>.npz: time (hours since 1900-01-01), z500 (m), u250, v250 (m/s), shape [ntime, 14, 64].
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "jet_trough"))
import extract_fields as E

if __name__ == "__main__":
    out, y0, y1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    os.makedirs(out, exist_ok=True)
    tot = 0
    for s in range(y1, y0 - 1, -1):
        path = f"{out}/fields_{s}.npz"
        if os.path.exists(path):
            continue
        t0 = time.time()
        hrs = E.season_hours(s)
        o, nb = E.wb2_times(hrs)
        tmp = path + ".tmp.npz"
        np.savez_compressed(tmp, time=hrs, **o)
        os.replace(tmp, path)
        tot += nb
        print(f"{s} {len(hrs)} times  {nb / 1e9:.2f} GB  {time.time() - t0:.0f}s  nan z500 {np.isnan(o['z500']).mean():.4f}", flush=True)
    print(f"streamed about {tot / 1e9:.1f} GB this run")
