"""Monthly ERA5 SST box means for the ENSO flavor indices (resumable).

4 samples a month (days 4, 11, 18, 25 at 12 UTC), 1979-01..2026-02, one global
1.37 MB chunk each (about 3 GB). Writes sst_boxes_raw.csv: month, box means (K).
usage: sst_indices.py OUT_CSV
"""
import sys, os
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
from common import get, tindex, LAT, LON

# (lat0, lat1, lon0, lon1) in degrees east, 0..360
BOXES = {
    "nino3": (-5, 5, 210, 270), "nino4": (-5, 5, 160, 210), "nino34": (-5, 5, 190, 240),
    "emi_a": (-10, 10, 165, 220), "emi_b": (-15, -5, 250, 290), "emi_c": (-10 + 0, 20, 125, 145),
}
# EMI box C is 10S-20N
BOXES["emi_c"] = (-10, 20, 125, 145)
W = np.cos(np.radians(LAT))[:, None] * np.ones((1, len(LON)))


def boxmean(sst, b):
    r = (LAT >= b[0]) & (LAT <= b[1]); c = (LON >= b[2]) & (LON <= b[3])
    s = sst[np.ix_(r, c)]; w = W[np.ix_(r, c)]
    ok = np.isfinite(s)
    return float((s[ok] * w[ok]).sum() / w[ok].sum())


def one(args):
    y, m = args
    vals = {k: [] for k in BOXES}
    for d in (4, 11, 18, 25):
        a = get("sea_surface_temperature", tindex(f"{y}-{m:02d}-{d:02d}T12"))
        for k, b in BOXES.items():
            vals[k].append(boxmean(a, b))
    return dict(year=y, month=m, **{k: float(np.mean(v)) for k, v in vals.items()})


if __name__ == "__main__":
    out = sys.argv[1]
    months = [(y, m) for y in range(1979, 2027) for m in range(1, 13) if (y, m) <= (2026, 2)]
    with ThreadPoolExecutor(8) as ex:
        rows = list(ex.map(one, months))
    pd.DataFrame(rows).to_csv(out, index=False)
    print(len(rows), "months")
