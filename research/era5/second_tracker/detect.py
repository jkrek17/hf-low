"""Run detection for tracker M (MSLP, 1979-2022, one file per calendar year) or V (850 hPa vorticity, Oct-Apr plus
4 days margin, 2004-05..2021-22, one file per season). Resumable; writes WORK/{M,V}/<unit>.csv (time,lat,lon,val).
usage: detect.py M|V FIRST LAST [PROCS]       (M: calendar years; V: season start years)
"""
import os, sys, datetime as dt
import numpy as np, pandas as pd
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor
import common as C

WORK = os.environ.get("ERA5_WORK", os.path.join(os.path.dirname(os.path.abspath(__file__)), "work"))
VAR = {"M": "mean_sea_level_pressure", "V": "vorticity"}


def window(tracker, u):
    if tracker == "M":
        return dt.datetime(u, 1, 1), dt.datetime(u + 1, 1, 1)
    return dt.datetime(u, 9, 26), dt.datetime(u + 1, 5, 4)      # 1 Oct - 30 Apr with 4 days margin each side


def unit(args):
    tracker, u = args
    path = f"{WORK}/{tracker}/{u}.csv"
    if os.path.exists(path):
        return u
    t0, t1 = window(tracker, u)
    c0, _ = C.chunk_of(t0)
    c1, _ = C.chunk_of(t1 - dt.timedelta(hours=6))
    cs = list(range(c0, c1 + 1))
    det = C.detect_m if tracker == "M" else C.detect_v
    rows = []
    with ThreadPoolExecutor(6) as ex:
        for c, arr in zip(cs, ex.map(lambda c: C.fetch(VAR[tracker], c), cs)):
            for k in range(8):
                t = C.stamp(c, k)
                tt = dt.datetime.strptime(str(t), "%Y%m%d%H")
                if not (t0 <= tt < t1):
                    continue
                for la, lo, v in det(arr[k]):
                    rows.append((t, la, lo, v))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(rows, columns=["time", "lat", "lon", "val"]).to_csv(path + ".tmp", index=False)
    os.replace(path + ".tmp", path)
    return u


if __name__ == "__main__":
    tr, a, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    n = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    with Pool(n) as p:
        for u in p.imap_unordered(unit, [(tr, y) for y in range(a, b + 1)]):
            print("done", tr, u, flush=True)
