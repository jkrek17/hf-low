"""Per-month exceedance counts from the compact gust cache (extract_cache.py). ERA5 PROXY, not observation.

usage: python3 -I monthly_counts.py CACHE_DIR OUT.npz

For each calendar month and each cutoff in CUTOFFS_KT: the number of sampled times (00 and 12 UTC, 2004-05 to 2025-26) at which the
gust was at or above the cutoff, per cell; plus the number of sampled times per month and per (season, month). Cache values are
0.5 kt steps (uint8 = kt * 2), so a cutoff is compared as value >= ceil(cutoff * 2).
Hours per season in a month = count * 12 h / 22 seasons.
"""
import glob, os, sys, datetime
import numpy as np, numcodecs

CUTOFFS_KT = (34.0, 38.0, 50.0, 54.0, 64.0, 71.7, 80.0, 90.0)
T0 = datetime.datetime(1900, 1, 1)
_c = numcodecs.Blosc()


def main(cache, out):
    files = sorted(glob.glob(os.path.join(cache, "gust_*.npz")))
    cnt = np.zeros((12, len(CUTOFFS_KT), 361, 1440), np.uint16)
    ntimes = np.zeros(12, np.int64)
    amax = np.zeros((12, 361, 1440), np.uint8)
    lim = [int(np.ceil(c * 2 - 1e-9)) for c in CUTOFFS_KT]
    for f in files:
        z = np.load(f)
        shp = tuple(z["shape"])
        arr = np.frombuffer(_c.decode(z["data"].tobytes()), np.uint8).reshape(shp)
        for i, h in enumerate(z["hours_since_1900"]):
            m = (T0 + datetime.timedelta(hours=int(h))).month - 1
            ntimes[m] += 1
            g = arr[i]
            amax[m] = np.maximum(amax[m], g)
            for k, l in enumerate(lim):
                cnt[m, k] += g >= l
    np.savez_compressed(out, counts=cnt, n_times=ntimes, cutoffs_kt=np.array(CUTOFFS_KT), max_kt2=amax,
                        note="ERA5 PROXY gust; counts of 00/12 UTC samples per calendar month 2004-05..2025-26; lat 90..0 N rows, lon 0..359.75")
    print("times per month", ntimes.tolist(), "files", len(files))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
