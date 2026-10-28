"""ASCAT (RSS v02.1 daily bytemaps) -> 1 degree cell wind maxima and coverage for each 00 UTC valid time, both basins.

Same cell grid and ocean mask as ens_fetch.py. For valid time V: all sensors' ocean 0.25 degree cells with time within +-3 h of V,
scatterometer rain flag (rain byte bit 0) dropped (rule of research/era5/ascat_check/PREREGISTRATION.md). Per 1 degree cell:
max wind (uint8, 2 x kt) and number of 0.25 degree ocean sub-cells with an observation (0..16) and number of ocean sub-cells (0..16).
Resumable, one npz per valid day. Needs the ascat_check adapter for download/parse. Work dir ENS_WORK."""
import os, sys, gzip, time, datetime as dt, concurrent.futures as cf
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "era5", "ascat_check"))
import rss_adapter as R, ens_fetch as F

INV = os.environ.get("RSS_INV", "/mnt/project-files/ascat-check/rss_inventory.csv")
KT = R.KT
HALF = 3.0


def day_cells(raw, date):
    """-> hours-since-epoch, ilat, ilon, wind_kt (float32) of valid, non-rain-flagged ocean cells in one daily file."""
    a = R.read_day(raw)
    day0 = np.datetime64(date).astype("datetime64[h]").astype("int64")
    tm = a[:, 0].astype(np.int16); wind = a[:, 1]; rain = a[:, 3]
    ok = (tm <= 250) & (wind <= 250) & ((rain & 1) == 0)
    seg, ii, jj = np.nonzero(ok)
    return day0 + tm[seg, ii, jj] * 0.1, ii, jj, wind[seg, ii, jj].astype(np.float32) * 0.2 * KT


def grid_global(h, ii, jj, w, V):
    """cells within +-3 h of V on the 0.25 degree (721 x 1440, north-up GEFS grid) -> max wind and obs flag."""
    s = np.abs(h - V) <= HALF
    # RSS lat index 0 = -89.875 (south-up, 720 rows); GEFS row i = lat 90 - 0.25 i. Use centre-of-cell lats: RSS row r -> GEFS row 358-r... map by 0.125 offset
    gi = 719 - ii[s]          # RSS cell centre lat = -89.875 + .25 r  ->  north-up row 719 - r (cell between GEFS rows; consistent offset for every cell)
    mx = np.zeros((720, 1440), np.float32); ob = np.zeros((720, 1440), bool)
    np.maximum.at(mx, (gi, jj[s]), w[s]); ob[gi, jj[s]] = True
    return mx, ob


def to_cells(mx, ob):
    """0.25 degree (720 x 1440; north-up rows = cell centres 89.875 ... -89.875) -> {basin: (wind uint8, covered, ocean)} on 1 degree cells."""
    land = F.land_mask()[:720]      # GEFS point rows 0..719 ~ cell rows (half-cell offset ignored, 0.125 degree)
    out = {}
    for k, b in F.BOX.items():
        iy0 = (90 - b["n"]) * 4; ny = b["n"] - b["s"]
        ix = (np.arange((b["e"] - b["w"]) * 4) + b["w"] * 4) % 1440
        sl = (slice(iy0, iy0 + ny * 4), ix)
        m = mx[sl[0]][:, ix]; o = ob[sl[0]][:, ix]; oc = ~land[sl[0]][:, ix]
        red = lambda x: x.reshape(ny, 4, -1, 4)
        wind = np.clip(np.round(red(np.where(oc, m, 0)).max(axis=(1, 3)) * 2), 0, 255).astype(np.uint8)
        out[k] = (wind, red(o & oc).sum(axis=(1, 3)).astype(np.uint8), red(oc).sum(axis=(1, 3)).astype(np.uint8))
    return out


def day_job(date, files):
    """date = valid day (YYYYMMDD). files: inventory rows for that day and the day before."""
    p = os.path.join(F.WORK, "ascat", date + ".npz")
    if os.path.exists(p): return p
    V = np.datetime64(f"{date[:4]}-{date[4:6]}-{date[6:]}").astype("datetime64[h]").astype("int64")
    mx = np.zeros((720, 1440), np.float32); ob = np.zeros((720, 1440), bool); n = 0
    for r in files:
        raw = R.get(r.url)
        if raw is None: continue
        d = r.date; h, ii, jj, w = day_cells(raw, f"{d[:4]}-{d[4:6]}-{d[6:]}")
        m2, o2 = grid_global(h, ii, jj, w, V); mx = np.maximum(mx, m2); ob |= o2; n += 1
    c = to_cells(mx, ob)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    np.savez_compressed(p + ".tmp.npz", nfiles=n, **{f"{k}_{x}": c[k][i] for k in c for i, x in enumerate(("wind", "cov", "oc"))})
    os.replace(p + ".tmp.npz", p); return p


def valid_days(seasons):
    out = []
    for y in seasons:
        for d in F.run_dates(y):
            r = dt.datetime.strptime(d, "%Y%m%d")
            out += [(r + dt.timedelta(days=k)).strftime("%Y%m%d") for k in (1, 2, 3, 4)]
    return sorted(set(out))


def jobs(seasons):
    inv = pd.read_csv(INV, dtype={"date": str})
    by = {d: g for d, g in inv.groupby("date")}
    out = []
    for v in valid_days(seasons):
        prev = (dt.datetime.strptime(v, "%Y%m%d") - dt.timedelta(days=1)).strftime("%Y%m%d")
        rows = [r for d in (prev, v) if d in by for r in by[d].itertuples()]
        out.append((v, rows))
    return out


if __name__ == "__main__":
    seasons = [int(x) for x in sys.argv[1:]]
    J = jobs(seasons); print(len(J), "valid days;", sum(len(r) for _, r in J), "file reads", flush=True)
    F.land_mask(); t0 = time.time()
    with cf.ThreadPoolExecutor(3) as ex:
        for i, p in enumerate(ex.map(lambda a: day_job(*a), J)):
            if i % 50 == 0: print(i, round(time.time() - t0), "s", flush=True)
    print("done", flush=True)
