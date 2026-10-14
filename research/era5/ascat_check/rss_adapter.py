"""Source adapter: RSS ASCAT v02.1 daily 0.25 degree bytemaps (data.remss.com/ascat) -> cell table for collocate.py.

Each daily file (one per satellite per UTC day) holds, per ascending/descending segment, 0.25 degree cells with the
time of day (6 min steps), 10 m wind speed (0.2 m/s steps; codes > 250 are land/ice/no data), wind direction, the rain
byte (bit 0 = scatterometer rain flag, bit 1 = radiometer rain data present, bits 2-7 = radiometer rain code) and the
sum-of-squares residual. Format: RSS readme_ascat_v02.1.txt, ascat_bmap_support/python/ascat_daily.py.
Satellites: Metop-A (metopa, v02.1, 2007-10 to 2021-11; the reprocessed v02.1.1 folder was empty before 2016 so v02.1 is used for all three), Metop-B (metopb, v02.1, 2012-), Metop-C (metopc, v02.1, 2019-).

Subcommands
  inventory OUT.csv            crawl month directories for the daily file names (small)
  extract  INV.csv OUTDIR      download the files needed by the 00/12 UTC pipeline A sample, keep only cells within 800 km of
                               a sample fix and within +-1.5 h of it, write OUTDIR/cells_<sensor>_<year>.npz (resumable)
Nothing is kept of a daily file except those cells. The compact cell table is the input to collocate.py (loaded with load_cells)."""
import sys, os, re, gzip, io, time, urllib.request, concurrent.futures as cf
import numpy as np, pandas as pd

BASE = "https://data.remss.com/ascat"
SENSORS = {"A": ("metopa", "bmaps_v02.1", range(2007, 2022)),
           "B": ("metopb", "bmaps_v02.1", range(2012, 2027)),
           "C": ("metopc", "bmaps_v02.1", range(2019, 2027))}
WIN_H = 1.5
RADIUS = 800.0
KT = 1.94384449
RE = 6371.0


def get(url, tries=20, timeout=120):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            err = e
        except Exception as e:
            err = e
        time.sleep(min(2 ** k, 30))
    raise RuntimeError(f"{url}: {err}")


def crawl(sensor):
    d, ver, years = SENSORS[sensor]
    def one(ym):
        y, m = ym
        b = get(f"{BASE}/{d}/{ver}/y{y}/m{m:02d}/")
        if b is None:
            return []
        out = []
        txt = b.decode()
        for mm in re.finditer(r'(\d+)\s+<A HREF="([^"]+?(\d{8})_?v[0-9.]+\.gz)">', txt):
            out.append(dict(sensor=sensor, date=mm.group(3), size=int(mm.group(1)), url="https://data.remss.com" + mm.group(2)))
        return out
    with cf.ThreadPoolExecutor(3) as ex:
        res = list(ex.map(one, [(y, m) for y in years for m in range(1, 13)]))
    return [r for x in res for r in x]


def sample_windows():
    """00/12 UTC pipeline A sample storms with stored index > 0: (time_utc_hours_since_1970, lat, lon)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import collocate as C
    S = C.storm_sample()
    t = (S.time.values.astype("datetime64[s]").astype("int64") / 3600.0)
    return S, t


def gc_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


LATS = -89.875 + 0.25 * np.arange(720)
LONS = 0.125 + 0.25 * np.arange(1440)


def read_day(raw):
    a = np.frombuffer(gzip.decompress(raw), np.uint8)
    if a.size != 2 * 5 * 720 * 1440:
        raise ValueError("unexpected file size")
    return a.reshape(2, 5, 720, 1440)


def cells_for_day(a, date, st, sth, slat, slon):
    """a: (2,5,720,1440) bytes; date: datetime64[D]; st/sth/slat/slon: arrays of the sample fixes (hours since 1970, lat, lon)
    whose +-1.5 h window touches this day. Returns unique rows (lat_idx, lon_idx, minutes_since_1970, wind_byte, rain_byte)."""
    day0 = date.astype("datetime64[h]").astype("int64")
    tm = a[:, 0].astype(np.int16)
    ok = tm <= 250
    seg, ii, jj = np.nonzero(ok)
    tcell = day0 + tm[seg, ii, jj] * 6.0 / 60.0              # hours since 1970 (6 min steps)
    wind = a[seg, 1, ii, jj]
    rain = a[seg, 3, ii, jj]
    out = []
    for h, la, lo in zip(sth, slat, slon):
        sel = np.abs(tcell - h) <= WIN_H + 1e-9
        if not sel.any():
            continue
        la_c, lo_c = LATS[ii[sel]], LONS[jj[sel]]
        near = gc_km(la_c, lo_c, la, lo) <= RADIUS
        k = np.nonzero(sel)[0][near]
        if k.size:
            out.append(np.column_stack([ii[k], jj[k], np.rint(tcell[k] * 60).astype(np.int64), wind[k], rain[k]]))
    return np.unique(np.vstack(out), axis=0) if out else np.zeros((0, 5), np.int64)


def extract(inv_path, outdir, workers=5):
    os.makedirs(outdir, exist_ok=True)
    inv = pd.read_csv(inv_path, dtype={"date": str})
    S, t = sample_windows()
    S = S.assign(th=t)
    th = np.sort(np.unique(np.round(t).astype(np.int64)))
    # storms by hour for fast lookup: all fixes with the same cache time share windows but differ in position
    fixes = S[["th", "peak_lat", "peak_lon"]].reset_index(drop=True)
    fixes_h = fixes.th.values
    def job(row):
        date = np.datetime64(f"{row.date[:4]}-{row.date[4:6]}-{row.date[6:]}")
        d0 = date.astype("datetime64[h]").astype("int64")
        m = (fixes_h > d0 - WIN_H - 1e-9) & (fixes_h < d0 + 24 + WIN_H + 1e-9)
        if not m.any():
            return row.sensor, row.date, None
        raw = get(row.url)
        if raw is None:
            return row.sensor, row.date, None
        a = read_day(raw)
        f = fixes[m]
        rows = cells_for_day(a, date, None, f.th.values, f.peak_lat.values, f.peak_lon.values)
        return row.sensor, row.date, rows
    done = set(os.listdir(outdir))
    for (sen, yr), g in inv.assign(yr=inv.date.str[:4]).groupby(["sensor", "yr"]):
        fn = f"cells_{sen}_{yr}.npz"
        if fn in done:
            continue
        t0 = time.time()
        with cf.ThreadPoolExecutor(workers) as ex:
            res = list(ex.map(job, [r for r in g.itertuples()]))
        have = [r for r in res if r[2] is not None and len(r[2])]
        dates = [r[1] for r in res if r[2] is not None]
        arr = np.vstack([r[2] for r in have]) if have else np.zeros((0, 5))
        np.savez_compressed(f"{outdir}/{fn}", rows=arr.astype(np.int64), dates_read=np.array(dates))
        print(sen, yr, len(g), "files", len(dates), "read", len(arr), "cell rows", round(time.time() - t0), "s", flush=True)


def load_cells(outdir, rain_rule="scat"):
    """Cell table for collocate.py: sensor, t (hours since 1900), lat, lon, wind_kt, ok. rain_rule: 'scat' = drop cells with the
    scatterometer rain flag (primary, registered); 'scat+rad' = also drop cells with radiometer rain > 0 (descriptive sensitivity)."""
    parts = []
    off = (np.datetime64("1970-01-01T00", "h") - np.datetime64("1900-01-01T00", "h")) / np.timedelta64(1, "h")
    for fn in sorted(os.listdir(outdir)):
        if not fn.startswith("cells_"):
            continue
        z = np.load(f"{outdir}/{fn}")["rows"]
        if not len(z):
            continue
        sen = fn.split("_")[1]
        rain = z[:, 4].astype(np.int64)
        ok = (rain & 1) == 0
        if rain_rule == "scat+rad":
            ok &= ~(((rain >> 2) & 63) > 1)
        parts.append(pd.DataFrame(dict(sensor=sen, t=z[:, 2] / 60.0 + off, lat=LATS[z[:, 0]], lon=LONS[z[:, 1]],
                                       wind_kt=z[:, 3] * 0.2 * KT, ok=ok.astype(np.int8))))
    return pd.concat(parts, ignore_index=True).sort_values("t").reset_index(drop=True)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "inventory":
        rows = []
        for s in SENSORS:
            r = crawl(s)
            print(s, len(r), "daily files", sum(x["size"] for x in r) / 1e9, "GB")
            rows += r
        pd.DataFrame(rows).to_csv(sys.argv[2], index=False)
    elif cmd == "extract":
        extract(sys.argv[2], sys.argv[3])
