#!/usr/bin/env python3
"""Weekly job for the 'pattern loading' note on the live outlook page: score the hemispheric Z500 / U250 pattern from the newest GEFS control run.

    python3 tools/pattern_loading.py [--docs docs] [--flat flat] [--today YYYY-MM-DD] [--force]

What it does (same method as research/era5/basin_bomb_score/ PR 277 and gfs_test/b_index.py; ERA5 proxy, pipeline A, experimental):
  * Cycle dates are fixed on the calendar: 00 UTC on 23 Sep + 7k days, k = 0..29, because the model climatology was built on exactly those dates.
    The job takes the newest such date that is not after today (UTC) and does nothing when the data file already holds it or outside k = 0..29.
  * Pulls the GEFS control (gec00, 0.5 degree) HGT 500 mb and UGRD 250 mb at forecast hours 12, 36, ..., 156 (12 UTC of days 1..7) from the NOAA open-data bucket
    (about 14 byte ranges of ~0.25 MB), decodes them with tools/grib2_t53.py, regrids to the 12 x 64 hemispheric grid (area-weighted, as hemispheric/fields.py),
    averages the seven steps, subtracts the model climatology of week k, projects on the frozen PR 41 EOFs and applies the frozen betas: index = the z500 and u250 PC part of model P.
  * The week it speaks about starts 8 days after the run (1 Oct + 7k days) and lasts seven days. Class: low (bottom fifth of 2004-05 to 2018-19 weeks), typical, high (top fifth).
  * Writes docs/pattern-loading-data.js and flat/pattern-loading-data.js (identical): window.HF_PATTERN = {run, week_start, made, k, basins: {atl|pac: {idx, cls, ref}}}.
    ref carries the class table (HF lows and bombs per week, share of weeks with at least one) from data/pattern_loading/params.json, so the page needs no other file.
Frozen inputs: data/pattern_loading/frozen.npz and params.json (written by tools/pattern_loading_freeze.py). Standard library plus numpy. Exits 1 if a pull fails.
GEFS v12 reforecast climatology is applied to the operational control run."""
import argparse, datetime as dt, json, os, struct, sys, time, urllib.error, urllib.request
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import grib2_t53 as G

BUCKET = "https://noaa-gefs-pds.s3.amazonaws.com"; UA = "hf-low-pattern-loading/0.1 (github.com/jkrek17/hf-low)"
HOURS = [12 + 24 * i for i in range(7)]; NWEEK = 30; K = 10
DATA_NAME = "pattern-loading-data.js"
FZ = np.load(os.path.join(ROOT, "data", "pattern_loading", "frozen.npz")); PR = json.load(open(os.path.join(ROOT, "data", "pattern_loading", "params.json")))

# --- regrid 0.25 degree (721 x 1440) -> 5.625 degree, 12 northern rows x 64 columns; copy of hemispheric/fields.py _lat_weights / _lon_weights ---
LAT32 = -87.1875 + 5.625 * np.arange(32)
def _lat_weights():
    s = np.arange(721); a = np.clip(90 - 0.25 * s - 0.125, -90, 90); b = np.clip(90 - 0.25 * s + 0.125, -90, 90)
    e0 = np.maximum(LAT32 - 2.8125, -90); e1 = np.minimum(LAT32 + 2.8125, 90); W = np.zeros((32, 721))
    for k in range(32):
        lo = np.maximum(a, e0[k]); hi = np.minimum(b, e1[k]); W[k] = np.maximum(0, np.sin(np.radians(hi)) - np.sin(np.radians(lo)))
    return W / W.sum(1, keepdims=True)
def _lon_weights():
    s = np.arange(1440); a = 0.25 * s - 0.125; b = 0.25 * s + 0.125; W = np.zeros((64, 1440))
    for j in range(64):
        for sh in (-360, 0, 360):
            W[j] += np.maximum(0, np.minimum(b, 5.625 * j + 2.8125 + sh) - np.maximum(a, 5.625 * j - 2.8125 + sh))
    return W / W.sum(1, keepdims=True)
WLAT = _lat_weights()[20:32]; WLON = _lon_weights()
def regrid(a): return WLAT @ a @ WLON.T
def to_quarter(a):
    """0.5 degree (361 x 720) -> 0.25 degree (721 x 1440), linear in lat and (periodic) lon; original points kept exactly (as gfs_test/b_index.py)."""
    lat = np.empty((721, a.shape[1])); lat[0::2] = a; lat[1::2] = 0.5 * (a[:-1] + a[1:])
    out = np.empty((721, 1440)); out[:, 0::2] = lat
    out[:, 1:-1:2] = 0.5 * (lat[:, :-1] + lat[:, 1:]); out[:, -1] = 0.5 * (lat[:, -1] + lat[:, 0])
    return out

# --- fetching ---
class Missing(Exception):
    pass
def get(u, rng=None, tries=6):
    h = {"User-Agent": UA}; want = None
    if rng:
        h["Range"] = rng; a, b = rng.split("=")[1].split("-"); want = int(b) - int(a) + 1
    err = None; n404 = 0
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=120) as r:
                b = r.read(); cl = r.headers.get("Content-Length")
            if (cl is not None and len(b) != int(cl)) or (want is not None and len(b) != want): raise IOError("short read")
            return b
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                n404 += 1
                if n404 >= 2: raise Missing(u)
            err = e
        except Exception as e:
            err = e
        time.sleep(min(2 ** k, 20))
    raise RuntimeError(f"fetch failed {u}: {err}")

def record(base, var, lev, hour):
    """One GRIB2 message (HGT 500 mb / UGRD 250 mb at `hour` h) of the control file, by byte range from its .idx (control rows only)."""
    L = [l.split(":") for l in get(base + ".idx").decode().strip().splitlines()]
    for i, p in enumerate(L):
        if p[3] == var and p[4] == lev and p[5] == f"{hour} hour fcst":
            if i + 1 >= len(L): raise RuntimeError("wanted record is last in idx")
            return get(base, f"bytes={int(p[1])}-{int(L[i + 1][1]) - 1}")
    raise Missing(base + f" no {var} {lev}")

def cycle_for(today):
    """(run date, k) of the newest fixed cycle date 23 Sep + 7k not after today, or None outside k = 0..29."""
    y = today.year if today >= dt.date(today.year, 9, 23) else today.year - 1
    n = (today - dt.date(y, 9, 23)).days // 7
    return (dt.date(y, 9, 23) + dt.timedelta(days=7 * n), n) if 0 <= n < NWEEK else None

def window_maps(fetch_msg, decode):
    out = {}
    for v in ("z500", "u250"):
        acc = [regrid(to_quarter(decode(fetch_msg(v, h)))) for h in HOURS]
        out[v] = np.mean(acc, axis=0)
    return out

def project(maps, k):
    res = {}
    for b, P in PR["basins"].items():
        tot = 0.0
        for v, beta in (("z500", P["beta_z"]), ("u250", P["beta_u"])):
            x = ((maps[v] - FZ["clim_" + v][k]) * FZ["w"]).reshape(-1)
            sc = ((x - FZ["mean_" + v]) @ FZ["vec_" + v].T) / FZ["sd_" + v]
            tot += float(sc @ np.asarray(beta))
        cls = "low" if tot < P["low_below"] else ("high" if tot >= P["high_from"] else "typical")
        res[b] = dict(idx=round(tot, 4), cls=cls, ref=P["table"])
    return res

def live_fetch(D):
    base = f"{BUCKET}/gefs.{D:%Y%m%d}/00/atmos/pgrb2ap5/gec00.t00z.pgrb2a.0p50.f"
    spec = {"z500": ("HGT", "500 mb"), "u250": ("UGRD", "250 mb")}
    return lambda v, h: record(base + f"{h:03d}", spec[v][0], spec[v][1], h)
def live_decode(msg): return G.decode_message(msg, grid=(720, 361)).astype(np.float64)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docs", default=os.path.join(ROOT, "docs")); ap.add_argument("--flat", default=os.path.join(ROOT, "flat"))
    ap.add_argument("--today"); ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    today = dt.date.fromisoformat(a.today) if a.today else dt.datetime.now(dt.timezone.utc).date()
    c = cycle_for(today)
    if c is None: print("Outside the 30 weekly cycles (23 Sep to 14 Apr); nothing to do."); return 0
    D, k = c; run = f"{D:%Y%m%d}00"; path = os.path.join(a.docs, DATA_NAME)
    if os.path.exists(path) and f'"run":"{run}"' in open(path, encoding="utf-8").read() and not a.force:
        print(f"Run {run} already in {DATA_NAME}; nothing to do."); return 0
    try:
        base = f"{BUCKET}/gefs.{D:%Y%m%d}/00/atmos/pgrb2ap5/gec00.t00z.pgrb2a.0p50.f{HOURS[-1]:03d}.idx"; get(base, tries=3)
    except Missing:
        print(f"GEFS run {run} not posted yet; try again later."); return 0
    res = project(window_maps(live_fetch(D), live_decode), k)
    S = dt.date(D.year if D.month >= 9 else D.year - 1, 10, 1) + dt.timedelta(days=7 * k)
    doc = {"run": run, "k": k, "week_start": S.isoformat(), "made": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%MZ"), "basins": res}
    txt = "window.HF_PATTERN=" + json.dumps(doc, separators=(",", ":")) + ";\n"
    for d in (a.docs, a.flat):
        os.makedirs(d, exist_ok=True); open(os.path.join(d, DATA_NAME), "w", encoding="utf-8").write(txt)
    print(run, {b: (r["idx"], r["cls"]) for b, r in res.items()})
    return 0

if __name__ == "__main__":
    sys.exit(main())
