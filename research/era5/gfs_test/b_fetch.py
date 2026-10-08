"""Stage B, phase 1: pull GEFS v12 reforecast control (c00) Z500 and U250 for the forecast-day-1..7 12 UTC steps.

Bucket noaa-gefs-retrospective, keys GEFSv12/reforecast/YYYY/YYYYMMDD00/c00/Days:1-10/{hgt,ugrd}_pres_abv700mb_YYYYMMDD00_c00.grib2
with a .idx beside each. Only the HGT 500 mb and UGRD 250 mb records at forecast hours 12, 36, 60, 84, 108, 132, 156
(12 UTC of forecast days 1..7 from a 00 UTC initialisation) are read, by HTTP byte range. Each record is saved as the raw GRIB2
message in WORK/raw/YYYYMMDD/{z500,u250}_fHHH.grb2 (outside the repo). Initialisations: 00 UTC on S-8 for S = 1 Oct + 7k, k = 0..29,
seasons 2004..2018. Resumable (a finished init has WORK/raw/YYYYMMDD/DONE), at most 4 parallel requests, retries with backoff.
Every byte downloaded (idx files and records) is appended to WORK/bytes.log.

usage: b_fetch.py [WORK] [--probe]      (--probe: first init only, print bytes)
"""
import os, sys, time, json, threading, datetime as dt
from concurrent.futures import ThreadPoolExecutor
import requests

BUCKET = "https://noaa-gefs-retrospective.s3.amazonaws.com"
HOURS = [12 + 24 * i for i in range(7)]            # 12 UTC of forecast days 1..7
SEASONS = range(2004, 2019)
NWEEK = 30
VARS = {"z500": ("hgt_pres_abv700mb", "HGT", "500 mb"), "u250": ("ugrd_pres_abv700mb", "UGRD", "250 mb")}
MAXPAR = 4

_lock = threading.Lock()
_sem = threading.Semaphore(MAXPAR)


def inits():
    out = []
    for s in SEASONS:
        for k in range(NWEEK):
            S = dt.date(s, 10, 1) + dt.timedelta(days=7 * k)
            out.append((s, k, S, S - dt.timedelta(days=8)))
    return out


def key(i0, var):
    return f"GEFSv12/reforecast/{i0.year}/{i0:%Y%m%d}00/c00/Days:1-10/{VARS[var][0]}_{i0:%Y%m%d}00_c00.grib2"


class Missing(Exception):
    pass


def get(work, url, rng=None, tag=""):
    h = {"Range": f"bytes={rng[0]}-{rng[1]}"} if rng else {}
    err = None
    for a in range(7):
        try:
            with _sem:
                r = requests.get(url, headers=h, timeout=120)
            if r.status_code in (403, 404):
                # the first full run saw spurious 404s on files that exist (checked by listing); retry before calling a file missing
                if a >= 4:
                    raise Missing(f"{r.status_code} {url}")
                time.sleep(2 ** a)
                continue
            r.raise_for_status()
            n = len(r.content)
            if rng and n != rng[1] - rng[0] + 1:
                raise IOError(f"short read {n} for {rng}")
            with _lock:
                open(os.path.join(work, "bytes.log"), "a").write(f"{dt.datetime.now(dt.timezone.utc).isoformat()}Z {n} {tag} {url.split('amazonaws.com/')[1]}\n")
            return r.content
        except Missing:
            raise
        except Exception as e:
            err = e
            time.sleep(min(60, 2 ** a))
    raise err


def records(idx_text, hgt_or_ugrd, lev):
    """{hour: (start, end)} for the records var/level at the wanted hours; end = next record start - 1."""
    L = idx_text.strip().split("\n")
    starts = [int(l.split(":")[1]) for l in L]
    out = {}
    for i, l in enumerate(L):
        p = l.split(":")
        if p[3] == hgt_or_ugrd and p[4] == lev and p[5].endswith(" hour fcst") and "-" not in p[5]:
            h = int(p[5].split()[0])
            if h in HOURS:
                if i + 1 >= len(L):
                    raise RuntimeError("wanted record is the last in the idx; end unknown")
                out[h] = (starts[i], starts[i + 1] - 1)
    return out


def fetch_init(work, i0):
    d = os.path.join(work, "raw", f"{i0:%Y%m%d}")
    if os.path.exists(os.path.join(d, "DONE")):
        return "cached"
    os.makedirs(d, exist_ok=True)
    for var, (_, vn, lev) in VARS.items():
        url = f"{BUCKET}/{key(i0, var)}"
        idx = get(work, url + ".idx", tag="idx").decode()
        rec = records(idx, vn, lev)
        miss = [h for h in HOURS if h not in rec]
        if miss:
            raise Missing(f"{i0} {var} hours missing from idx: {miss}")
        for h in HOURS:
            fn = os.path.join(d, f"{var}_f{h:03d}.grb2")
            if os.path.exists(fn) and os.path.getsize(fn) == rec[h][1] - rec[h][0] + 1:
                continue
            b = get(work, url, rec[h], tag=f"{var}_f{h:03d}")
            open(fn + ".part", "wb").write(b)
            os.replace(fn + ".part", fn)
    open(os.path.join(d, "DONE"), "w").write("ok\n")
    return "ok"


def main():
    work = next((a for a in sys.argv[1:] if not a.startswith("--")), "/home/claude/gfs_work_b")
    os.makedirs(work, exist_ok=True)
    todo = inits()
    if "--probe" in sys.argv:
        todo = todo[:1]
    status = {}
    def job(t):
        s, k, S, i0 = t
        try:
            return t, fetch_init(work, i0)
        except Missing as e:
            return t, f"missing: {e}"
        except Exception as e:
            return t, f"error: {e!r}"
    t0 = time.time()
    with ThreadPoolExecutor(MAXPAR) as ex:
        for n, (t, st) in enumerate(ex.map(job, todo)):
            status[f"{t[3]:%Y%m%d}"] = dict(season=t[0], week=t[1], S=str(t[2]), status=st)
            if n % 25 == 0 or not st in ("ok", "cached"):
                print(n, t[3], st, f"{time.time()-t0:.0f}s", flush=True)
    json.dump(status, open(os.path.join(work, "fetch_status.json"), "w"), indent=1)
    bad = {k: v for k, v in status.items() if v["status"] not in ("ok", "cached")}
    tot = sum(int(l.split()[1]) for l in open(os.path.join(work, "bytes.log")))
    print(f"inits {len(status)}, not fetched {len(bad)}, total bytes in bytes.log {tot/1e9:.3f} GB")
    for k, v in bad.items():
        print(k, v)


if __name__ == "__main__":
    main()
