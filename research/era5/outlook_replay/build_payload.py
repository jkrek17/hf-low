"""Replay payload: for each fixed case (replay_cases.csv) the gust band (0..6) of every 1-degree ocean point for forecast days 1-4,
the OPC HF/S fixes of the 96 h, and the leave-one-season-out band tables (built without that case's season).
GFS gust is a model forecast (PROXY for sustained wind); OPC's call is an analyst judgement. usage: build_payload.py DMAX_DIR"""
import os, sys, json
import numpy as np, pandas as pd
sys.argv = [sys.argv[0], sys.argv[1], "/home/user/hf-low/data/hf_lows", "/dev/null"]
import tabulate as T                                                       # reuses the grid, mask and max-within-radius
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results"); dm = sys.argv[1]
cases = pd.read_csv(os.path.join(RES, "replay_cases.csv"), dtype={"cycle": str})
pdays = pd.read_csv(os.path.join(RES, "point_days.csv"), dtype={"cycle": str})
pdays["season"] = [int(c[:4]) if int(c[4:6]) >= 10 else int(c[:4]) - 1 for c in pdays.cycle]
ok = T.OK; bas = T.BAS
flat = np.nonzero(ok.ravel())[0]
lat1 = (T.PL.ravel()[flat]).round(2).tolist(); lon1 = (T.PN.ravel()[flat]).round(2).tolist()
def tables(excl):
    out = {}
    for b, name in ((1, "atl"), (2, "pac")):
        for day in (1, 2, 3, 4):
            x = pdays[(pdays.basin == name) & (pdays.day == day) & (pdays.season != excl)].groupby("band")[["n", "hf", "storm"]].sum().reindex(range(7)).fillna(0)
            y = pdays[(pdays.basin == name) & (pdays.day == day) & (pdays.season != excl) & (pdays.hf > 0)].groupby("band").cycle.nunique().reindex(range(7)).fillna(0)
            out[f"{name}{day}"] = {"hf": [round(100 * a / b_, 2) if b_ else 0 for a, b_ in zip(x.hf, x.n)], "st": [round(100 * a / b_, 2) if b_ else 0 for a, b_ in zip(x.storm, x.n)],
                                   "n": [int(v) for v in x.n], "c": [int(v) for v in y]}
    return out
def rle(a):
    s = []; prev = a[0]; n = 1
    for v in a[1:]:
        if v == prev: n += 1
        else: s.append(f"{prev}{n},"); prev = v; n = 1
    s.append(f"{prev}{n},"); return "".join(s)
seasons = sorted(set(cases.season))
payload = {"grid": {"lat": lat1, "lon": lon1, "basin": bas.ravel()[flat].tolist()}, "radius_nm": T.RD, "truth_nm": T.TRUTH_NM,
           "bands_kt": ["<40", "40-48", "48-56", "56-64", "64-72", "72-80", ">=80"], "tables": {str(s): tables(s) for s in seasons}, "tables_all": tables(-1), "cases": []}
for r in cases.itertuples():
    t0 = pd.Timestamp(r.cycle); d = np.load(os.path.join(dm, "cycles", f"{r.cycle}.npz"))["d"]; days = []
    for k in range(4):
        pred = T.maxwithin(d[k].astype(np.float32) / 2, T.RD[k])
        band = np.digitize(pred, T.BANDS) - 1
        days.append(rle(band.ravel()[flat].astype(int).tolist()))
    fx = T.OPC[(T.OPC.t > t0) & (T.OPC.t <= t0 + pd.Timedelta(hours=96))]
    fixes = [[round(float(a), 1), round(float(b), 1), c, int((t - t0).total_seconds() // 3600)] for a, b, c, t in zip(fx.lat, fx.lon, fx.Category, fx.t)]
    payload["cases"].append({"cycle": r.cycle, "season": int(r.season), "why": r.why, "days": days, "fixes": fixes})
json.dump(payload, open(os.path.join(RES, "replay_payload.json"), "w"), separators=(",", ":"))
print(len(flat), "points", os.path.getsize(os.path.join(RES, "replay_payload.json")) / 1e6, "MB")
