"""Raw GFS gust against OPC HF designations out to 96 h (12-hourly leads). Extends research/era5/gfs_raw_gust (PR 182).
Deviations from PR 182, set after the card of 2026-10-14: leads every 12 h to 96 h (PR 182: every 6 h to 48 h, scored at 12 h steps);
linker MAXD 1800 km per 12 h step (900 per 6 h) and 2400 km gate (1200); MIN_LEN 3 fixes (24 h span) instead of 4 at 6 h.
usage: score96.py LOWS_DIR ARCHIVE_DIR OUT_PREFIX     (pipeline A detector on GFS = PROXY fixes; truth = OPC archive Category HF)"""
import os, sys, glob, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/gfs_test"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/hf_history")
import a_tracks as A
import track as T

lows_dir, arch, outp = sys.argv[1:4]
MIN_LEN = 3
rows = []
for f in sorted(glob.glob(os.path.join(lows_dir, "*.csv"))):
    d = pd.read_csv(f, dtype={"cycle": str})
    if not len(d): continue
    cyc = os.path.basename(f)[:10]
    t = A.link_cycle(d); t["cycle"] = cyc
    t["valid_time"] = [A.stamp(cyc, l) for l in t["lead"]]
    t["basin"] = [T.basin(a, b) for a, b in zip(t.lat, t.lon)]
    t["season"] = A.season_of(cyc)
    n = t.groupby("track").size()
    rows.append(t[t.track.isin(n[n >= MIN_LEN].index) & (t.basin != "")])
F = pd.concat(rows, ignore_index=True)
F["G"] = F.g800
print("fixes", len(F), "tracks", F.groupby(["cycle", "track"]).ngroups)

def haversine(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

parts = []
for tag in ("Atl", "Pac"):
    d = pd.read_csv(os.path.join(arch, f"HF_Data_-_{tag}.csv"), dtype={"date": str})
    ok = d.date.str.fullmatch(r"\d{10}", na=False) & pd.to_datetime(d.date, format="%Y%m%d%H", errors="coerce").notna()
    d = d[ok].copy(); d["t"] = pd.to_datetime(d.date, format="%Y%m%d%H"); parts.append(d)
Aa = pd.concat(parts, ignore_index=True); Aa["lon"] = Aa.Longitude % 360; Aa["lat"] = Aa.Latitude
H = Aa[Aa.Category == "HF"].copy(); H["th"] = H.t.values.astype("datetime64[h]").astype(np.int64); H = H.sort_values("th")
th, la, lo = H.th.values, H.lat.values, H.lon.values
v = pd.to_datetime(F.valid_time.astype(str), format="%Y%m%d%H").values.astype("datetime64[h]").astype(np.int64)
a = np.searchsorted(th, v - 3, side="left"); b = np.searchsorted(th, v + 3, side="right")
y = np.zeros(len(F), int); flat, flon = F.lat.values, F.lon.values
for i in np.where(b > a)[0]:
    y[i] = int((haversine(flat[i], flon[i], la[a[i]:b[i]], lo[a[i]:b[i]]) <= 600).any())
F["Y0"] = y
F.to_csv(outp + "_fixes.csv.gz", index=False)

def match_T(g, yy):
    k = int(yy.sum()); gs = np.sort(g)[::-1]
    return float(gs[k]) if 0 <= k < len(gs) else float("-inf")

SE = sorted(F.season.unique()); res = {}
Tm = {(bb, s): match_T(*(lambda d: (d.G.values, d.Y0.values))(F[(F.basin == bb) & (F.lead == 0) & (F.season != s)])) for bb in ("atl", "pac") for s in SE}
F["T"] = [Tm[(bb, s)] for bb, s in zip(F.basin, F.season)]

def score(sel, hit):
    H_ = int((hit & (sel.Y0 == 1)).sum()); FA = int((hit & (sel.Y0 == 0)).sum()); M = int((~hit & (sel.Y0 == 1)).sum()); CN = len(sel) - H_ - FA - M
    N = len(sel); e = ((H_ + FA) * (H_ + M) + (M + CN) * (FA + CN)) / N
    return dict(n=N, hits=H_, fa=FA, miss=M, cn=CN, pod=H_ / max(H_ + M, 1), agree=H_ / max(H_ + FA, 1), hss=((H_ + CN) - e) / (N - e))

for bb in ("atl", "pac"):
    for l in sorted(F.lead.unique()):
        s = F[(F.basin == bb) & (F.lead == l)]
        res[f"{bb}_{l}"] = dict(matched=score(s, (s.G > s["T"]).values), kt64=score(s, (s.G > 64).values), kt72=score(s, (s.G > 72).values))
bands = [(60, 64), (64, 68), (68, 72), (72, 76), (76, 80), (80, 85), (85, 90), (90, 999)]
bt = {}
for l in sorted(F.lead.unique()):
    s = F[F.lead == l]
    bt[int(l)] = {f"{lo_}-{hi_}": (int(((s.G >= lo_) & (s.G < hi_)).sum()), int(s[(s.G >= lo_) & (s.G < hi_)].Y0.sum())) for lo_, hi_ in bands}
json.dump(dict(T={f"{k[0]}_{k[1]}": v_ for k, v_ in Tm.items()}, scores=res, bands=bt, n_fixes=len(F)), open(outp + ".json", "w"), indent=1, default=int)
for k, v_ in res.items(): print(k, {m: (round(x["pod"], 3), round(x["agree"], 3), round(x["hss"], 3)) for m, x in v_.items()})
