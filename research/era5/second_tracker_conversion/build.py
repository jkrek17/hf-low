"""RA-30 step 2: per-track tables for trackers M and V (ERA5 proxy, 1.5 degrees): basin, first-fix date, peak position,
deepening flags (dp12 <= -3.6 and -2.4 hPa per 12 h at 00/12 UTC fixes), and pipeline A HF labels carried through the PR 116
match at distance D, D/2 and 2D. Plan: PREREGISTRATION.md.
usage: python3 build.py TRACKS_DIR V_MSLP500 OUT_DIR       (TRACKS_DIR holds M_fixes.csv.gz and V_fixes.csv.gz, not committed)
"""
import os, sys
import numpy as np, pandas as pd

TD, VP, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
ERA = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
HFH = f"{ERA}/hf_history/results"
sys.path.insert(0, f"{ERA}/second_tracker")
import common as C
os.makedirs(OUT, exist_ok=True)
DKM = {"M": 500.0, "V": 700.0}
EXPECT = {"M": 83038, "V": 50040}


def basin_of(lat, lon):
    return np.where((lat >= 30) & (lat <= 67) & ((lon >= 262) | (lon <= 10)), "atl",
                    np.where((lat >= 27) & (lat <= 67) & (lon >= 135) & (lon <= 240), "pac", ""))


# ---- pipeline A HF events (same set as PR 116: role event, seasons 2004..2021)
cat = pd.read_csv(f"{HFH}/era5_hf_catalog.csv")
ev = cat[(cat.role == "event") & cat.season.between(2004, 2021)]
fx = pd.read_csv(f"{HFH}/era5_hf_catalog_tracks.csv")
HF = fx[fx.track.isin(ev.track) & (fx.g800 >= 71.7)]
HFX = {k: g[["time", "lat", "lon"]].values for k, g in HF.groupby("track")}


def match_labels(BYT, D):
    """set of tracker tracks matched to at least one A HF event (rule of PR 116: half of the HF fixes within D km, same time)."""
    out = {}
    for k, f in HFX.items():
        hits = {}
        for tm, la, lo in f:
            a = BYT.get(int(tm))
            if a is None:
                continue
            d = C.dist(la, lo, a[:, 1], a[:, 2])
            for t in a[d <= D, 0]:
                hits[int(t)] = hits.get(int(t), 0) + 1
        if hits:
            t = max(hits, key=hits.get)
            if hits[t] >= max(1, 0.5 * len(f)):
                out[k] = t
    return out


rep = []
for tr in ("M", "V"):
    d = pd.read_csv(f"{TD}/{tr}_fixes.csv.gz")
    d["time"] = d.time.astype(np.int64)
    n = d.track.nunique()
    assert n == EXPECT[tr], (tr, n)
    print(tr, "tracks", n, "fixes", len(d), "== expected", EXPECT[tr], flush=True)
    d["hh"] = d.time % 100
    if tr == "M":
        d["p"] = d.val
    else:
        v = pd.read_csv(VP)
        d = d.merge(v, on=["track", "time"], how="left")
        d["p"] = d.mslp500
    # intensity-peak fix: lowest MSLP (M) or largest vorticity (V)
    pk = d.loc[(d.val if tr == "M" else -d.val).groupby(d.track).idxmin()]
    g = d.groupby("track")
    T = pd.DataFrame(dict(first=g.time.min(), nfix=g.size()))
    T = T.join(pk.set_index("track")[["lat", "lon", "time"]].rename(columns=dict(lat="peak_lat", lon="peak_lon", time="peak_time")))
    T["basin"] = basin_of(T.peak_lat.values, T.peak_lon.values)
    T["date"] = pd.to_datetime(T["first"] // 100, format="%Y%m%d")
    # 12 h pressure change at 00/12 UTC fixes whose predecessor 12 h earlier is a fix of the same track
    s = d[d.hh.isin([0, 12])][["track", "time", "p"]].copy()
    s["t"] = pd.to_datetime(s.time.astype(str), format="%Y%m%d%H")
    prev = s[["track", "t", "p"]].copy()
    prev["t"] = prev.t + pd.Timedelta(hours=12)
    m = s.merge(prev, on=["track", "t"], suffixes=("", "_prev"))
    m["dp12"] = m.p - m.p_prev
    mn = m.groupby("track").dp12.min()
    T["min_dp12"] = mn
    T["deep"] = (T.min_dp12 <= -3.6).astype(int)
    T["deep24"] = (T.min_dp12 <= -2.4).astype(int)
    T["deep60"] = (T.min_dp12 <= -6.0).astype(int)
    BYT = {k: g2[["track", "lat", "lon"]].values for k, g2 in d.groupby("time")}
    for lab, f in (("hf", 1.0), ("hf_half", 0.5), ("hf_dbl", 2.0)):
        ml = match_labels(BYT, DKM[tr] * f)
        T[lab] = T.index.isin(set(ml.values())).astype(int)
        rep.append((tr, lab, len(ml), int(T[lab].sum())))
        if lab == "hf":                                     # reproduce PR 116's matched tracks
            e = pd.read_csv(f"{ERA}/second_tracker/results/events_{tr}.csv")
            e = e[e.matched]
            ok = set(e.tr_track) == set(ml.values())
            print(tr, "events matched", len(ml), "vs PR 116", len(e), "same track set:", ok, flush=True)
            rep.append((tr, "repro_PR116_same_tracks", int(ok), len(e)))
    T.reset_index().to_csv(f"{OUT}/tracks_{tr}.csv.gz", index=False)
    print(tr, T.groupby("basin")[["deep", "deep24", "hf"]].agg(["sum", "size"]), flush=True)
pd.DataFrame(rep, columns=["tracker", "item", "a", "b"]).to_csv(f"{OUT}/build_report.csv", index=False)
