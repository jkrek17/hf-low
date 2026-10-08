"""Match archive HF events (docs/data/hf-lows.json) to pipeline A tracks in all_tracks.csv.gz.

all_tracks has one row per track (peak position/time of the gust index, start/end), not the fixes,
so the match is at track level.  Rule (fixed in PREREGISTRATION.md before any skill was computed):
  candidates = same basin, track [start,end] within 12 h of the archive [start,end];
  distance   = great-circle km from the track's peak position to the archive position at the
               peak time (archive fixes interpolated; clamped to the nearest end outside the span),
               plus 25 km per 6 h of clamping;
  match      = the candidate with the smallest distance, if <= MATCH_KM (500).
Neither the gust value nor the pressure is used to choose.  Each archive event takes one track;
a track may take several archive events only if it is the nearest for each.
"""
import json, os, sys
import numpy as np, pandas as pd

RE = 6371.0
MATCH_KM = 500.0
PEN_KM_PER_6H = 25.0


def hav(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def hours(t):
    t = np.asarray(t, dtype=np.int64)
    d = pd.to_datetime((t // 100).astype(str), format="%Y%m%d")
    return (d.view("int64") // 3_600_000_000_000 + (t % 100)).astype(float) if hasattr(d, "view") else None


def to_h(v):
    v = int(v)
    return pd.Timestamp(year=v // 1000000, month=v // 10000 % 100, day=v // 100 % 100, hour=v % 100).timestamp() / 3600.0


def load_archive(path):
    d = json.load(open(path))
    ev = [dict(zip(d["lowFields"], r)) for r in d["lows"]]
    for e in ev:
        e["fixes"] = [dict(zip(d["fixFields"], f)) for f in e["fixes"]]
    return ev


def archive_pos(e, t):
    """interpolated (lat, lon, clamp hours) at hour t; lon in 0..360"""
    F = [f for f in e["fixes"] if f["lat"] is not None and f["lon"] is not None]
    if not F:
        return None
    th = np.array([to_h(f["date"]) for f in F]); o = np.argsort(th); th = th[o]
    la = np.array([F[i]["lat"] for i in o]); lo = np.array([F[i]["lon"] % 360 for i in o])
    clamp = max(th[0] - t, t - th[-1], 0.0)
    tt = min(max(t, th[0]), th[-1])
    lo_u = np.unwrap(np.radians(lo)); lo_u = np.degrees(lo_u)
    return np.interp(tt, th, la), np.interp(tt, th, lo_u) % 360, clamp


def match_all(T, ev, lo_season=1979):
    """T: all_tracks frame (needs start,end,peak_time,peak_lat,peak_lon,basin).  Returns list of dict."""
    T = T.copy()
    T["s_h"] = [to_h(v) for v in T.start]; T["e_h"] = [to_h(v) for v in T.end]; T["p_h"] = [to_h(v) for v in T.peak_time]
    out = []
    byb = {b: g for b, g in T.groupby("basin")}
    for e in ev:
        if e["season"] < lo_season:
            continue
        g = byb[e["basin"]]
        a0, a1 = to_h(e["start"]), to_h(e["end"])
        c = g[(g.e_h >= a0 - 12) & (g.s_h <= a1 + 12)]
        best, bd = None, 1e9
        if len(c):
            dists = np.full(len(c), 1e9)
            for i, (ph, pla, plo) in enumerate(zip(c.p_h.values, c.peak_lat.values, c.peak_lon.values)):
                r = archive_pos(e, ph)
                if r is None:
                    continue
                la, lo, cl = r
                dists[i] = hav(la, lo, pla, plo) + PEN_KM_PER_6H * cl / 6.0
            j = int(np.argmin(dists))
            if dists[j] <= MATCH_KM:
                best, bd = int(c.track.values[j]), float(dists[j])
        out.append(dict(key=e["basin"] + ":" + e["id"], basin=e["basin"], season=e["season"], cls=e["cls"], peak=e["peak"],
                        arch_minp=e["minP"], track=best, dist_km=bd if best is not None else np.nan))
    return pd.DataFrame(out)
