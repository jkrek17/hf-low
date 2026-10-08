"""Forecast fixes and their outcomes from the ERA5 track points.

One row per in-domain fix at 00 or 12 UTC, with the storm-state predictors
known at that time and the outcomes 24 and 48 h later on the same track.

  msl        central pressure (hPa)
  dp12       msl change over the previous 12 h (hPa), nan if the track is younger
  age        hours since the track's first fix
  speed      translation speed over the previous 6 h (kt), nan at the first fix
  heading    direction of motion over the previous 6 h (deg, toward), for Hart's B
  g800       current gust index (kt): max ocean gust within 800 km
  dp24       msl(t+24) - msl(t) (hPa), nan if the track ends first
  ndr24      normalised deepening rate in Bergerons, -dp24 * sin60 / sin(mean lat) / 24
  hf24,hf48  the gust index reaches 71.7 kt at any fix in (t, t+24] / (t, t+48]
  hf_now     the gust index is already >= 71.7 kt at t

A track that ends inside the window counts as not reaching HF afterwards; the
low has filled or merged. ERA5 tracks are linked 6-hourly (track.py), so a gap
inside a track cannot occur.

usage: fixes.py TRACK_POINTS_CSV OUT_CSV
"""
import sys
import numpy as np, pandas as pd

HF = 71.7
SIN60 = np.sin(np.radians(60))


def build(P):
    P = P.sort_values(["track", "time"]).reset_index(drop=True)
    P["dt"] = pd.to_datetime(P.time, format="%Y%m%d%H")
    out = []
    for tid, g in P.groupby("track", sort=False):
        g = g.reset_index(drop=True)
        n = len(g)
        h = ((g.dt - g.dt.iloc[0]).dt.total_seconds() / 3600).values.astype(int)
        lat, lon, msl = g.lat.values, g.lon.values, g.msl.values
        gi = np.where(g.basin.fillna("").values != "", g.g800.values, 0.0)
        # motion over the previous 6 h only: a centred difference would use the
        # next fix, which is not known at forecast time
        i0 = np.r_[0, np.arange(n - 1)]
        dlat = lat - lat[i0]
        dlon = ((lon - lon[i0] + 180) % 360 - 180) * np.cos(np.radians(lat))
        dh = np.maximum(h - h[i0], 1)
        heading = np.degrees(np.arctan2(dlon, dlat)) % 360
        speed = np.hypot(dlat, dlon) * 60 / dh            # deg -> nmi per h = kt
        heading[0], speed[0] = np.nan, np.nan
        for i in range(n):
            if not isinstance(g.basin[i], str) or g.basin[i] == "" or g.time[i][-2:] not in ("00", "12"):
                continue
            j12 = i - 2 if i >= 2 else None
            j24 = i + 4 if i + 4 < n else None
            ml = np.radians(0.5 * (lat[i] + lat[j24])) if j24 is not None else np.nan
            dp24 = msl[j24] - msl[i] if j24 is not None else np.nan
            out.append(dict(
                track=tid, time=g.time[i], basin=g.basin[i], season=g.season[i],
                lat=lat[i], lon=lon[i], msl=msl[i],
                dp12=msl[i] - msl[j12] if j12 is not None else np.nan,
                age=h[i], speed=speed[i], heading=heading[i], g800=gi[i],
                dp24=dp24, ndr24=-dp24 / 24 * SIN60 / np.sin(ml) if j24 is not None else np.nan,
                hf_now=gi[i] >= HF,
                hf24=bool((gi[i + 1:i + 5] >= HF).any()),
                hf48=bool((gi[i + 1:i + 9] >= HF).any()),
                left24=n - 1 - i < 4))
    return pd.DataFrame(out)


if __name__ == "__main__":
    P = pd.read_csv(sys.argv[1], dtype={"time": str})
    if "season" not in P:
        y, m = P.time.str[:4].astype(int), P.time.str[4:6].astype(int)
        P["season"] = np.where(m >= 6, y, y - 1)
    F = build(P)
    F.to_csv(sys.argv[2], index=False)
    print(len(F), "fixes", F.track.nunique(), "tracks")
