"""Shared loaders for RA-12 (parent-daughter HF pairs). Pipeline A (research/era5/hf_history), ERA5 PROXY.

Tracks: HF event = pipeline A track with gust800_kt >= 71.7, season >= 2004 (decision 1).
Genesis proxy = first fix in intensity/results/fixes_2004.csv.gz (00/12 UTC fixes below 1010 hPa); it is the tracker's
first detection, not true genesis. Observed genesis = lifecycle_events.gen_obs (first fix >= 1000 hPa and north of 21N).
"""
import os
import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
E5 = os.path.join(REPO, "research/era5")
RE = 6371.0
THR = 71.7


def tparse(s):
    return pd.to_datetime(pd.Series(np.asarray(s)).astype(str), format="%Y%m%d%H")


def hav(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def bearing(lat1, lon1, lat2, lon2):
    """Initial great-circle bearing (deg clockwise from north) from point 1 to point 2."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    y = np.sin(dl) * np.cos(p2)
    x = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return (np.degrees(np.arctan2(y, x)) + 360) % 360


def load():
    T = pd.read_csv(os.path.join(E5, "hf_history/results/all_tracks.csv.gz"), dtype={"start": str, "end": str, "peak_time": str})
    T = T[T.basin.isin(["atl", "pac"]) & (T.gust800_kt >= THR) & (T.season >= 2004)].copy()
    L = pd.read_csv(os.path.join(E5, "hf_history/results/lifecycle_events.csv"))[["track", "gen_obs", "tc"]]
    T = T.merge(L, on="track", how="left")
    F = pd.read_csv(os.path.join(E5, "intensity/results/fixes_2004.csv.gz"))
    F = F[F.track.isin(set(T.track))].copy()
    F["t"] = (tparse(F.time) - pd.Timestamp("2000-01-01")).dt.total_seconds().values / 3600
    F = F.sort_values(["track", "t"])
    T["tpeak"] = (tparse(T.peak_time) - pd.Timestamp("2000-01-01")).dt.total_seconds().values / 3600
    first = F.groupby("track").first()[["t", "lat", "lon"]].rename(columns={"t": "t0", "lat": "lat0", "lon": "lon0"})
    last = F.groupby("track").last()[["t"]].rename(columns={"t": "t1"})
    T = T.merge(first, on="track").merge(last, on="track")
    T["tstart"] = (tparse(T.start) - pd.Timestamp("2000-01-01")).dt.total_seconds().values / 3600
    T["lag0"] = T.t0 - T.tstart   # hours between tracker start and first 00/12 fix
    return T.reset_index(drop=True), {k: v for k, v in F.groupby("track")}


def pos_at(fx, t):
    """Linear-interpolated (lat, lon, heading) of a track at time t (inside [t0, t1]); None outside."""
    tt = fx.t.values
    if t < tt[0] or t > tt[-1]:
        return None
    lat, lon = fx.lat.values, np.unwrap(np.radians(fx.lon.values)) * 180 / np.pi
    la, lo = np.interp(t, tt, lat), np.interp(t, tt, lon) % 360
    i = min(max(np.searchsorted(tt, t) , 1), len(tt) - 1) if len(tt) > 1 else 0
    hd = bearing(lat[i - 1], lon[i - 1] % 360, lat[i], lon[i] % 360) if len(tt) > 1 else np.nan
    return la, lo, hd
