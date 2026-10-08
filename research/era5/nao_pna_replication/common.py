"""Shared pieces for the RA-8 replication (ERA5 PROXY; pipeline A detector and linker, MSLP-only run of PR 38).
Plan: PREREGISTRATION.md. Estimators reused from tele_combos/core.py."""
import datetime as dt
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "research/era5/hf_history"))
sys.path.insert(0, os.path.join(REPO, "research/era5/tele_combos"))
import core            # noqa: E402  (tele_combos estimators)
import track as A      # noqa: E402  (pipeline A basin rule)

KUR = os.path.join(REPO, "research/era5/enso_kuroshio/data")
fs = core.fs


def indices():
    out = {}
    for k in ("NAO", "PNA"):
        s = fs.read_cpc(os.path.join(core.CPC, f"norm.daily.{k.lower()}.index.b500101.current.ascii"))
        out[k + "_lag"] = s.rolling(7, min_periods=7).mean().shift(4)           # days -10..-4
        out[k + "_now"] = s.rolling(7, min_periods=7, center=True).mean()       # days -3..+3
    return pd.DataFrame(out)


def load_tracks():
    """One row per track: first-fix date, and the fix of minimum MSLP (lon, lat, basin there and at genesis)."""
    t = pd.read_csv(os.path.join(KUR, "tracks.csv.gz"), dtype={"t0": str})
    f = pd.read_csv(os.path.join(KUR, "fixes.csv.gz"), dtype={"time": str})
    i = f.groupby("tid").msl.idxmin()        # first occurrence of the minimum
    m = f.loc[i, ["tid", "lat", "lon", "msl"]].rename(columns={"lat": "mlat", "lon": "mlon", "msl": "mmsl"})
    t = t.merge(m, on="tid")
    t["gen"] = pd.to_datetime(t.t0.str[:8], format="%Y%m%d")
    for tag, la, lo in (("min", t.mlat, t.mlon), ("gen", t.lat0, t.lon0)):
        t["basin_" + tag] = [A.basin(a, b) for a, b in zip(la, lo)]
    return t


def table(T, I, basin, A_, B_, s0, s1, now=False, basin_by="min"):
    """Return (D, E) for winters s0..s1: analysis days (3 Dec - 31 Mar) and events with standardised predictors."""
    suf = "_now" if now else "_lag"
    rows = []
    for s in range(s0, s1 + 1):
        d0 = dt.date(s, 12, 1)
        n = (dt.date(s + 1, 4, 1) - d0).days
        for k in range(2, n):                        # drop 1 and 2 Dec
            rows.append((s, k, pd.Timestamp(d0 + dt.timedelta(days=k))))
    D = pd.DataFrame(rows, columns=["season", "day", "date"])
    raw = I.reindex(D.date)[[A_ + suf, B_ + suf]].reset_index(drop=True)
    ok = raw.notna().all(axis=1).values
    D, raw = D[ok].reset_index(drop=True), raw[ok].reset_index(drop=True)
    z = (raw - raw.mean()) / raw.std()
    D["A"], D["B"] = z.iloc[:, 0].values, z.iloc[:, 1].values
    D["AB"] = D.A * D.B
    X = T[T["basin_" + basin_by] == basin].copy()
    X["season"] = X.winter
    X["day"] = (X.gen - pd.to_datetime(X.winter.astype(str) + "-12-01")).dt.days
    X["lon"] = np.where(X.mlon < 180, X.mlon + 360, X.mlon) if basin == "atl" else X.mlon
    X["lat"] = X.mlat
    E = X.merge(D[["season", "day", "A", "B", "AB"]], on=["season", "day"], how="inner")
    E["month"] = E.gen.dt.month
    E["era"] = 0.0
    return D, E
