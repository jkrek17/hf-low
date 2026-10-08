"""Trough and jet predictors at every framework fix. Outcome-free: reads only fix time and position.

ERA5 proxy. Inputs: work/fields/fields_<season>.npz (extract_fields.py) and the framework's fix table
(only track, time, lat, lon are used).

  trough_up      minus the lowest Z500 anomaly (m) among 5.625 deg cells 500-2500 km from the fix at
                 bearings 225-315 deg. Anomaly = Z500 minus the calendar-month mean at the cell over the
                 fit seasons 2004-05..2014-15 (00 and 12 UTC pooled). Larger = deeper trough.
  trough_loc     minus the mean Z500 anomaly of cells within 1000 km (nearest cell if none), m.
  jet_up         maximum 250 hPa wind (kt) among the same upstream cells.
  trough_up_eddy as trough_up with the eddy height (Z500 minus its zonal mean, same latitude and time).
  nsec           number of cells in the upstream sector (0 -> trough_up and jet_up are missing).

usage: features.py FIXES_CSV.GZ FIELDS_DIR OUT_CSV.GZ
"""
import sys, glob
import numpy as np, pandas as pd

LAT = -87.1875 + 5.625 * np.arange(18, 32)
LON = 5.625 * np.arange(64)
FIT_LAST = 2014
KT = 1 / 0.514444
RE = 6371.0


def load(fdir):
    T, Z, U, V = [], [], [], []
    for f in sorted(glob.glob(f"{fdir}/fields_*.npz")):
        d = np.load(f)
        T.append(d["time"]); Z.append(d["z500"]); U.append(d["u250"]); V.append(d["v250"])
    return np.concatenate(T), np.concatenate(Z), np.concatenate(U), np.concatenate(V)


def month_of(hours):
    return pd.to_datetime(np.datetime64("1900-01-01T00") + hours.astype("timedelta64[h]")).month.values


def main(fx, fdir, out):
    t, Z, U, V = load(fdir)
    mon = month_of(t)
    season = np.where(mon >= 6, pd.to_datetime(np.datetime64("1900-01-01T00") + t.astype("timedelta64[h]")).year.values,
                      pd.to_datetime(np.datetime64("1900-01-01T00") + t.astype("timedelta64[h]")).year.values - 1)
    fit = season <= FIT_LAST
    clim = np.stack([Z[fit & (mon == m)].mean(0) for m in range(1, 13)])      # [12, 14, 64]
    A = Z - clim[mon - 1]
    E = Z - Z.mean(2, keepdims=True)
    W = np.hypot(U, V) * KT
    F = pd.read_csv(fx, usecols=["track", "time", "lat", "lon"], dtype={"time": str})
    tt = pd.to_datetime(F.time, format="%Y%m%d%H").values.astype("datetime64[h]")
    hf = ((tt - np.datetime64("1900-01-01T00")) / np.timedelta64(1, "h")).astype(int)
    pos = pd.Series(np.arange(len(t)), index=t)
    ti = pos.reindex(hf).values
    ok_t = ~np.isnan(ti)
    ti = np.where(ok_t, ti, 0).astype(int)
    phi2, lam2 = np.radians(LAT)[:, None], np.radians(LON)[None, :]
    res = {k: np.full(len(F), np.nan) for k in ("trough_up", "trough_loc", "jet_up", "trough_up_eddy")}
    nsec = np.zeros(len(F), int)
    for a in range(0, len(F), 4000):
        b = min(a + 4000, len(F))
        la, lo = np.radians(F.lat.values[a:b]), np.radians(F.lon.values[a:b])
        phi1, lam1 = la[:, None, None], lo[:, None, None]
        dl = lam2[None] - lam1
        d = 2 * RE * np.arcsin(np.sqrt(np.sin((phi2[None] - phi1) / 2) ** 2 +
                                       np.cos(phi1) * np.cos(phi2[None]) * np.sin(dl / 2) ** 2))
        brg = np.degrees(np.arctan2(np.sin(dl) * np.cos(phi2[None]),
                                    np.cos(phi1) * np.sin(phi2[None]) - np.sin(phi1) * np.cos(phi2[None]) * np.cos(dl))) % 360
        sec = (d >= 500) & (d <= 2500) & (brg >= 225) & (brg <= 315)
        loc = d <= 1000
        # nearest cell if none within 1000 km
        near = d == d.reshape(len(d), -1).min(1)[:, None, None]
        loc = np.where(loc.any((1, 2))[:, None, None], loc, near)
        Ab, Eb, Wb = A[ti[a:b]], E[ti[a:b]], W[ti[a:b]]
        sa = np.where(sec, Ab, np.inf).reshape(b - a, -1).min(1)
        se = np.where(sec, Eb, np.inf).reshape(b - a, -1).min(1)
        sw = np.where(sec, Wb, -np.inf).reshape(b - a, -1).max(1)
        has = sec.any((1, 2)) & ok_t[a:b]
        res["trough_up"][a:b] = np.where(has, -sa, np.nan)
        res["trough_up_eddy"][a:b] = np.where(has, -se, np.nan)
        res["jet_up"][a:b] = np.where(has, sw, np.nan)
        res["trough_loc"][a:b] = np.where(ok_t[a:b], -(Ab * loc).sum((1, 2)) / loc.sum((1, 2)), np.nan)
        nsec[a:b] = sec.sum((1, 2))
    O = F[["track", "time"]].copy()
    for k, v in res.items():
        O[k] = np.round(v, 2)
    O["nsec"] = nsec
    O.to_csv(out, index=False)
    print(f"{len(O)} fixes; missing time {np.mean(~ok_t):.4f}; no sector cells {np.mean(nsec == 0):.4f}")
    print(O[["trough_up", "trough_loc", "jet_up", "trough_up_eddy", "nsec"]].describe().round(2).to_string())


if __name__ == "__main__":
    main(*sys.argv[1:4])
