"""Shared data preparation for the Greenland-high share study. ERA5 PROXY, pipeline A.
Plan: PREREGISTRATION.md (committed first, 7a74004). Same sample, indices and Poisson machinery as hf-low PR 63."""
import os, sys, warnings
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "freq_split"))
import split as S  # noqa: E402

warnings.filterwarnings("ignore")
S0, S1, WIN = 2004, 2025, "octapr"
SEED = 20261008
BOX = (55.0, 67.0, -50.0, -15.0)  # PR 58's deep-low region: lat lo/hi, lon lo/hi (signed)


def sl(lon):
    return np.where(np.asarray(lon) > 180, np.asarray(lon) - 360, np.asarray(lon))


def prep(tracks, cpc, repo, ghfile, fixes_file, env_file):
    T = S.load_tracks(tracks)
    T = T[T.basin == "atl"].reset_index(drop=True)
    I = S.indices(cpc, repo)
    gh = pd.read_csv(ghfile, parse_dates=["date"]).set_index("date").GH.asfreq("D")
    I["GH_lag"] = gh.rolling(7, min_periods=7).mean().shift(4)
    D = S.window_days(WIN, S0, S1)
    raw = I.reindex(D.date)[["NAO_lag", "GH_lag"]].reset_index(drop=True)
    ok = raw.notna().all(axis=1).values
    D, raw = D[ok].reset_index(drop=True), raw[ok].reset_index(drop=True)
    mu, sd = raw.mean(), raw.std(ddof=0)
    z = (raw - mu) / sd
    Tw = S.assign(T, WIN, S0, S1)
    Tw = Tw.merge(D[["season", "day"]].assign(_in=1), on=["season", "day"], how="inner")
    zi = ((I.reindex(Tw.gen)[["NAO_lag", "GH_lag"]].reset_index(drop=True)) - mu) / sd
    Tw["zNAO"], Tw["zGH"] = zi.NAO_lag.values, zi.GH_lag.values
    Tw["mon"] = Tw.gen.dt.month
    # ---- fix-level to storm-level
    F = pd.read_csv(fixes_file)
    F = F[F.basin == "atl"].sort_values(["track", "time"])
    E = pd.read_csv(env_file)
    F = F.merge(E, on=["track", "time"], how="left")
    first = F.groupby("track").head(1).set_index("track")
    imin = F.loc[F.groupby("track").msl.idxmin()].set_index("track")
    st = pd.DataFrame({
        "lat0": first.lat, "lon0": sl(first.lon), "jet0": first.jet250, "eady0": first.eady,
        "sstgrad0": first.sstgrad, "tcwv0": first.tcwv,
        "latm": imin.lat, "lonm": sl(imin.lon), "ndr_max": F.groupby("track").ndr24.max()})
    Tw = Tw.merge(st, left_on="track", right_index=True, how="left")
    Tw["lonp"] = sl(Tw.peak_lon)
    Tw["box"] = ((Tw.latm >= BOX[0]) & (Tw.latm <= BOX[1]) & (Tw.lonm >= BOX[2]) & (Tw.lonm <= BOX[3])).astype(float)
    Tw.loc[Tw.latm.isna(), "box"] = np.nan
    return dict(T=T, Tw=Tw.reset_index(drop=True), D=D, z=z, raw=raw, mu=mu, sd=sd, F=F)
