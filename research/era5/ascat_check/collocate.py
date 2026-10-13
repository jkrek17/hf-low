"""Collocate scatterometer cells with ERA5 gust around pipeline A storm fixes.

Input scatterometer cell table (any source adapter must write this): csv/csv.gz with columns
  sensor, time (UTC, ISO), lat, lon (-180..360 ok), wind_kt, ok (1 = passes rain/ambiguity/quality flags)
usage: collocate.py cells.csv.gz out_prefix [--sample-rows N]
Writes out_prefix_footprints.csv (one row per storm x sensor pass) and out_prefix_pairs.csv (matched cells)."""
import sys, numpy as np, pandas as pd
import era5side as E

WINDOW_H = 1.5
MIN_COVER = 0.40


def storm_sample():
    S = E.sample()
    S = S[S.gust800_kt > 0].copy()                  # PREREGISTRATION deviation 1
    return S


def collocate(cells, S=None, ocean=None, idx=None, with_pairs=True):
    S = storm_sample() if S is None else S
    ocean = E.lsm_ocean() if ocean is None else ocean
    idx = E.cache_index() if idx is None else idx
    S = S[S.h1900.isin(idx)].copy()
    S["fr"] = [idx[h] for h in S.h1900]
    S = S.sort_values("fr")
    cells = cells[cells.ok == 1].copy()
    cells["t"] = pd.to_datetime(cells.time).astype("datetime64[s]").astype("int64") / 3600.0 + \
        E.hours1900("1970-01-01T00")           # hours since 1900
    cells = cells.sort_values("t")
    tt = cells.t.values
    foot, pairs = [], []
    for _, row in S.iterrows():
        a, b = np.searchsorted(tt, [row.h1900 - WINDOW_H, row.h1900 + WINDOW_H + 1e-9])
        if b <= a:
            continue
        c = cells.iloc[a:b]
        d = E.gc_km(c.lat.values, c.lon.values, row.peak_lat, row.peak_lon)
        c = c[d <= 800]
        if len(c) == 0:
            continue
        R, C = E.disc_cells(row.peak_lat, row.peak_lon)
        oc = ocean[R, C]
        n_oc = int(oc.sum())
        if n_oc == 0:
            continue
        g = E.field(*row.fr)
        # ocean cells only: drop scatterometer cells whose nearest ERA5 cell is land
        ri = np.clip(np.rint((90 - c.lat.values) / 0.25).astype(int), 0, 360)
        ci = np.rint((c.lon.values % 360) / 0.25).astype(int) % 1440
        c = c[ocean[ri, ci]]
        ri, ci = ri[ocean[ri, ci]], ci[ocean[ri, ci]]
        if len(c) == 0:
            continue
        gv = E.bilinear(g, c.lat.values, c.lon.values)
        covered = len(set(zip(ri.tolist(), ci.tolist())))
        for sen, cs in c.assign(g=gv).groupby("sensor"):
            cov = len(set(zip(np.rint((90 - cs.lat.values) / 0.25).astype(int).tolist(),
                              (np.rint((cs.lon.values % 360) / 0.25).astype(int) % 1440).tolist()))) / n_oc
            foot.append(dict(track=row.track, basin=row.basin, season=row.season, peak_time=row.peak_time,
                             minp=row.minp, idx=row.gust800_kt, sensor=sen, n_cells=len(cs), cover=cov,
                             scat_max=cs.wind_kt.max(), era5_max_fp=cs.g.max(),
                             scat_p90=np.percentile(cs.wind_kt, 90), era5_p90_fp=np.percentile(cs.g, 90)))
            if with_pairs:
                pairs.append(cs.assign(track=row.track, basin=row.basin, season=row.season)
                             [["track", "basin", "season", "sensor", "lat", "lon", "wind_kt", "g"]])
    F = pd.DataFrame(foot)
    P = pd.concat(pairs) if pairs else pd.DataFrame()
    return F, P


if __name__ == "__main__":
    cells = pd.read_csv(sys.argv[1])
    F, P = collocate(cells)
    F["accepted"] = F.cover >= MIN_COVER
    F.to_csv(sys.argv[2] + "_footprints.csv", index=False)
    P.to_csv(sys.argv[2] + "_pairs.csv.gz", index=False)
    print(len(F), "footprints,", int(F.accepted.sum()), "with coverage >=", MIN_COVER)
