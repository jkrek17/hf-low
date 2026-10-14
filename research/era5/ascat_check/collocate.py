"""Collocate scatterometer cells with ERA5 gust around pipeline A storm fixes.

Input scatterometer cell table (any source adapter must write this): csv/csv.gz with columns
  sensor, time (UTC, ISO) or t (hours since 1900), lat, lon (-180..360 ok), wind_kt, ok (1 = passes rain/ambiguity/quality flags)
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


def scat_ocean_grid():
    """720x1440 scatterometer-grid ocean mask (cell centres lat -89.875+0.25i, lon 0.125+0.25j; True where the mean of the four
    surrounding ERA5 land-sea-mask nodes is < 0.5; southern hemisphere False). Fixes the half-cell offset between the two grids."""
    lsm = np.load(f"{E.CACHE}/lsm_nh_0p25.npy").astype(np.float32)            # 361 x 1440, rows 90N..0N, cols 0..359.75
    k = 719 - np.arange(360, 720)                                              # scat row i (NH) -> upper ERA5 row
    m = np.arange(1440)
    mean = (lsm[k][:, m] + lsm[k + 1][:, m] + lsm[k][:, (m + 1) % 1440] + lsm[k + 1][:, (m + 1) % 1440]) / 4.0
    out = np.zeros((720, 1440), bool)
    out[360:] = mean < 0.5
    return out


def disc_scat_ocean(lat, lon, ocs, radius=800.0):
    """number of scatterometer-grid ocean cells with centres within radius km of (lat, lon)."""
    la = -89.875 + 0.25 * np.arange(720)
    lo = 0.125 + 0.25 * np.arange(1440)
    dlat = np.degrees(radius / E.RE)
    i = np.nonzero((la >= lat - dlat) & (la <= lat + dlat))[0]
    dlon = min(180.0, dlat / max(np.cos(np.radians(min(abs(lat) + dlat, 89))), 0.05))
    j0 = int(round(((lon % 360) - 0.125) / 0.25))
    n = int(np.ceil(dlon / 0.25)) + 1
    j = (j0 + np.arange(-n, n + 1)) % 1440
    LA, LO = np.meshgrid(la[i], lo[j], indexing="ij")
    d = E.gc_km(LA, LO, lat, lon)
    return int(((d <= radius) & ocs[np.ix_(i, j)]).sum())


def collocate(cells, S=None, ocean=None, idx=None, with_pairs=True):
    S = storm_sample() if S is None else S
    ocs = scat_ocean_grid()
    idx = E.cache_index() if idx is None else idx
    S = S[S.h1900.isin(idx)].copy()
    S["fr"] = [idx[h] for h in S.h1900]
    S = S.sort_values("fr")
    cells = cells[cells.ok == 1].copy()
    if "t" not in cells:                       # adapters may supply t (hours since 1900) directly: rss_adapter.load_cells does
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
        n_oc = disc_scat_ocean(row.peak_lat, row.peak_lon, ocs)
        if n_oc == 0:
            continue
        g = E.field(*row.fr)
        # ocean cells only: the scatterometer cell (0.25 deg, centres on x.125) must be ocean in the 4 surrounding ERA5 nodes' mean
        si = np.rint((c.lat.values + 89.875) / 0.25).astype(int)
        sj = np.rint(((c.lon.values % 360) - 0.125) / 0.25).astype(int) % 1440
        keep = ocs[si, sj]
        c, si, sj = c[keep], si[keep], sj[keep]
        if len(c) == 0:
            continue
        gv = E.bilinear(g, c.lat.values, c.lon.values)
        for sen, cs in c.assign(g=gv).groupby("sensor"):
            cov = len(set(zip(np.rint((cs.lat.values + 89.875) / 0.25).astype(int).tolist(),
                              (np.rint(((cs.lon.values % 360) - 0.125) / 0.25).astype(int) % 1440).tolist()))) / n_oc
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
