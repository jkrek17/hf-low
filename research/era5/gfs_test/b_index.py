"""Stage B, phase 1: GEFS reforecast Z500/U250 -> 5.625 degree window-mean anomaly maps (and the PR 41 index when its EOFs exist).
NO SCORING: nothing here reads archive counts or the ERA5 index.

Per initialisation (00 UTC on S-8) the raw GRIB messages saved by b_fetch.py (HGT 500 mb and UGRD 250 mb at forecast hours
12, 36, 60, 84, 108, 132, 156 = 12 UTC of valid days S-7 .. S-1) are decoded (0.5 degree, 361 x 720, lat 90 .. -90, lon 0 .. 359.5),
linearly interpolated onto the 0.25 degree 721 x 1440 grid that hemispheric/fields.py `regrid` expects, and regridded by that function
(area-weighted, 12 NH rows x 64 lon). Z500 in m (GRIB gh is gpm), U250 in m/s. The window mean is the plain mean of the seven 12 UTC steps
(no 00 UTC steps are used, as PR 41's fields are 12 UTC values).

Climatology: because every init is on the fixed calendar day Sep 23 + 7k and the same 7 leads, the model climatology for week k and
lead l is the mean over the (up to 15) seasons of the same-week same-lead forecast; as the window mean is linear it is stored as the
mean of the window-mean maps, one 12 x 64 map per (week k, variable), one value per season per cell, nothing smoothed. Leap days shift
the calendar date of a window by one day for k beyond the end of February in leap-year seasons; ignored. Anomaly = window map minus that mean.
No linear detrend (PR 41's detrend is a per-cell trend of Oct-Apr seasonal means over 1979-2021; the equivalent here would be a trend over
15 values per cell and would absorb part of the signal; the preregistration fixes the anomaly as the model-climatology anomaly only).

Index: needs the discovery-fit EOFs (mean, vectors, sd of the ERA5 weekly anomaly maps 2004-2014), which frozen_primary.json does NOT
carry (see the phase 1 report). `apply_frozen` below implements the application once an `eof` dict is supplied; without it the
index columns of gefs_index.csv are NaN and index_status says why.

usage: b_index.py [WORK] [OUT]        |        b_index.py --selftest
"""
import os, sys, types, json, glob, datetime as dt
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HEM = os.path.join(HERE, "..", "hemispheric")
HOURS = [12 + 24 * i for i in range(7)]
SEASONS = list(range(2004, 2019))
NWEEK = 30
K = 10
LAT = (-87.1875 + 5.625 * np.arange(32))[20:32]


def _regrid():
    """fields.regrid, imported with the cloud-IO modules stubbed (they are not installed and regrid does not use them)."""
    for m in ("numcodecs", "gcsfs"):
        if m not in sys.modules:
            try:
                __import__(m)
            except ImportError:
                sys.modules[m] = types.ModuleType(m)
    sys.path.insert(0, HEM)
    import fields
    return fields.regrid


def decode(path):
    import eccodes as ec
    with open(path, "rb") as f:
        g = ec.codes_grib_new_from_file(f)
        ni, nj = ec.codes_get(g, "Ni"), ec.codes_get(g, "Nj")
        assert (ni, nj) == (720, 361) and ec.codes_get(g, "jScansPositively") == 0
        assert ec.codes_get(g, "latitudeOfFirstGridPointInDegrees") == 90.0
        a = ec.codes_get_values(g).reshape(nj, ni).astype(np.float64)
        info = (ec.codes_get(g, "shortName"), ec.codes_get(g, "level"), ec.codes_get(g, "stepRange"), ec.codes_get(g, "validityTime"))
        ec.codes_release(g)
    return a, info


def to_quarter(a):
    """0.5 degree (361 x 720) -> 0.25 degree (721 x 1440), linear in lat and (periodic) lon; original points kept exactly."""
    lat = np.empty((721, a.shape[1])); lat[0::2] = a; lat[1::2] = 0.5 * (a[:-1] + a[1:])
    out = np.empty((721, 1440)); out[:, 0::2] = lat
    out[:, 1:-1:2] = 0.5 * (lat[:, :-1] + lat[:, 1:]); out[:, -1] = 0.5 * (lat[:, -1] + lat[:, 0])
    return out


def window_map(work, i0, var, regrid):
    acc = []
    for h in HOURS:
        a, (sn, lev, step, vt) = decode(os.path.join(work, "raw", f"{i0:%Y%m%d}", f"{var}_f{h:03d}.grb2"))
        assert (lev, int(step), vt) == ((500 if var == "z500" else 250), h, 1200), (lev, step, vt)
        acc.append(regrid(to_quarter(a)))
    return np.mean(acc, axis=0)


# --------------------------------------------------------------- the PR 41 index, once EOFs are available
def apply_frozen(eof, beta, z_anom, u_anom, sst_anom=None):
    """PR 41 `P` model pattern index = penalised-PC part of the linear predictor, X[:, nu:] @ beta[nu:].
    eof: dict(w=[12,64] sqrt(cos lat) weights, mean={z500,u250,sst}, vec={..}[10,N], sd={..}[10]) from the primary-split EOFs;
    beta: the 38 coefficients of models.P.beta in frozen_primary.json (7 month dummies, log1p(prev), then 10 z, 10 u, 10 sst PC coefficients);
    z_anom, u_anom: [n,12,64] anomaly maps. SST has no GEFS counterpart: if sst_anom is None the SST block is dropped (a change of the index)."""
    b = np.asarray(beta)[8:]
    out = 0.0
    for i, (v, M) in enumerate((("z500", z_anom), ("u250", u_anom))):
        X = (M * eof["w"]).reshape(len(M), -1)
        sc = ((X - eof["mean"][v]) @ eof["vec"][v].T) / eof["sd"][v]
        out = out + sc @ b[i * K:(i + 1) * K]
    return out


def selftest():
    """No network. (1) The PC-part formula applied to the committed primary-split PC scores (weekly_table.csv.gz) with the committed
    frozen betas gives the discovery-fit index; (2) it is compared with the committed leave-one-season-out oos_index files. These are
    different fits (22 fold-specific EOFs/lambdas versus one primary fit), so exact reproduction is NOT possible from committed files;
    the check reports the correlation and fails if the exact-reproduction criterion is not met."""
    R = os.path.join(HEM, "results")
    t = pd.read_csv(os.path.join(R, "weekly_table.csv.gz"))
    f = json.load(open(os.path.join(R, "frozen_primary.json")))
    cols = [f"{v}{i}" for v in "zus" for i in range(1, 11)]
    ok = True
    for b in ("atl", "pac"):
        be = np.array(f["basins"][b]["models"]["P"]["beta"])
        assert len(be) == 38
        d = t.assign(disc=t[cols].values @ be[8:])
        m = d.merge(pd.read_csv(os.path.join(R, f"oos_index_{b}.csv")), on=["season", "week"])
        r = np.corrcoef(m.disc, m.idx)[0, 1]
        mx = np.abs(m.disc - m.idx).max()
        print(f"{b}: discovery-fit index vs committed oos_index: corr {r:.3f}, max abs diff {mx:.3f} (n={len(m)})")
        ok &= bool(mx < 1e-3)
    print("EXACT REPRODUCTION:", "PASS" if ok else "FAIL (expected: oos_index is leave-one-season-out, not the discovery fit)")
    print("EOF vectors/mean/sd in frozen_primary.json: ABSENT; keys =", sorted(f.keys()), "->  maps cannot be projected")
    return ok


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else "/home/claude/gfs_work_b"
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    regrid = _regrid()
    rows, Z, U = [], [], []
    miss = {s: 0 for s in SEASONS}
    for s in SEASONS:
        for k in range(NWEEK):
            S = dt.date(s, 10, 1) + dt.timedelta(days=7 * k)
            i0 = S - dt.timedelta(days=8)
            done = os.path.exists(os.path.join(work, "raw", f"{i0:%Y%m%d}", "DONE"))
            rows.append(dict(season=s, week=k, init_date=str(i0), week_start=str(S), have=done))
            if done:
                Z.append(window_map(work, i0, "z500", regrid)); U.append(window_map(work, i0, "u250", regrid))
            else:
                miss[s] += 1
                Z.append(np.full((12, 64), np.nan)); U.append(np.full((12, 64), np.nan))
    Z, U = np.array(Z).reshape(len(SEASONS), NWEEK, 12, 64), np.array(U).reshape(len(SEASONS), NWEEK, 12, 64)
    Zc, Uc = np.nanmean(Z, axis=0), np.nanmean(U, axis=0)             # [30, 12, 64], mean of the available seasons
    nseas = np.isfinite(Z[:, :, 0, 0]).sum(0)
    Za, Ua = Z - Zc, U - Uc
    np.savez_compressed(os.path.join(out, "gefs_climatology.npz"), z500=Zc.astype(np.float32), u250=Uc.astype(np.float32),
                        n_seasons=nseas, lat=LAT)
    fn = os.path.join(out, "gefs_anom_maps.npz")
    np.savez_compressed(fn, z500=Za.reshape(-1, 12, 64).astype(np.float32), u250=Ua.reshape(-1, 12, 64).astype(np.float32),
                        season=np.repeat(SEASONS, NWEEK), week=np.tile(np.arange(NWEEK), len(SEASONS)), lat=LAT)
    if os.path.getsize(fn) > 10e6:
        np.savez_compressed(fn, z500=Za.reshape(-1, 12, 64).astype(np.float16), u250=Ua.reshape(-1, 12, 64).astype(np.float16),
                            season=np.repeat(SEASONS, NWEEK), week=np.tile(np.arange(NWEEK), len(SEASONS)), lat=LAT)
    df = pd.DataFrame(rows)
    df["atl_index"] = np.nan
    df["pac_index"] = np.nan
    df["index_status"] = np.where(df.have, "anomaly map built; index NOT computed: EOF vectors/mean/sd absent from frozen_primary.json",
                                  "init missing")
    df.drop(columns="have").to_csv(os.path.join(out, "gefs_index.csv"), index=False)
    print("rows", len(df), "missing inits per season", {s: n for s, n in miss.items() if n}, "total", sum(miss.values()))
    print("anom z500 std (m) over all cells", np.nanstd(Za), " u250 std (m/s)", np.nanstd(Ua))


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    main()
