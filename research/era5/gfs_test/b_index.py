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

Index (phase 1b): the primary-split (discovery seasons 2004-2014, lag-1 window S-7..S-1) EOFs of PR 41 are REBUILT deterministically from
WeatherBench2 Z500/U250 (`pull_era5`, `rebuild_eofs`; hemispheric/hemlib.py is called unchanged, with hemlib.VARS narrowed to z500, u250 --
SST/MSLP are not needed for the z/u EOFs, which are computed per variable) and checked against the committed weekly_table.csv.gz scores
(`eof_check`). `index_zu` = frozen betas of models.P (basins.<b>.models.P.beta, lambda 1000; layout: 7 month dummies Oct..Apr, log1p(prev), 10 z500
PC, 10 u250 PC, 10 SST PC = 38) beta[8:18] . z-PC + beta[18:28] . u-PC, PCs being unit-variance (discovery sd) scores of the ERA5 or GEFS anomaly
maps; NOT standardised further, no intercept/month/prev term, SST block dropped (a different quantity from the full PR 41 pattern index).
GEFS maps are projected unchanged (model-climatology anomalies, no detrend). Nothing here reads archive counts; no score is computed.

usage: b_index.py [GEFSWORK] [OUT]      GEFS anomaly maps (+ index_zu.csv when the EOF file exists)
       b_index.py --pull ERA5WORK       fetch the needed WeatherBench2 chunks (needs gcsfs, numcodecs; 3 threads)
       b_index.py --eofs ERA5WORK       rebuild pr41_eofs_zu.npz, check against weekly_table.csv.gz, write index_zu.csv
       b_index.py --selftest
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


# --------------------------------------------------------------- ERA5 (WeatherBench2) pull and EOF rebuild
EOFS = os.path.join(HERE, "results", "pr41_eofs_zu.npz")
WB2_T0 = dt.datetime(1959, 1, 1)


def season_dates(s):
    """Days of season s needed. 1991..2019 (climatology seasons): hemlib's full 5 Sep .. 10 May (the 15-day running mean and its end
    truncation are defined on that calendar). Other trend seasons (1979-90, 2020-21): 1 Oct .. 30 Apr (the trend uses Oct-Apr means only).
    2004-2018 additionally reach back to 24 Sep (lag-1 window of week 0)."""
    if 1991 <= s <= 2019:
        a, b = dt.date(s, 9, 5), dt.date(s + 1, 5, 10)
    else:
        a, b = dt.date(s, 10, 1), dt.date(s + 1, 4, 30)
    if 2004 <= s <= 2018:
        a = min(a, dt.date(s, 9, 24))
    return [a + dt.timedelta(i) for i in range((b - a).days + 1)]


def _wb2_index(d):
    return (dt.datetime(d.year, d.month, d.day) - WB2_T0).days * 4 + 2        # 12 UTC step of the 6-hourly store


def needed_chunks():
    return sorted({_wb2_index(d) // 100 for s in range(1979, 2022) for d in season_dates(s)})


def pull_era5(work, nth=3):
    """Fetch Z500 and U250 WB2 chunks (100 six-hourly steps each, all 13 levels are inside the object, 7-10 MB) with hemispheric/fields.wb2_chunk
    (needs gcsfs + numcodecs); store the (100, 12, 64) NH arrays as work/chunks/<var>_<c>.npy; bytes in work/chunks/bytes.log."""
    from concurrent.futures import ThreadPoolExecutor
    sys.path.insert(0, HEM)
    import fields
    os.makedirs(os.path.join(work, "chunks"), exist_ok=True)

    def job(vc):
        v, c = vc
        fn = os.path.join(work, "chunks", f"{v}_{c}.npy")
        if os.path.exists(fn):
            return 0
        a, n = fields.wb2_chunk(v, c)
        np.save(fn + ".tmp.npy", a); os.rename(fn + ".tmp.npy", fn)
        with open(os.path.join(work, "chunks", "bytes.log"), "a") as f:
            f.write(f"{v} {c} {n}\n")
        return n
    tot = 0
    with ThreadPoolExecutor(nth) as ex:
        for n in ex.map(job, [(v, c) for c in needed_chunks() for v in ("z500", "u250")]):
            tot += n
    print(f"pulled {tot/1e9:.3f} GB this run")


def era5_weekly_maps(work):
    """hemlib anomalies (1991-2020 calendar climatology with 15-day running mean, per-cell linear trend of Oct-Apr means 1979-2021 removed)
    and lag-1 window means for seasons 2004-2018. Returns (maps dict z500/u250 [15*30,12,64], seasons)."""
    sys.path.insert(0, HEM)
    import hemlib as H
    cache = {}

    def get(v, c):
        if (v, c) not in cache:
            cache[(v, c)] = np.load(os.path.join(work, "chunks", f"{v}_{c}.npy"))
        return cache[(v, c)]
    F = {}
    for s in range(1979, 2022):
        ds = season_dates(s)
        idx = [_wb2_index(d) for d in ds]
        o = {v: np.stack([get(v, i // 100)[i % 100] for i in idx]) for v in ("z500", "u250")}
        o["mslp"] = np.zeros_like(o["z500"]); o["sst"] = np.zeros_like(o["z500"])      # unused placeholders (hemlib loops over all four)
        F[s] = (ds, o)
    A = H.anomalies(F)
    old = H.VARS
    H.VARS = ("z500", "u250")
    try:
        M = H.weekly_maps(F, A, SEASONS, -7, -1)
    finally:
        H.VARS = old
    return M


def rebuild_eofs(M, path=EOFS):
    """Primary split: EOFs of the 2004-2014 maps (first 11*30 rows), exactly hemlib.EOFs, z500 and u250 only."""
    sys.path.insert(0, HEM)
    import hemlib as H
    ntr = 11 * NWEEK
    old = H.VARS
    H.VARS = ("z500", "u250")
    try:
        eof = H.EOFs({v: M[v][:ntr] for v in H.VARS}, None)
        sc = eof.scores({v: M[v][:ntr] for v in H.VARS})
    finally:
        H.VARS = old
    np.savez_compressed(path, w=eof.w, **{f"{k}_{v}": getattr(eof, k)[v] for k in ("mean", "vec", "sd") for v in ("z500", "u250")})
    return eof, sc


def load_eofs(path=EOFS):
    z = np.load(path)
    return dict(w=z["w"], mean={v: z[f"mean_{v}"] for v in ("z500", "u250")}, vec={v: z[f"vec_{v}"] for v in ("z500", "u250")},
                sd={v: z[f"sd_{v}"] for v in ("z500", "u250")})


def eof_check(sc):
    """Rebuilt discovery PC scores vs committed weekly_table.csv.gz (seasons 2004-2014). Returns per-PC correlation (z1..z10,u1..u10) and max |diff| after sign."""
    t = pd.read_csv(os.path.join(HEM, "results", "weekly_table.csv.gz"))
    t = t[(t.season >= 2004) & (t.season <= 2014)].sort_values(["season", "week"])
    r, d = [], []
    for i, v in enumerate("zu"):
        for j in range(K):
            c = t[f"{v}{j+1}"].values
            rr = np.corrcoef(sc[:, i * K + j], c)[0, 1]
            r.append(rr); d.append(np.abs(sc[:, i * K + j] - np.sign(rr) * c).max())
    return np.array(r), np.array(d)


def betas():
    f = json.load(open(os.path.join(HEM, "results", "frozen_primary.json")))
    out = {}
    for b in ("atl", "pac"):
        be = np.array(f["basins"][b]["models"]["P"]["beta"])
        assert len(be) == 38
        out[b] = be
    return out


def index_zu_table(eof, Mera5, gefs_npz):
    """Rows season x week 2004-2018 (era5 has all; GEFS NaN where an init is missing)."""
    be = betas()
    g = np.load(gefs_npz)
    gz, gu = g["z500"].astype(float), g["u250"].astype(float)
    assert (g["season"] == np.repeat(SEASONS, NWEEK)).all() and (g["week"] == np.tile(np.arange(NWEEK), len(SEASONS))).all()
    df = pd.DataFrame(dict(season=np.repeat(SEASONS, NWEEK), week=np.tile(np.arange(NWEEK), len(SEASONS))))
    ok = np.isfinite(gz).all((1, 2)) & np.isfinite(gu).all((1, 2))
    for b in ("atl", "pac"):
        df[f"{b}_era5"] = apply_frozen(eof, be[b], Mera5["z500"].astype(float), Mera5["u250"].astype(float))
        gi = np.full(len(df), np.nan)
        gi[ok] = apply_frozen(eof, be[b], gz[ok], gu[ok])
        df[f"{b}_gefs"] = gi
    return df[["season", "week", "atl_era5", "pac_era5", "atl_gefs", "pac_gefs"]]


def run_eofs(work):
    M = era5_weekly_maps(work)
    eof, sc = rebuild_eofs(M)
    r, d = eof_check(sc)
    names = [f"{v}{j}" for v in "zu" for j in range(1, 11)]
    print("EOF check vs weekly_table.csv.gz (2004-2014, n=330): min corr", r.min(), " max |diff| after sign", d.max())
    print("per-PC |corr|:", dict(zip(names, np.round(np.abs(r), 6))))
    fz = json.load(open(os.path.join(HEM, "results", "frozen_primary.json")))["eof_variance_explained"]
    print("variance explained (rebuilt vs committed):", {v: round(float(eof.var_explained[v]), 4) for v in ("z500", "u250")}, fz)
    if abs(r).min() < 0.999:
        print("EOF CHECK FAILED; index_zu.csv NOT written"); return False
    df = index_zu_table(load_eofs(), M, os.path.join(HERE, "results", "gefs_anom_maps.npz"))
    df.to_csv(os.path.join(HERE, "results", "index_zu.csv"), index=False, float_format="%.6g")
    print("wrote index_zu.csv", df.shape, "GEFS rows present", int(df.atl_gefs.notna().sum()))
    return True


def selftest():
    """No network. (1) the committed pr41_eofs_zu.npz is consistent (orthonormal vectors, sd > 0, weights = sqrt(cos lat)); (2) apply_frozen equals the
    explicit PC formula on random maps; (3) the discovery-fit index from the committed PC scores (weekly_table.csv.gz) and the frozen betas is finite and
    reproduces what apply_frozen would give on the same PCs (beta[8:28] on z,u PCs). The leave-one-season-out oos_index files are a different fit
    and are NOT expected to match (that earlier criterion is dropped; see PREREGISTRATION deviation 5). Full EOF reproduction is `--eofs`."""
    ok = True
    e = load_eofs()
    for v in ("z500", "u250"):
        G = e["vec"][v] @ e["vec"][v].T
        ok &= bool(np.allclose(G, np.eye(K), atol=1e-8)) and bool((e["sd"][v] > 0).all())
    ok &= bool(np.allclose(e["w"], np.sqrt(np.cos(np.radians(LAT)))[:, None] * np.ones((1, 64))))
    rng = np.random.default_rng(1)
    z, u = rng.normal(size=(5, 12, 64)) * 50, rng.normal(size=(5, 12, 64)) * 5
    be = betas()["atl"]
    man = 0
    for i, (v, M) in enumerate((("z500", z), ("u250", u))):
        sc = ((M * e["w"]).reshape(5, -1) - e["mean"][v]) @ e["vec"][v].T / e["sd"][v]
        man = man + sc @ be[8 + 10 * i: 18 + 10 * i]
    ok &= bool(np.allclose(man, apply_frozen(e, be, z, u)))
    t = pd.read_csv(os.path.join(HEM, "results", "weekly_table.csv.gz"))
    cols = [f"{v}{i}" for v in "zu" for i in range(1, 11)]
    d = t[cols].values @ be[8:28]
    ok &= bool(np.isfinite(d).all())
    print("SELFTEST:", "PASS" if ok else "FAIL")
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
    df["index_status"] = np.where(df.have, "anomaly map built", "init missing")
    if os.path.exists(EOFS):
        be, eof = betas(), load_eofs()
        okm = np.isfinite(Z).all((2, 3)).ravel() & np.isfinite(U).all((2, 3)).ravel()
        for b in ("atl", "pac"):
            v = np.full(len(df), np.nan)
            v[okm] = apply_frozen(eof, be[b], Za.reshape(-1, 12, 64)[okm].astype(float), Ua.reshape(-1, 12, 64)[okm].astype(float))
            df[f"{b}_index"] = v
        df["index_status"] = np.where(df.have, "z500+u250 index_zu (frozen models.P betas, SST block dropped)", "init missing")
    df.drop(columns="have").to_csv(os.path.join(out, "gefs_index.csv"), index=False)
    print("rows", len(df), "missing inits per season", {s: n for s, n in miss.items() if n}, "total", sum(miss.values()))
    print("anom z500 std (m) over all cells", np.nanstd(Za), " u250 std (m/s)", np.nanstd(Ua))


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    if "--pull" in sys.argv:
        pull_era5(sys.argv[sys.argv.index("--pull") + 1]); sys.exit()
    if "--eofs" in sys.argv:
        sys.exit(0 if run_eofs(sys.argv[sys.argv.index("--eofs") + 1]) else 1)
    main()
