"""HRRR-Alaska analysis (PREREGISTRATION_hrrr.md): H1 HRRR vs ERA5 zone gust, H2 distance of the peak-gust cell to land, H3 buoys, H4 rate table by regime and ingredient band.

HRRR values are a 3 hour MODEL forecast; ERA5 is a PROXY. Reads $ERA5_WORK/goa_jets/hrrr/*.npz (hrrr_pull.py), results/stage0_zone_gust.csv, results/buoy_era5_gust.csv and the NDBC files
in work/goa_jets/buoys. Writes results/hrrr_features.csv, hrrr_h1.csv, hrrr_h2.csv, hrrr_h3.csv, hrrr_h4.csv, hrrr_analysis.txt. rng seed 16.
usage: hrrr_analyse.py features | hrrr_analyse.py analyse
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RES = os.path.join(HERE, "results")
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "goa_jets", "hrrr")
R = 6371.0; KT = 1.943844
ZONES = {"Z1": (56, 61, -150, -138), "Z2": (53, 58, -140, -130), "Z3": (55, 60, -158, -148)}
rng = np.random.default_rng(16)
out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))


def xyz(la, lo):
    la, lo = np.radians(la), np.radians(lo)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]) * R


def regime(zn, wd):
    if zn in ("Z1", "Z2"):
        if 100 <= wd < 170: return "BARRIER"
        if 170 <= wd < 260: return "ONSHORE"
        return "OFFSHORE"
    if wd >= 300 or wd < 60: return "GAP"
    if 90 <= wd < 200: return "INFLOW"
    return "OTHER"


def season(t):
    return t.year if t.month >= 10 else t.year - 1


def load_grid():
    c = np.load(os.path.join(WORK, "crop.npz")); la, lo = c["lat"].astype(float), c["lon"].astype(float)
    land = np.load(os.path.join(WORK, "land.npy")).astype(float) >= 0.5
    return la, lo, land


def dland_grid(la, lo, land):
    tree = cKDTree(xyz(la[land], lo[land])); d = np.full(la.shape, np.nan)
    d[~land] = tree.query(xyz(la[~land], lo[~land]))[0]
    return d, tree


def features():
    la, lo, land = load_grid(); ocean = ~land; dl, _ = dland_grid(la, lo, land)
    box = lambda s, n, w, e: (la >= s) & (la <= n) & (lo >= w) & (lo <= e)
    ZM = {z: box(*b) & ocean for z, b in ZONES.items()}
    ZC = {z: box(44, 66, b[2], b[3]) & ocean for z, b in ZONES.items()}
    M_INT = box(62, 66, -155, -145); M_EXIT = box(56, 59, -156, -150) & ocean
    rows = []; miss = 0
    from hrrr_pull import vtimes
    for t in vtimes():
        fn = os.path.join(WORK, f"{t:%Y%m%d%H}.npz")
        if not os.path.exists(fn): miss += 1; continue
        z = np.load(fn); g = z["gust"].astype(float); u = z["u"].astype(float); v = z["v"].astype(float); p = z["msl"].astype(float) / 10.0   # hPa
        interior_p = p[M_INT].mean(); exit_p = p[M_EXIT].mean()
        for zn in ZONES:
            m = ZM[zn]; gz = np.where(m, g, -1.0); k = np.unravel_index(np.argmax(gz), gz.shape)
            um, vm = u[m].mean(), v[m].mean(); wdir = np.degrees(np.arctan2(-um, -vm)) % 360
            pile = np.nan
            if zn in ("Z1", "Z2"):
                a = p[ZC[zn] & (dl <= 100)]; b = p[ZC[zn] & (dl >= 300) & (dl <= 400)]
                pile = a.mean() - b.mean() if len(a) and len(b) else np.nan
            rows.append(dict(time=t, zone=zn, hrrr_gust=float(gz[k]), glat=la[k], glon=lo[k], g_dland_km=dl[k], wmax_kt=float(np.hypot(u, v)[m].max()),
                             wspd_mean_kt=float(np.hypot(um, vm)), wdir=wdir, COASTPILE=pile, INTERIOR_DP=(interior_p - exit_p) if zn == "Z3" else np.nan))
    df = pd.DataFrame(rows); df.to_csv(os.path.join(RES, "hrrr_features.csv"), index=False)
    print(len(df), "rows;", miss, "valid times missing")


def boot_ci(f, seasons_data, n=2000):
    ss = list(seasons_data); r = []
    for _ in range(n):
        r.append(f(pd.concat([seasons_data[s] for s in rng.choice(ss, len(ss))])))
    r = np.array(r, float); return np.nanpercentile(r, [2.5, 97.5])


def analyse():
    f = pd.read_csv(os.path.join(RES, "hrrr_features.csv"), parse_dates=["time"])
    s0 = pd.read_csv(os.path.join(RES, "stage0_zone_gust.csv"), parse_dates=["time"])
    e = pd.concat([s0[["time", f"{z}_gust_kt", f"{z}_lat", f"{z}_lon"]].set_axis(["time", "era5_gust", "era5_lat", "era5_lon"], axis=1).assign(zone=z) for z in ZONES])
    f = f.merge(e, on=["time", "zone"], how="left"); assert f.era5_gust.notna().all()
    f["regime"] = [regime(z, w) for z, w in zip(f.zone, f.wdir)]; f["season"] = [season(t) for t in f.time]
    for thr in (54, 64, 71.7): f[f"h{int(thr)}"] = f.hrrr_gust >= thr; f[f"e{int(thr)}"] = f.era5_gust >= thr
    NT = f.time.nunique(); P("valid times with a field:", NT, "of 3396; per season:", f.groupby("season").time.nunique().to_dict())
    # ---- H1
    rows = []
    for zn in ZONES:
        d = f[f.zone == zn]
        for rg in ["(all)"] + sorted(d.regime.unique()):
            x = d if rg == "(all)" else d[d.regime == rg]
            if len(x) < 30: continue
            sd = {s: y for s, y in x.groupby("season")}
            ratio_f = lambda y: float((y.hrrr_gust / y.era5_gust).median())
            ci = boot_ci(ratio_f, sd, 1000)
            row = dict(zone=zn, regime=rg, n=len(x), spearman=float(spearmanr(x.hrrr_gust, x.era5_gust)[0]), median_ratio=ratio_f(x), ratio_lo=ci[0], ratio_hi=ci[1],
                       hrrr_gust_median=float(x.hrrr_gust.median()), era5_gust_median=float(x.era5_gust.median()))
            for thr in (54, 64):
                a, b = x[f"h{thr}"], x[f"e{thr}"]
                row.update({f"both{thr}": int((a & b).sum()), f"hrrr_only{thr}": int((a & ~b).sum()), f"era5_only{thr}": int((~a & b).sum()), f"hrrr_total{thr}": int(a.sum()), f"era5_total{thr}": int(b.sum())})
            rows.append(row)
    H1 = pd.DataFrame(rows); H1.to_csv(os.path.join(RES, "hrrr_h1.csv"), index=False)
    P("\nH1 HRRR vs ERA5 zone-maximum gust at the same 00/12 UTC times (regime from HRRR wind). ratio = median HRRR/ERA5; event counts for 54 and 64 kt")
    P(H1.round(3).to_string(index=False))
    # ---- H2
    la, lo, land = load_grid(); dl, tree = dland_grid(la, lo, land)
    f["era5_dland_km"] = tree.query(xyz(f.era5_lat.values, f.era5_lon.values))[0]
    rows = []
    for zn, rg in (("Z1", "BARRIER"), ("Z2", "BARRIER"), ("Z3", "GAP")):
        d = f[(f.zone == zn) & (f.regime == rg)]
        for lab, ev, col in (("HRRR events (54 kt+), HRRR peak cell", d.h54, "g_dland_km"), ("ERA5 events (54 kt+), ERA5 peak cell", d.e54, "era5_dland_km")):
            x = d[ev]
            if len(x) < 5: continue
            rows.append(dict(zone=zn, regime=rg, subset=lab, n=len(x), median_km=float(x[col].median()), share_within_10=float((x[col] <= 10).mean()), share_within_25=float((x[col] <= 25).mean()),
                             share_within_50=float((x[col] <= 50).mean())))
    H2 = pd.DataFrame(rows); H2.to_csv(os.path.join(RES, "hrrr_h2.csv"), index=False)
    P("\nH2 distance from the peak-gust cell to the nearest land (km), events in the barrier / gap regime (HRRR land mask for both)"); P(H2.round(3).to_string(index=False))
    # ---- H3 buoys
    import buoy as B
    c = np.load(os.path.join(WORK, "crop.npz")); lat, lon = c["lat"].astype(float), c["lon"].astype(float)
    eg = pd.read_csv(os.path.join(RES, "buoy_era5_gust.csv"), parse_dates=["time"]); eg["station"] = eg.station.astype(str)
    egd = {(a, b): c for a, b, c in zip(eg.time, eg.station, eg.era5_gust_kt)}
    hz = {(a, b): c for a, b, c in zip(f.time, f.zone, f.hrrr_gust)}
    zone_of = lambda la_, lo_: next((z for z in ("Z1", "Z3", "Z2") if ZONES[z][0] <= la_ <= ZONES[z][1] and ZONES[z][2] <= lo_ <= ZONES[z][3]), None)
    times = sorted(f.time.unique()); fields = {}
    rows = []
    for s, (bla, blo) in B.ST.items():
        j, i = np.unravel_index(np.argmin((lat - bla) ** 2 + ((lon - blo) * np.cos(np.radians(bla))) ** 2), lat.shape)
        dist = float(np.hypot(lat[j, i] - bla, (lon[j, i] - blo) * np.cos(np.radians(bla))) * 111.19)
        lf = float(land[j, i]); b = B.read_buoy(s); zn = zone_of(bla, blo)
        P(f"buoy {s}: nearest HRRR cell {lat[j, i]:.3f}N {lon[j, i]:.3f}E is {dist:.1f} km away, land={bool(lf)}; zone {zn}")
        for t in times:
            t = pd.Timestamp(t); w = b[(b.t - t).abs() <= pd.Timedelta(minutes=40)]
            if w.empty or np.isnan(w.gs).all(): continue
            w = w.iloc[(w.t - t).abs().argmin()]
            z = fields.setdefault(t, np.load(os.path.join(WORK, f"{t:%Y%m%d%H}.npz")))
            hg = float(z["gust"][j, i]); hw = float(np.hypot(z["u"][j, i], z["v"][j, i]))
            rows.append(dict(time=t, station=s, hrrr_cell_land=bool(lf), cell_km=dist, buoy_ws_kt=w.ws * KT, buoy_gust_kt=w.gs * KT, hrrr_wind_kt=hw, hrrr_gust_kt=hg,
                             era5_gust_kt=egd.get((t, s), np.nan), hrrr_zone_gust=hz.get((t, zn), np.nan)))
        fields.clear()
    BP = pd.DataFrame(rows); BP.to_csv(os.path.join(RES, "hrrr_buoy_pairs.csv"), index=False)
    out3 = []
    for s, d in BP.groupby("station"):
        for lab, m in (("buoy gust >= 35 kt", d.buoy_gust_kt >= 35), ("buoy gust >= 45 kt", d.buoy_gust_kt >= 45)):
            x = d[m]
            if len(x) < 5: continue
            q = lambda a: np.nanpercentile(a, [25, 50, 75])
            out3.append(dict(station=s, subset=lab, n=len(x), hrrr_cell_land=bool(d.hrrr_cell_land.iloc[0]), cell_km=float(d.cell_km.iloc[0]),
                             gust_ratio_med=q(x.hrrr_gust_kt / x.buoy_gust_kt)[1], gust_ratio_q25=q(x.hrrr_gust_kt / x.buoy_gust_kt)[0], gust_ratio_q75=q(x.hrrr_gust_kt / x.buoy_gust_kt)[2],
                             wind_ratio_med=q(x.hrrr_wind_kt / x.buoy_ws_kt)[1], era5_gust_ratio_med=float(np.nanmedian(x.era5_gust_kt / x.buoy_gust_kt)),
                             share_hrrr_zone_ge54=float((x.hrrr_zone_gust >= 54).mean())))
    H3 = pd.DataFrame(out3); H3.to_csv(os.path.join(RES, "hrrr_h3.csv"), index=False)
    P("\nH3 buoys: ratio HRRR (nearest cell) / buoy at times chosen on the buoy gust; ERA5 gust ratio on the same times; share of those times with HRRR zone gust >= 54 kt")
    P(H3.round(3).to_string(index=False))
    # ---- H4
    BANDS = [("bottom half", 0, 50), ("50th-80th", 50, 80), ("80th-95th", 80, 95), ("top 5%", 95, 100.001)]
    rows = []
    for zn in ZONES:
        d = f[f.zone == zn]; ing = "INTERIOR_DP" if zn == "Z3" else "COASTPILE"
        for rg in sorted(d.regime.unique()):
            dr = d[d.regime == rg]
            sets = [("(all times in regime)", "all", dr)]
            if len(dr) >= 80:
                cuts = np.nanpercentile(dr[ing], [50, 80, 95])
                for bn, lo_, hi_ in BANDS:
                    qlo = -np.inf if lo_ == 0 else cuts[[50, 80, 95].index(lo_)]; qhi = np.inf if hi_ > 100 else cuts[[50, 80, 95].index(hi_)]
                    sets.append((ing, bn, dr[(dr[ing] >= qlo) & (dr[ing] < qhi)]))
            for ingn, bn, x in sets:
                if len(x) == 0: continue
                sd = {s: y for s, y in x.groupby("season")}
                ci = boot_ci(lambda y: 100 * y.h54.mean(), sd, 1000); ps = 100 * x.groupby("season").h54.mean()
                rows.append(dict(zone=zn, regime=rg, ingredient=ingn, band=bn, n=len(x), hrrr_events=int(x.h54.sum()), hrrr_strong=int(x.h64.sum()), hrrr_hf=int(x.h71.sum()),
                                 pct_hrrr54=100 * x.h54.mean(), lo=ci[0], hi=ci[1], pct_hrrr64=100 * x.h64.mean(), pct_era5_54_same_times=100 * x.e54.mean(), pct_era5_64_same_times=100 * x.e64.mean(),
                                 season_min=ps.min(), season_max=ps.max()))
    H4 = pd.DataFrame(rows); H4.to_csv(os.path.join(RES, "hrrr_h4.csv"), index=False)
    P("\nH4 chance (%) that a 00/12 UTC time has an HRRR zone gust of 54 kt or more (64 kt in pct_hrrr64), with the ERA5 rate at the same times; exact over all times, season-block interval")
    P(H4[H4.band == "all"].round(2).to_string(index=False))
    for zn in ZONES:
        ing = "INTERIOR_DP" if zn == "Z3" else "COASTPILE"; rg = "GAP" if zn == "Z3" else "BARRIER"
        P(f"\nH4 {zn} {rg} by band of {ing}:"); P(H4[(H4.zone == zn) & (H4.regime == rg) & (H4.ingredient == ing)][["band", "n", "hrrr_events", "pct_hrrr54", "lo", "hi", "pct_hrrr64", "pct_era5_54_same_times", "season_min", "season_max"]].round(1).to_string(index=False))
    open(os.path.join(RES, "hrrr_analysis.txt"), "w").write("\n".join(out)); print("\n".join(out))


if __name__ == "__main__":
    features() if sys.argv[1] == "features" else analyse()
