"""Descriptive (non-verdict) outputs of the ASCAT check, plus labelled post hoc checks.
usage: descriptives.py footprints.csv pairs.csv.gz cells_dir"""
import sys, numpy as np, pandas as pd
import analyse as A, collocate as C, rss_adapter as R, era5side as E

F = pd.read_csv(sys.argv[1]); F = F[F.cover >= 0.40]
F1 = F.sort_values("cover", ascending=False).drop_duplicates("track")            # the analysis set (one pass per storm)
print("== analysis set (one best-coverage pass per storm, cover >= 0.40)")
print(F1.groupby("basin").agg(storms=("track", "size"), seasons=("season", "nunique"), cover_mean=("cover", "mean"),
                              hf_index=("idx", lambda x: int((x >= 71.7).sum()))).to_string())
print("sensor of chosen pass:", F1.sensor.value_counts().to_dict(), " first/last season:", int(F1.season.min()), int(F1.season.max()))
print("\n== footprint maxima (kt): ERA5 vs scatterometer, by ERA5 index class")
F1 = F1.assign(cls=pd.cut(F1.idx, [0, 54, 64, 71.7, 200], right=False, labels=["<54", "54-64", "64-71.7", ">=71.7"]))
print(F1.groupby("cls", observed=True).agg(n=("track", "size"), era5_fp_med=("era5_max_fp", "median"), scat_med=("scat_max", "median"),
      scat_p90=("scat_max", lambda x: np.percentile(x, 90)), share_scat_ge45=("scat_max", lambda x: (x >= 45).mean()),
      share_scat_ge50=("scat_max", lambda x: (x >= 50).mean())).round(2).to_string())
print("\n== matched-quantile conversion (descriptive): scatterometer max at the quantile where ERA5 footprint max = x kt")
for b, g in list(F1.groupby("basin")) + [("both", F1)]:
    e = np.sort(g.era5_max_fp.values); s = np.sort(g.scat_max.values)
    for x in (54, 64, 71.7):
        q = (e < x).mean()
        print(f"{b:4s} ERA5 fp max {x:5.1f} kt = quantile {q:.3f} -> scat {np.quantile(s, min(q, 1)):.1f} kt (ratio {np.quantile(s, min(q, 1)) / x:.2f}), n={len(g)}")
print("\n== P4 effect sizes: ratio scat/ERA5 footprint max, ERA5 fp max >= 54 kt, all accepted passes, season-block 95% CI")
G = F.copy(); G = G[G.era5_max_fp >= 54]; G["lr"] = np.log(G.scat_max / G.era5_max_fp)
rng = np.random.default_rng(7)
for s, g in G.groupby("sensor"):
    bs = A.boot_seasons(g, lambda d: float(np.exp(d.lr.mean())), 500)
    print(f"{s}: n={len(g)} ratio {np.exp(g.lr.mean()):.3f} [{np.percentile(bs, 2.5):.3f}, {np.percentile(bs, 97.5):.3f}]  seasons {g.season.min()}-{g.season.max()}")
print("ratio by season (all sensors):", G.groupby("season").lr.agg(lambda x: round(float(np.exp(x.mean())), 3)).to_dict())
print("slope of ratio on season (per decade):", round(float(np.polyfit(G.season, G.lr, 1)[0] * 10), 4), "log units  (descriptive)")
print("\n== pointwise (cell pair) comparison on the analysis set")
P = pd.read_csv(sys.argv[2])
P = P.merge(F1[["track", "sensor"]], on=["track", "sensor"])
S = C.storm_sample()[["track", "peak_lat", "peak_lon"]].drop_duplicates("track")
P = P.merge(S, on="track")
P["dist"] = E.gc_km(P.lat.values, P.lon.values, P.peak_lat.values, P.peak_lon.values)
print("cell pairs", len(P), "pearson", round(float(np.corrcoef(P.g, P.wind_kt)[0, 1]), 3), "slope wind on gust", round(float(np.polyfit(P.g, P.wind_kt, 1)[0]), 3),
      "mean wind/gust", round(float(P.wind_kt.mean() / P.g.mean()), 3))
P["latband"] = pd.cut(P.lat, [0, 40, 50, 60, 90]); P["dband"] = pd.cut(P.dist, [0, 200, 400, 600, 800])
print(P.groupby("latband", observed=True).apply(lambda d: pd.Series(dict(n=len(d), ratio=d.wind_kt.mean() / d.g.mean(), bias_kt=(d.wind_kt - d.g).mean()))).round(3).to_string())
print(P.groupby("dband", observed=True).apply(lambda d: pd.Series(dict(n=len(d), ratio=d.wind_kt.mean() / d.g.mean(), bias_kt=(d.wind_kt - d.g).mean()))).round(3).to_string())
print("\n== rain-flag losses (cells within 800 km and +-1.5 h of a sample fix)")
cells = R.load_cells(sys.argv[3], "scat"); cells2 = R.load_cells(sys.argv[3], "scat+rad")
atl = ((cells.lon % 360) > 250) | ((cells.lon % 360) < 60)
print("scatterometer rain flag: dropped", round(1 - cells.ok.mean(), 4), " Atlantic-side", round(1 - cells.ok[atl].mean(), 4), " Pacific-side", round(1 - cells.ok[~atl].mean(), 4),
      "(basin by cell longitude, approximate)")
print("scatterometer + radiometer rain: dropped", round(1 - cells2.ok.mean(), 4))
print("\n== POST HOC (not registered): robustness of P1 to footprint size")
for b, g in F1.groupby("basin"):
    g = g.assign(cb=pd.qcut(g.n_cells, 3, labels=False))
    print(b, "rho within n_cells terciles:", [round(float(d.era5_max_fp.rank().corr(d.scat_max.rank())), 3) for _, d in g.groupby("cb")])
    ra, rb = g.era5_max_fp.rank(), g.scat_max.rank(); rn = g.n_cells.rank()
    res = lambda y: y - np.polyval(np.polyfit(rn, y, 1), rn)
    print(b, "partial rho given n_cells:", round(float(np.corrcoef(res(ra), res(rb))[0, 1]), 3))
