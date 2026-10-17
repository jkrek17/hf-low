"""Numbers quoted in README_stations.md that the registered H1/H2 outputs do not hold (added after the claim review; ALL POST HOC unless said). Writes results/stations_posthoc.txt. usage: stations_posthoc.py"""
import os
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, "results"); out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
ZS = {"Z1": ["pila2", "mdxa2"], "Z2": ["cdea2", "lcna2", "gixa2"], "Z3": ["pila2", "amaa2", "fila2", "auga2"]}
p = pd.read_csv(os.path.join(R, "stations_pairs.csv"), parse_dates=["time"]); hf = pd.read_csv(os.path.join(R, "hrrr_features.csv"), parse_dates=["time"]); s0 = pd.read_csv(os.path.join(R, "stage0_zone_gust.csv"), parse_dates=["time"])
P(f"station sample: pairs {len(p)}, distinct times {p.time.nunique()}, first {p.time.min()}, last {p.time.max()}")
for s, d in p.groupby("station"): P(f"  {s}: {len(d)} pairs, {d.time.min():%Y-%m-%d} to {d.time.max():%Y-%m-%d}")
P("\nH1 beside the registered nearest-cell ratio: 3x3-max HRRR/station median (obs >= 35 kt)")
for s, d in p[p.obs_kt >= 35].groupby("station"): P(f"  {s}: n {len(d)}, 3x3 {(d.hrrr_max3_kt / d.obs_kt).median():.3f}, nearest {(d.hrrr_nn_kt / d.obs_kt).median():.3f}, ERA5 {(d.era5_kt / d.obs_kt).median():.3f}")
P("\nselection check: the same ratios when the TIMES are chosen by the model's own value >= 35 kt instead of the station's (median model/station)")
for s, d in p.groupby("station"):
    a, b = d[d.hrrr_nn_kt >= 35], d[d.era5_kt >= 35]
    P(f"  {s}: HRRR n {len(a)} {(a.hrrr_nn_kt / a.obs_kt).median():.2f}; ERA5 n {len(b)} {(b.era5_kt / b.obs_kt).median():.2f}")
P("\nH2 classes: zone gust medians (kt) and times with a reporting zone station")
for z, sts in ZS.items():
    e = s0[["time", f"{z}_gust_kt"]].rename(columns={f"{z}_gust_kt": "era5"}); h = hf[hf.zone == z][["time", "hrrr_gust", "g_dland_km"]]
    d = h.merge(e, on="time"); q = p[p.station.isin(sts)].groupby("time").obs_kt.max().rename("om").reset_index(); d = d.merge(q, on="time"); d["hit"] = d.om >= 45
    d["cls"] = np.where((d.hrrr_gust >= 54) & (d.era5 < 54), "HRRR-only", np.where((d.hrrr_gust >= 54) & (d.era5 >= 54), "both", np.where((d.hrrr_gust < 54) & (d.era5 < 54), "neither", "ERA5-only")))
    P(f" {z}: times {len(d)}")
    for c in ["HRRR-only", "both", "neither"]:
        x = d[d.cls == c]; P(f"   {c}: n {len(x)}, hit {100 * x.hit.mean():.1f}/100, median ERA5 zone gust {x.era5.median():.1f}, HRRR zone gust {x.hrrr_gust.median():.1f}, HRRR peak-cell distance to land {x.g_dland_km.median():.1f} km")
    r = 100 * d[d.cls == "HRRR-only"].hit.mean() / (100 * d[d.cls == "both"].hit.mean()); P(f"   HRRR-only rate / both rate: {r:.2f}")
    # matched-ERA5 version: HRRR-only vs neither, restricted to ERA5 zone gust 40-54 kt
    m = d[(d.era5 >= 40) & (d.era5 < 54)]; a, b = m[m.cls == "HRRR-only"], m[m.cls == "neither"]
    P(f"   POST HOC matched ERA5 zone gust 40-54 kt: HRRR-only n {len(a)} hit {100 * a.hit.mean():.1f}; neither n {len(b)} hit {100 * b.hit.mean():.1f}")
d = hf[hf.zone == "Z1"].merge(s0[["time", "Z1_gust_kt"]], on="time"); q = p[p.station == "pila2"].groupby("time").obs_kt.max().rename("om").reset_index(); d = d.merge(q, on="time")
a = d[(d.hrrr_gust >= 54) & (d.Z1_gust_kt < 54)]; b = d[(d.hrrr_gust < 54) & (d.Z1_gust_kt < 54)]
P(f"\nPilot Rock alone in Z1: HRRR-only n {len(a)} hit {100 * (a.om >= 45).mean():.1f}; neither n {len(b)} hit {100 * (b.om >= 45).mean():.1f}")
open(os.path.join(R, "stations_posthoc.txt"), "w").write("\n".join(out) + "\n"); print("\n".join(out))
