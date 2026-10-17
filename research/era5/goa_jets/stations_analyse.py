"""Coastal-station test, H1 and H2 as registered in PREREGISTRATION_stations.md. Reads results/stations_pairs.csv, hrrr_features.csv, stage0_zone_gust.csv.
Writes results/stations_h1.csv, stations_h2.csv, stations_analysis.txt. usage: stations_analyse.py"""
import os, math
import numpy as np, pandas as pd
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, "results")
ZS = {"Z1": ["pila2", "mdxa2"], "Z2": ["cdea2", "lcna2", "gixa2"], "Z3": ["pila2", "amaa2", "fila2", "auga2"]}
NAME = {"pila2": "Pilot Rock", "mdxa2": "Middleton Island", "cdea2": "Cape Decision", "lcna2": "Lincoln Island", "gixa2": "Guard Island", "amaa2": "East Amatuli", "fila2": "Flat Island", "auga2": "Augustine Island", "cspa2": "Cape Spencer"}
rng = np.random.default_rng(16); out = []; P = lambda *a: out.append(" ".join(str(x) for x in a))
season = lambda t: t.year if t.month >= 6 else t.year - 1
p = pd.read_csv(os.path.join(R, "stations_pairs.csv"), parse_dates=["time"]); p["season"] = [season(t) for t in p.time]
# ---- H1
P("H1: station gust >= 35 kt; median model/station gust ratio and Spearman, nearest cell (registered) and 3x3 max (beside it)")
rows = []
for s, d in p.groupby("station"):
    x = d[d.obs_kt >= 35]
    if len(x) < 2: continue
    r = dict(station=s, name=NAME[s], n_obs=len(d), n35=len(x), hrrr_ratio=(x.hrrr_nn_kt / x.obs_kt).median(), hrrr3_ratio=(x.hrrr_max3_kt / x.obs_kt).median(), era5_ratio=(x.era5_kt / x.obs_kt).median(),
             hrrr_rho=spearmanr(x.hrrr_nn_kt, x.obs_kt)[0], era5_rho=spearmanr(x.era5_kt, x.obs_kt)[0])
    r["hrrr_closer"] = abs(math.log(r["hrrr_ratio"])) < abs(math.log(r["era5_ratio"])); rows.append(r)
h1 = pd.DataFrame(rows).round(3); h1.to_csv(os.path.join(R, "stations_h1.csv"), index=False); P(h1.to_string(index=False))
el = h1[h1.n35 >= 30]; P(f"stations with >= 30 such times: {len(el)}; HRRR median ratio closer to 1 than ERA5 at {int(el.hrrr_closer.sum())} of them (rule: 5 or more -> 'HRRR is the better coastal reference')")
# ---- H2
hf = pd.read_csv(os.path.join(R, "hrrr_features.csv"), parse_dates=["time"]); s0 = pd.read_csv(os.path.join(R, "stage0_zone_gust.csv"), parse_dates=["time"])
rows = []
for z, sts in ZS.items():
    e = s0[["time", f"{z}_gust_kt"]].rename(columns={f"{z}_gust_kt": "era5"}); h = hf[hf.zone == z][["time", "hrrr_gust"]]
    d = h.merge(e, on="time"); d["season"] = [season(t) for t in d.time]
    q = p[p.station.isin(sts)].groupby("time").obs_kt.agg(["max", "size"]).reset_index().rename(columns={"max": "obs_max", "size": "n_st"})
    d = d.merge(q, on="time"); d["hit"] = d.obs_max >= 45
    d["cls"] = np.where((d.hrrr_gust >= 54) & (d.era5 < 54), "HRRR-only", np.where((d.hrrr_gust >= 54) & (d.era5 >= 54), "both", np.where((d.hrrr_gust < 54) & (d.era5 < 54), "neither", "ERA5-only")))
    rate = lambda x: 100 * x.hit.mean() if len(x) else np.nan
    ss = sorted(d.season.unique())
    def ratio(dd):
        a, b = dd[dd.cls == "HRRR-only"], dd[dd.cls == "neither"]
        return rate(a) / rate(b) if len(a) and len(b) and rate(b) > 0 else np.nan
    bs = []
    for _ in range(2000):
        dd = pd.concat([d[d.season == k] for k in rng.choice(ss, len(ss))]); bs.append(ratio(dd))
    bs = np.array(bs); lo, hi = (np.nanpercentile(bs, [2.5, 97.5]) if np.isfinite(bs).sum() > 100 else (np.nan, np.nan))
    n1 = int((d.cls == "HRRR-only").sum()); r0 = rate(d[d.cls == "neither"]) / 100
    # minimum detectable rate ratio (two-proportion z, alpha .05 two-sided, power .8) with n1 vs the neither group
    n0 = int((d.cls == "neither").sum()); mde = np.nan
    for rr in np.arange(1.05, 60, 0.05):
        p1 = min(r0 * rr, 0.999); pb = (n1 * p1 + n0 * r0) / (n1 + n0); se0 = math.sqrt(pb * (1 - pb) * (1 / max(n1, 1) + 1 / n0)); se1 = math.sqrt(p1 * (1 - p1) / max(n1, 1) + r0 * (1 - r0) / n0)
        if se1 > 0 and (p1 - r0 - 1.96 * se0) / se1 >= 0.84: mde = round(float(rr), 2); break
    rt = ratio(d)
    verdict = "cannot tell" if n1 < 30 else ("supported" if rt >= 3 else "partly supported" if rt >= 1.5 else "not supported")
    for c in ["HRRR-only", "both", "neither", "ERA5-only"]:
        x = d[d.cls == c]; rows.append(dict(zone=z, cls=c, n=len(x), hit_per_100=round(rate(x), 1) if len(x) else np.nan))
    rows.append(dict(zone=z, cls="RATIO HRRR-only/neither", n=n1, hit_per_100=round(rt, 2), ci_lo=round(lo, 2), ci_hi=round(hi, 2), mde_ratio=mde, verdict=verdict))
h2 = pd.DataFrame(rows); h2.to_csv(os.path.join(R, "stations_h2.csv"), index=False)
P("\nH2: share of times with a zone-station gust >= 45 kt, by zone class (times with at least one zone station reporting); ratio HRRR-only / neither with 8-season bootstrap interval, minimum detectable ratio, registered verdict")
P(h2.to_string(index=False))
open(os.path.join(R, "stations_analysis.txt"), "w").write("\n".join(out) + "\n"); print("\n".join(out))
