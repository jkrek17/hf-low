"""Test 1, stage A: converters vs non-converters among deepening storms, environment table at t0, t0-12 h, t0-24 h (no pull). See PREREGISTRATION.md.
ERA5 proxy, pipeline A. usage: stage_a.py REPO_ROOT [B]. Writes results/pairs.csv, stage_a_tests.csv, stage_a_change.csv, stage_a_latmatched.csv, stage_a_confounds.csv, stage_a_summary.txt"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_vs_storm"))
from strat import cells, boot, bh
root = sys.argv[1]; B = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
E5 = os.path.join(root, "research/era5"); OUT = os.path.join(E5, "hf_conversion/results"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261012); SEASONS = np.arange(2004, 2026)
F = pd.read_csv(os.path.join(E5, "intensity/results/fixes_2004.csv.gz"))
E = pd.read_csv(os.path.join(E5, "intensity/results/env_2004.csv.gz"))
F = F.merge(E, on=["track", "time"], how="left")
F["tm"] = pd.to_datetime(F.time.astype(str), format="%Y%m%d%H"); F["mon"] = F.tm.dt.month
ENV = ["jet250", "eady", "sstgrad", "sst_t500", "tcwv", "flux", "vadv500", "div300"]
el = F[(F.g800 < 71.7) & (F.ndr24 >= 1.0) & F.ndr24.notna()].copy()
el["conv"] = el.hf24.astype(bool)
lines = [f"eligible fixes {len(el)} on {el.track.nunique()} storms; converting {int(el.conv.sum())} fixes / {el[el.conv].track.nunique()} storms"]
pick = el.sample(frac=1, random_state=20261012).drop_duplicates(["track", "conv"])    # one random eligible fix per storm per group
both = set(pick[pick.conv].track) & set(pick[~pick.conv].track)
lines.append(f"one fix per storm per group: converters {int(pick.conv.sum())}, non-converters {int((~pick.conv).sum())}; storms in both groups {len(both)}")
pairs = []
for (b, m), g in pick.groupby(["basin", "mon"]):
    c = g[g.conv].sample(frac=1, random_state=int(m) + 100 * (b == "pac")); n = g[~g.conv].copy()
    for r in c.itertuples():
        if len(n) == 0: break
        d = (n.ndr24 - r.ndr24).abs()
        j = d.idxmin()
        if d[j] < 0.25:
            pairs.append((r.Index, j, b, m)); n = n.drop(j)
P = pd.DataFrame(pairs, columns=["ic", "in", "basin", "mon"])
lines.append(f"matched pairs {len(P)} (Atlantic {int((P.basin=='atl').sum())}, Pacific {int((P.basin=='pac').sum())}); converters unmatched {int(pick.conv.sum()) - len(P)}; mean |d ndr24| {np.abs(pick.loc[P.ic, 'ndr24'].values - pick.loc[P['in'], 'ndr24'].values).mean():.3f}")
sel = pd.concat([pick.loc[P.ic].assign(pair=range(len(P)), grp=0), pick.loc[P["in"]].assign(pair=range(len(P)), grp=1)]).reset_index(drop=True)
sel[["pair", "grp", "track", "time", "basin", "season", "mon", "lat", "lon", "msl", "ndr24", "age", "speed", "heading"]].to_csv(os.path.join(OUT, "pairs.csv"), index=False)
idx = F.set_index(["track", "time"])
LAGS = [0, 12, 24]
def lagged(var):
    out = np.full((len(sel), 3), np.nan)
    for j, L in enumerate(LAGS):
        t = (sel.tm - pd.Timedelta(hours=L)).dt.strftime("%Y%m%d%H").astype("int64")
        key = pd.MultiIndex.from_arrays([sel.track, t])
        out[:, j] = idx[var].reindex(key).values
    return out
VALS = {v: lagged(v) for v in ENV + ["lat", "age", "msl", "speed", "dp12"]}
grp = sel.grp.values; season = sel.season.values; stratum = sel.mon.values - 1
lines.append("missing share at t0: " + ", ".join(f"{v} {np.isnan(VALS[v][:, 0]).mean():.2f}" for v in ENV))
def run(vals, strat, nstr, label, lagnames):
    rows = []
    for b in ("atl", "pac"):
        m = (sel.basin.values == b)
        for v, x in vals.items():
            sm, cn = cells(x[m], grp[m], season[m], strat[m], SEASONS, nstr)
            obs, reps, p = boot(sm, cn, B, rng)
            se = reps.std(0); sd = np.sqrt((np.nanvar(x[m][grp[m] == 0], 0) + np.nanvar(x[m][grp[m] == 1], 0)) / 2)
            for j, ln in enumerate(lagnames):
                rows.append(dict(basin=b, var=v, lag=ln, conv=obs[1][j], non=obs[2][j], diff=obs[0][j], lo=np.percentile(reps[:, j], 2.5), hi=np.percentile(reps[:, j], 97.5),
                                 se=se[j], p=p[j], mde80=2.8 * se[j], sd=sd[j], d_std=obs[0][j] / sd[j], mde_std=2.8 * se[j] / sd[j],
                                 n_conv=int((~np.isnan(x[m][grp[m] == 0, j])).sum()), n_non=int((~np.isnan(x[m][grp[m] == 1, j])).sum())))
    R = pd.DataFrame(rows); R["q"] = bh(R.p.values)[0]; R["family"] = label
    return R
prim = run({v: VALS[v] for v in ENV}, stratum, 12, "primary", ["t0", "t0-12h", "t0-24h"])
prim.to_csv(os.path.join(OUT, "stage_a_tests.csv"), index=False)
chg = {v: (VALS[v][:, [0]] - VALS[v][:, [2]]) for v in ENV}
change = run(chg, stratum, 12, "change_24h", ["t0 minus t0-24h"]); change.to_csv(os.path.join(OUT, "stage_a_change.csv"), index=False)
latb = np.clip(((sel.lat.values - 20) // 5).astype(int), 0, 10)
latm = run({v: VALS[v] for v in ENV}, latb, 11, "latitude_matched", ["t0", "t0-12h", "t0-24h"]); latm.to_csv(os.path.join(OUT, "stage_a_latmatched.csv"), index=False)
conf = run({v: VALS[v][:, [0]] for v in ["lat", "age", "msl", "speed", "dp12"]}, stratum, 12, "confounds_descriptive", ["t0"]); conf.to_csv(os.path.join(OUT, "stage_a_confounds.csv"), index=False)
for nm, R in (("primary", prim), ("change", change), ("latitude-matched", latm)):
    lines.append(f"{nm}: {int((R.q < 0.05).sum())} of {len(R)} pass BH q<0.05")
open(os.path.join(OUT, "stage_a_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
pd.set_option("display.width", 250)
print(prim[["basin", "var", "lag", "conv", "non", "diff", "d_std", "mde_std", "n_conv", "n_non", "q"]].round(2).to_string())
print(conf[["basin", "var", "conv", "non", "diff", "d_std", "q"]].round(2).to_string())
