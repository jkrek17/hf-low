"""Test 3, stage A: HF vs storm-force-only life cycle from the 00/12 UTC track table (no pull). See PREREGISTRATION.md.
ERA5 proxy, pipeline A. usage: stage_a.py REPO_ROOT [B]"""
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_vs_storm"))
from strat import cells, boot, bh
root = sys.argv[1]; B = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
OUT = os.path.join(root, "research/era5/hf_lifecycle/results"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261012)
SEASONS = np.arange(2004, 2026); LAGS = [-48, -36, -24, -12, 0, 12, 24]
S = pd.read_csv(os.path.join(root, "research/era5/hf_vs_storm/results/storms.csv"))
F = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"), usecols=["track", "time", "lat", "lon", "msl", "dp12", "g800", "speed", "age"])
F["tm"] = pd.to_datetime(F.time.astype(str), format="%Y%m%d%H")
F = F[F.track.isin(S.track)]
nfix = F.groupby("track").size()
keep = nfix[nfix >= 2].index
lines = [f"storms {len(S)} (HF {int((S.grp=='HF').sum())}, SF {int((S.grp=='SF').sum())}); with 00/12 fixes {F.track.nunique()}; with >=2 fixes {len(keep)}"]
F = F[F.track.isin(keep)].sort_values(["track", "tm"])
anc = F.loc[F.groupby("track").g800.idxmax()].set_index("track")   # ties: first (earlier) fix
VARS = ["msl", "dp12", "g800", "speed", "lat", "age"]
Fi = F.set_index(["track", "tm"])
rows = []
tracks = anc.index.values
val = {v: np.full((len(tracks), len(LAGS)), np.nan) for v in VARS}
for i, tr in enumerate(tracks):
    t0 = anc.loc[tr, "tm"]
    g = F[F.track == tr].set_index("tm")
    for j, L in enumerate(LAGS):
        t = t0 + pd.Timedelta(hours=L)
        if t in g.index:
            for v in VARS:
                val[v][i, j] = g.loc[t, v]
D = S.set_index("track").loc[tracks].reset_index()
grp = (D.grp == "SF").astype(int).values; season = D.season.values; stratum = D.mon.values - 1
res = []; first = {}
for sel_name, selmask in (("all", np.ones(len(D), bool)), ("complete", ~np.isnan(val["msl"]).any(1) & ~np.isnan(val["dp12"]).any(1))):
    for b in ("atl", "pac"):
        m = (D.basin.values == b) & selmask
        for v in VARS:
            x = val[v][m]
            sm, cn = cells(x, grp[m], season[m], stratum[m], SEASONS, 12)
            obs, reps, p = boot(sm, cn, B, rng)
            se = reps.std(0)
            nh = np.array([(~np.isnan(x[grp[m] == 0, j])).sum() for j in range(len(LAGS))]); ns = np.array([(~np.isnan(x[grp[m] == 1, j])).sum() for j in range(len(LAGS))])
            sd = np.sqrt((np.nanvar(x[grp[m] == 0], 0) + np.nanvar(x[grp[m] == 1], 0)) / 2)
            for j, L in enumerate(LAGS):
                res.append(dict(sel=sel_name, basin=b, var=v, lag=L, hf=obs[1][j], sf=obs[2][j], diff=obs[0][j], lo=np.percentile(reps[:, j], 2.5),
                                hi=np.percentile(reps[:, j], 97.5), se=se[j], p=p[j], mde80=2.8 * se[j], sd=sd[j], d_std=obs[0][j] / sd[j], mde_std=2.8 * se[j] / sd[j], n_hf=nh[j], n_sf=ns[j]))
            first[(sel_name, b, v)] = (obs[0], reps, se)
R = pd.DataFrame(res)
fam = R[(R["sel"] == "all") & R["var"].isin(["msl", "dp12"])].copy()
fam["q"] = bh(fam.p.values)[0]
R = R.merge(fam[["basin", "var", "lag", "q"]], on=["basin", "var", "lag"], how="left").assign()
R.loc[R["sel"] != "all", "q"] = np.nan
R.to_csv(os.path.join(OUT, "stage_a_tests.csv"), index=False)
lines.append(f"primary family: {int((fam.q < 0.05).sum())} of {len(fam)} pass BH q<0.05")
# first persistent lag: earliest lag from which q<0.05 holds at every later lag through 0 h (lag index <= 4)
def first_lag(sig):
    idx = [j for j in range(5) if all(sig[k] for k in range(j, 5))]
    return LAGS[idx[0]] if idx else None
out = {}
for b in ("atl", "pac"):
    for v in ("msl", "dp12"):
        sg = (fam[(fam.basin == b) & (fam["var"] == v)].sort_values("lag").q.values < 0.05)
        pt = first_lag(sg)
        d, reps, se = first[("all", b, v)]
        z = np.abs(reps) / se
        sgn = np.sign(reps) == np.sign(d)
        fl = []
        for r in range(len(reps)):
            s_ = (z[r] > 1.96) & sgn[r]
            fl.append(first_lag(s_))
        vals = [x if x is not None else 99 for x in fl]
        dist = {str(k): float(np.mean(np.array(vals) == k)) for k in [-48, -36, -24, -12, 0, 99]}
        lo, hi = np.percentile(vals, [2.5, 97.5])
        out[f"{b}_{v}"] = dict(first_persistent_lag=pt, resample_dist=dist, resample_p2_5=float(lo), resample_p97_5=float(hi))
json.dump(out, open(os.path.join(OUT, "stage_a_first_lag.json"), "w"), indent=1)
lines.append(json.dumps(out))
open(os.path.join(OUT, "stage_a_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
pd.set_option("display.width", 250)
print(R[(R["sel"] == "all") & R["var"].isin(["msl", "dp12"])][["basin", "var", "lag", "hf", "sf", "diff", "lo", "hi", "d_std", "mde_std", "n_hf", "n_sf", "q"]].round(2).to_string())
