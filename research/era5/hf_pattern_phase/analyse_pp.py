"""Test 2 analysis: structure of HF lows under top vs bottom tercile of the hemispheric-pattern index. ERA5 proxy, pipeline A.
usage: analyse_pp.py ERA5_WORK [B]. Reads stats from $WORK/hf_vs_storm (800 onsets+peaks) and $WORK/hf_pattern_phase (new onsets). Writes results/pp_*.csv, pp_summary.txt"""
import sys, os, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_vs_storm"))
from strat import cells, boot, bh
W = sys.argv[1]; B = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "results")
rng = np.random.default_rng(20261012); SEASONS = np.arange(2004, 2026)
st = pd.concat([pd.read_csv(f) for d in ("hf_vs_storm", "hf_pattern_phase") for f in glob.glob(os.path.join(W, d, "stats", "*.csv"))])
H = pd.read_csv(os.path.join(OUT, "hf_index.csv")); H = H[H.in_scope].copy()
lines = [f"in-scope HF storms {len(H)} (Atl {int((H.basin=='atl').sum())}, Pac {int((H.basin=='pac').sum())}); out of scope {int((~pd.read_csv(os.path.join(OUT,'hf_index.csv')).in_scope).sum())}"]
SC = ["pc", "msl_grad", "gmax", "wsmax", "a_g48", "gmax_r", "g48_rmax", "g48_right", "d2m_500", "gust_factor", "a_hf"]
def build(anchor):
    a = st[st.anchor == anchor].drop_duplicates(["track", "anchor"])
    return H.merge(a[["track"] + SC + ["lat", "match_km", "g800", "g800_cat", "gmax_near_grn_ice"]], on="track", how="inner")
for anchor in ("HF_onset", "HF_peak"):
    d = build(anchor)
    lines.append(f"{anchor}: joined {len(d)}; match_km>25 {int((d.match_km>25).sum())}; g800 off by >0.5 kt {int(((d.g800-d.g800_cat).abs()>0.5).sum())} of {int(d.g800_cat.notna().sum())}")
O = build("HF_onset"); O = O[(O.match_km <= 25)]
lines.append(f"onsets analysed {len(O)}; tercile counts " + ", ".join(f"{b}/{t} {int(((O.basin==b)&(O.tercile==t)).sum())}" for b in ("atl","pac") for t in ("bottom","middle","top")))
assert O[SC].notna().any().all(), "all-NaN scalar"
def tern(D, strat_fn, nstr, label, rows=None):
    res = []
    for b in ("atl", "pac"):
        m = (D.basin == b) & D.tercile.isin(["top", "bottom"]); x = D[m]
        grp = (x.tercile == "bottom").astype(int).values            # 0 = top, 1 = bottom
        vals = x[SC].values.astype(float)
        sm, cn = cells(vals, grp, x.season.values, strat_fn(x), SEASONS, nstr)
        obs, reps, p = boot(sm, cn, B, rng)
        se = reps.std(0); sd = x[SC].std().values
        for j, v in enumerate(SC):
            res.append(dict(family=label, basin=b, var=v, top=obs[1][j], bottom=obs[2][j], diff=obs[0][j], lo=np.percentile(reps[:, j], 2.5), hi=np.percentile(reps[:, j], 97.5),
                            se=se[j], p=p[j], mde80=2.8 * se[j], sd=sd[j], d_std=obs[0][j] / sd[j], mde_std=2.8 * se[j] / sd[j],
                            n_top=int(((grp == 0) & x[v].notna().values).sum()), n_bot=int(((grp == 1) & x[v].notna().values).sum())))
    R = pd.DataFrame(res)
    if rows is not None: R = R[R["var"].isin(rows)]
    R["q"] = bh(R.p.values)[0]; return R
mon = lambda x: x.mon.values - 1
prim_vars = ["pc", "msl_grad", "gmax", "wsmax", "a_g48", "gmax_r", "g48_rmax", "g48_right", "d2m_500", "gust_factor", "a_hf"]
P = tern(O, lambda x: (x.mon.values - 1), 12, "primary"); P.to_csv(os.path.join(OUT, "pp_primary.csv"), index=False)
S1 = tern(O, lambda x: np.clip(((x.onset_lat.values - 20) // 5).astype(int), 0, 11), 12, "S1_latmatched"); S1.to_csv(os.path.join(OUT, "pp_S1.csv"), index=False)
O3 = O[~O.gmax_near_grn_ice.astype(bool)]
S3 = tern(O3, lambda x: x.mon.values - 1, 12, "S3_no_grn_ice"); S3 = S3[S3.basin == "atl"].copy(); S3["q"] = bh(S3.p.values)[0]; S3.to_csv(os.path.join(OUT, "pp_S3.csv"), index=False)
Pk = build("HF_peak"); Pk = Pk[Pk.match_km <= 25]
S4 = tern(Pk, lambda x: x.mon.values - 1, 12, "S4_peak_anchor"); S4.to_csv(os.path.join(OUT, "pp_S4.csv"), index=False)
# S2: slope of each scalar on the index z, month fixed effects, season-block bootstrap, all in-scope storms
def slope(D, idx_seasons):
    out = []
    for v in SC:
        x = D[["z", "mon", v]].dropna(); y = x[v] - x.groupby("mon")[v].transform("mean"); zz = x.z - x.groupby("mon").z.transform("mean")
        out.append((zz * y).sum() / (zz ** 2).sum())
    return np.array(out)
S2 = []
for b in ("atl", "pac"):
    D = O[(O.basin == b) & O.z.notna()]; obs = slope(D, None)
    by = {s: g for s, g in D.groupby("season")}; reps = []
    for _ in range(B):
        ch = rng.choice(SEASONS, len(SEASONS)); reps.append(slope(pd.concat([by[s] for s in ch if s in by]), None))
    reps = np.array(reps); se = reps.std(0); p = 2 * np.minimum((reps > 0).mean(0) + 1 / (B + 1), (reps < 0).mean(0) + 1 / (B + 1)).clip(max=1)
    for j, v in enumerate(SC):
        sd = D[v].std(); S2.append(dict(family="S2_slope", basin=b, var=v, slope_per_SD=obs[j], lo=np.percentile(reps[:, j], 2.5), hi=np.percentile(reps[:, j], 97.5), se=se[j], p=p[j], mde80=2.8 * se[j], sd=sd, d_std=obs[j] / sd, mde_std=2.8 * se[j] / sd, n=len(D)))
S2 = pd.DataFrame(S2); S2["q"] = bh(S2.p.values)[0]; S2.to_csv(os.path.join(OUT, "pp_S2.csv"), index=False)
for nm, R in (("primary", P), ("S1", S1), ("S2", S2), ("S3", S3), ("S4", S4)):
    lines.append(f"{nm}: {int((R.q < 0.05).sum())} of {len(R)} pass BH q<0.05")
open(os.path.join(OUT, "pp_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines)); pd.set_option("display.width", 250)
print(P[["basin", "var", "top", "bottom", "diff", "d_std", "mde_std", "n_top", "n_bot", "p", "q"]].round(3).to_string())
