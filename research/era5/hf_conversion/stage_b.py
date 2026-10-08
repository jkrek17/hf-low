"""Test 1, stage B: 0.25 degree environment scalars and storm-relative composites, converters vs non-converters. ERA5 proxy, pipeline A. See PREREGISTRATION.md.
usage: stage_b.py REPO_ROOT ERA5_WORK [B] [B_PIXEL]. Reads $WORK/hf_conversion_pull/{stats,boxes}. Writes results/stage_b_tests.csv, stage_b_coverage.csv, comp_<basin>_L<lag>.npz, stage_b_summary.txt"""
import sys, os, glob, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_vs_storm"))
from strat import cells, boot, bh
root, W = sys.argv[1], sys.argv[2]; B = int(sys.argv[3]) if len(sys.argv) > 3 else 2000; BP = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
OUT = os.path.join(root, "research/era5/hf_conversion/results"); rng = np.random.default_rng(20261012); SEASONS = np.arange(2004, 2026); LAGS = [0, 12, 24]
S = pd.read_csv(os.path.join(OUT, "conv_sample.csv"))
R = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(W, "hf_conversion_pull", "stats", "*.csv")))])
a = R.anchor.str.extract(r"g(\d)_p(\d+)_L(\d+)").astype(int); R["grp"], R["pair"], R["lag"] = a[0], a[1], a[2]
R["ok"] = R.snapped & ((R.lag != 0) | (R.match_km <= 25))
SC = ["jet_dist", "jet_cos", "thetae850", "stab", "baroc"]
lines = [f"lag rows pulled {len(R)}; snapped within 150 km {int(R.snapped.sum())}; usable ({'lag 0 within 25 km'}) {int(R.ok.sum())}"]
cov = R.groupby(["basin", "grp", "lag"]).agg(n=("ok", "size"), usable=("ok", "sum")).reset_index(); cov.to_csv(os.path.join(OUT, "stage_b_coverage.csv"), index=False)
M = S.rename(columns={"grp": "g"})[["track", "pair", "g", "season", "mon"]]
R = R.merge(M, left_on=["track", "pair", "grp"], right_on=["track", "pair", "g"], how="left")
assert R.season.notna().all() and R[SC].notna().any().all()
res = []
for b in ("atl", "pac"):
    for L in LAGS:
        x = R[(R.basin == b) & (R.lag == L) & R.ok]
        sm, cn = cells(x[SC].values.astype(float), x.grp.values, x.season.values.astype(int), x.mon.values.astype(int) - 1, SEASONS, 12)
        obs, reps, p = boot(sm, cn, B, rng); se = reps.std(0); sd = np.sqrt((x[x.grp == 0][SC].var().values + x[x.grp == 1][SC].var().values) / 2)
        for j, v in enumerate(SC):
            res.append(dict(basin=b, var=v, lag=-L, conv=obs[1][j], non=obs[2][j], diff=obs[0][j], lo=np.percentile(reps[:, j], 2.5), hi=np.percentile(reps[:, j], 97.5), se=se[j], p=p[j],
                            mde80=2.8 * se[j], sd=sd[j], d_std=obs[0][j] / sd[j], mde_std=2.8 * se[j] / sd[j], n_conv=int(x[x.grp == 0][v].notna().sum()), n_non=int(x[x.grp == 1][v].notna().sum())))
T = pd.DataFrame(res); T["q"] = bh(T.p.values)[0]; T.to_csv(os.path.join(OUT, "stage_b_tests.csv"), index=False)
lines.append(f"primary family: {int((T.q < 0.05).sum())} of {len(T)} pass BH q<0.05")
boxes = {}
for f in sorted(glob.glob(os.path.join(W, "hf_conversion_pull", "boxes", "*.npz"))):
    z = np.load(f)
    for k in z.files: boxes[k] = z[k]
FIELDS = ["msl", "ws250", "te", "stab", "d2m", "gt"]
for b in ("atl", "pac"):
    for L in (0, 24):
        x = R[(R.basin == b) & (R.lag == L) & R.ok].reset_index(drop=True)
        x = x[[f"{r.track}_{r.anchor}_msl" in boxes for r in x.itertuples()]].reset_index(drop=True)
        out = {"n_conv": int((x.grp == 0).sum()), "n_non": int((x.grp == 1).sum())}
        for fld in FIELDS:
            arr = np.stack([boxes[f"{r.track}_{r.anchor}_{fld}"].astype(np.float32).ravel() for r in x.itertuples()])
            sm, cn = cells(arr, x.grp.values, x.season.values.astype(int), x.mon.values.astype(int) - 1, SEASONS, 12)
            obs, _, p = boot(sm, cn, BP, rng, keep=False); qq, rej = bh(p)
            for nm, v in (("diff", obs[0]), ("conv", obs[1]), ("non", obs[2])): out[f"{fld}_{nm}"] = v.reshape(121, 121).astype(np.float32)
            out[f"{fld}_rej"] = rej.reshape(121, 121); out[f"{fld}_p"] = p.reshape(121, 121).astype(np.float32)
            lines.append(f"{b} lag -{L} h {fld}: {int(rej.sum())} of {rej.size} pixels pass BH q<0.05")
        np.savez_compressed(os.path.join(OUT, f"comp_{b}_L{L}.npz"), **out)
open(os.path.join(OUT, "stage_b_summary.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
pd.set_option("display.width", 250); print(T[["basin", "var", "lag", "conv", "non", "diff", "d_std", "mde_std", "n_conv", "n_non", "q"]].round(2).to_string())
