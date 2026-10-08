"""Test 3, stage B: structure scalars at 12 h steps around the peak, HF vs storm-force-only (250+250 per basin). ERA5 proxy, pipeline A. See PREREGISTRATION.md.
usage: stage_b.py REPO_ROOT ERA5_WORK [B]. Reads $WORK/hf_lifecycle/stats (6 new lags) and $WORK/hf_vs_storm/stats (lag 0). Writes results/stage_b_tests.csv, stage_b_first_lag.json, stage_b_coverage.csv, stage_b_summary.txt"""
import sys, os, glob, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "hf_vs_storm"))
from strat import cells, boot, bh
root, W = sys.argv[1], sys.argv[2]; B = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
OUT = os.path.join(root, "research/era5/hf_lifecycle/results"); rng = np.random.default_rng(20261012)
SEASONS = np.arange(2004, 2026); LAGS = [-48, -36, -24, -12, 0, 12, 24]
S = pd.read_csv(os.path.join(OUT, "lc_sample.csv"))
new = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(W, "hf_lifecycle", "stats", "*.csv"))])
old = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(W, "hf_vs_storm", "stats", "*.csv"))])
old = old[old.anchor.isin(["HF_peak", "SF_peak"])].merge(S[["track", "anchor"]], on=["track", "anchor"])
old["anchor"] = "lag+0"; old["snapped"] = old.match_km <= 25
new["lag"] = new.anchor.str.replace("lag", "").astype(int); old["lag"] = 0
A = pd.concat([new, old]); A["lag"] = A["lag"].astype(int)
MEAS = ["pc", "msl_grad", "g48_rmax", "gmax_r"]
chk = A[A.snapped & (A.lag != 0)]
lines = [f"sample {len(S)}; stats rows new {len(new)} (snapped {int(new.snapped.sum())}), lag 0 {len(old)} (match<=25 km {int(old.snapped.sum())})"]
# g800 reproduction at lags whose time is a 00/12 UTC fix (table g800 available)
F = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"), usecols=["track", "time", "g800"])
c = chk.assign(time=chk.time.astype("int64")).merge(F, on=["track", "time"], suffixes=("", "_tab")).dropna(subset=["g800_tab", "g800"])
nz0 = int((c.g800_tab == 0).sum()); c = c[c.g800_tab > 0]      # table g800 == 0 is "not computed" (30% of all fixes), not a measured zero
lines.append(f"g800 reproduction at 00/12 UTC lags: {len(c)} fixes with a table value > 0 ({nz0} more have table 0 = not computed, left out); within 0.5 kt {((c.g800 - c.g800_tab).abs() <= 0.5).mean():.4f}; within 1 kt {((c.g800 - c.g800_tab).abs() <= 1).mean():.4f}; max {(c.g800 - c.g800_tab).abs().max():.1f} kt")
assert A[MEAS].notna().any().all()
D = S.copy(); idx = {(t, a): i for i, (t, a) in enumerate(zip(D.track, D.anchor))}
val = {m: np.full((len(D), 7), np.nan) for m in MEAS}; snapd = np.zeros((len(D), 7), bool); expected = np.zeros((len(D), 7), bool)
T6 = pd.read_csv(os.path.join(OUT, "times_lc.csv"))
key = dict(zip(zip(S.track, S.basin), range(len(S))))
trk2row = {t: i for i, t in enumerate(D.track)}
for r in A.itertuples():
    i = trk2row[r.track]; j = LAGS.index(r.lag)
    expected[i, j] = True
    if r.snapped:
        snapd[i, j] = True
        for m in MEAS: val[m][i, j] = getattr(r, m)
for r in T6.itertuples():
    expected[trk2row[r.track], LAGS.index(int(r.anchor.replace("lag", "")))] = True
expected[:, 4] = True
cov = pd.DataFrame([dict(basin=b, grp="HF" if g == 0 else "SF", lag=L, n_storms=int(((D.basin == b) & (D.grp == g)).sum()), n_track_lag=int((expected[:, j] & (D.basin == b).values & (D.grp == g).values).sum()),
                         n_snapped=int((snapd[:, j] & (D.basin == b).values & (D.grp == g).values).sum()))
                    for b in ("atl", "pac") for g in (0, 1) for j, L in enumerate(LAGS)])
cov.to_csv(os.path.join(OUT, "stage_b_coverage.csv"), index=False)
grp = D.grp.values; season = D.season.values; stratum = D.mon.values - 1
res = []; keepfl = {}
for b in ("atl", "pac"):
    m = (D.basin.values == b)
    for v in MEAS:
        x = val[v][m]
        sm, cn = cells(x, grp[m], season[m], stratum[m], SEASONS, 12)
        obs, reps, p = boot(sm, cn, B, rng); se = reps.std(0)
        sd = np.sqrt((np.nanvar(x[grp[m] == 0], 0) + np.nanvar(x[grp[m] == 1], 0)) / 2)
        for j, L in enumerate(LAGS):
            res.append(dict(basin=b, var=v, lag=L, hf=obs[1][j], sf=obs[2][j], diff=obs[0][j], lo=np.percentile(reps[:, j], 2.5), hi=np.percentile(reps[:, j], 97.5), se=se[j], p=p[j],
                            mde80=2.8 * se[j], sd=sd[j], d_std=obs[0][j] / sd[j], mde_std=2.8 * se[j] / sd[j],
                            n_hf=int((~np.isnan(x[grp[m] == 0, j])).sum()), n_sf=int((~np.isnan(x[grp[m] == 1, j])).sum())))
        keepfl[(b, v)] = (obs[0], reps, se)
R = pd.DataFrame(res); R["q"] = bh(R.p.values)[0]; R.to_csv(os.path.join(OUT, "stage_b_tests.csv"), index=False)
lines.append(f"primary family: {int((R.q < 0.05).sum())} of {len(R)} pass BH q<0.05")
def first_lag(sig):
    ix = [j for j in range(5) if all(sig[k] for k in range(j, 5))]
    return LAGS[ix[0]] if ix else None
out = {}
for b in ("atl", "pac"):
    for v in MEAS:
        sg = R[(R.basin == b) & (R["var"] == v)].sort_values("lag").q.values < 0.05
        d, reps, se = keepfl[(b, v)]
        z = np.abs(reps) / se; sgn = np.sign(reps) == np.sign(d)
        fl = [first_lag((z[r] > 1.96) & sgn[r]) for r in range(len(reps))]
        vals = np.array([x if x is not None else 99 for x in fl])
        out[f"{b}_{v}"] = dict(first_persistent_lag=first_lag(sg), resample_dist={str(k): float(np.mean(vals == k)) for k in [-48, -36, -24, -12, 0, 99]}, resample_p2_5=float(np.percentile(vals, 2.5)), resample_p97_5=float(np.percentile(vals, 97.5)))
json.dump(out, open(os.path.join(OUT, "stage_b_first_lag.json"), "w"), indent=1)
lines.append(json.dumps({k: v["first_persistent_lag"] for k, v in out.items()}))
open(os.path.join(OUT, "stage_b_summary.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
pd.set_option("display.width", 250)
print(R[R["var"].isin(["pc", "msl_grad"])][["basin", "var", "lag", "hf", "sf", "diff", "d_std", "mde_std", "n_hf", "n_sf", "q"]].round(2).to_string())
