"""RA-20 time selection and counts (labels only; no ERA5 field is read here).

Writes results/times_armB.csv (the times to pull) and results/selection.txt.
Run: python3 -I research/era5/explosive_onset/select_times.py <repo_root>
"""
import sys, os, numpy as np, pandas as pd
root = sys.argv[1]
f = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"))
SEED, N_NONCASE = 20261008, 600
E = f[(f.g800 < 55) & (~f.hf_now)].copy()      # eligible: sub-55 kt gust index, not HF now (PR 12/68/76 convention: all fixes, left24 kept)
E["case"] = E.hf24.astype(bool)
tc = E.groupby("time").agg(n_elig=("case", "size"), n_case=("case", "sum")).reset_index()
case_times = tc[tc.n_case > 0].time.to_numpy()
non = tc[tc.n_case == 0].time.to_numpy()
rng = np.random.RandomState(SEED)
pick = np.sort(rng.choice(non, N_NONCASE, replace=False))
out = pd.concat([tc[tc.time.isin(case_times)].assign(kind="case_time"),
                 tc[tc.time.isin(pick)].assign(kind="noncase_sample", weight=len(non) / N_NONCASE)]).sort_values("time")
out["weight"] = out.weight.fillna(1.0)
out.to_csv(os.path.join(root, "research/era5/explosive_onset/results/times_armB.csv"), index=False)
c = E[E.case]
L = []
L.append(f"fixes in table {len(f)}; eligible (g800<55, not HF now) {len(E)}; hf24 positives among all fixes {int(f.hf24.sum())}")
L.append(f"cases whose track leaves the domain within 24 h (kept) {int(E[E.hf24].left24.sum())}")
L.append(f"CASES {len(c)} (Atlantic {int((c.basin=='atl').sum())}, Pacific {int((c.basin=='pac').sum())}), tracks {c.track.nunique()}, seasons {c.season.nunique()}")
L.append(f"cases per season: min {c.groupby('season').size().min()}, median {int(c.groupby('season').size().median())}, max {c.groupby('season').size().max()}")
L.append(f"case times {len(case_times)}; non-case times available {len(non)}; sampled {N_NONCASE} (seed {SEED}), weight {len(non)/N_NONCASE:.2f} each")
ct = E[E.time.isin(case_times)]
L.append(f"eligible fixes at case times {len(ct)}; non-case (controls) {int((~ct.case).sum())}; controls per case {(~ct.case).sum()/len(c):.1f}")
dctl = ct[(~ct.case) & (ct.dp12 <= -3.6)]
L.append(f"deepening (dp12<=-3.6) controls at case times {len(dctl)}; cases that are deepening {int((c.dp12<=-3.6).sum())}; cases with dp12 missing {int(c.dp12.isna().sum())}")
# matched 1:k feasibility: same time, |g800 diff|<=3 kt, |dp12 diff|<=2 hPa
nm = []
ctl = ct[~ct.case]
g = {t: d for t, d in ctl.groupby("time")}
for _, r in c.iterrows():
    d = g.get(r.time)
    if d is None: nm.append(0); continue
    ok = (abs(d.g800 - r.g800) <= 3) & (((abs(d.dp12 - r.dp12) <= 2)) if not pd.isna(r.dp12) else d.dp12.isna())
    nm.append(int(ok.sum()))
nm = np.array(nm)
L.append(f"cases with >=1 matched control (same time, gust within 3 kt, dp12 within 2 hPa, or both dp12 missing): {int((nm>0).sum())} of {len(c)}; median matches {int(np.median(nm))}")
L.append(f"stratum base rate (eligible, all fixes) {E.case.mean():.4f}")
open(os.path.join(root, "research/era5/explosive_onset/results/selection.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
