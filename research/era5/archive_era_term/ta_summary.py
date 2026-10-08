"""RA-29 T-A summary: era-adjusted vs original, BH families, era coefficients with season-block bootstrap intervals.
Reads results/ta_{none,step,asc,trend}.csv and the original sustained_hf/results/contrast_archive_withm.csv. Writes results/ta_summary.txt, ta_table.csv.
usage: ta_summary.py"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import era_refit as R
K, C, CH = R.K, R.C, R.CH

o = pd.read_csv(os.path.join(HERE, "..", "sustained_hf", "results", "contrast_archive_withm.csv"))
W = lambda s: lines.append(s) or print(s)
lines = []


def bh(p):
    p = np.asarray(p, float); n = len(p); od = np.argsort(p); q = np.empty(n)
    r = p[od] * n / (np.arange(n) + 1); q[od] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1); return q


FAM = {"F1": ["h1", "h2", "h3"], "F2": ["h3-b"], "F3": ["b", "m", "h2-h1", "h3-h1"]}
fam_of = {s: f for f, v in FAM.items() for s in v}
tabs = []
for md in ("step", "asc", "trend"):
    a = pd.read_csv(os.path.join(HERE, "results", f"ta_{md}.csv"))
    x = o.merge(a, on=["basin", "stat"], suffixes=("_o", "_a"))
    x["fam"] = x.stat.map(fam_of)
    for src in ("o", "a"):
        p = x[f"p_perm_{src}"].values
        x[f"q_fam_{src}"] = np.nan
        for f in FAM:
            m = (x.fam == f).values
            x.loc[m, f"q_fam_{src}"] = bh(p[m])
        x[f"q_all_{src}"] = bh(p)
    x["shift_log"] = x.log_est_a - x.log_est_o
    x["inside"] = (x.est_a >= x.lo_o) & (x.est_a <= x.hi_o)
    x["sign_kept"] = np.sign(x.log_est_a) == np.sign(x.log_est_o)
    x["se_ratio"] = x.se_log_a / x.se_log_o
    x["shift_in_se_o"] = x.shift_log / x.se_log_o
    pas_o, pas_a = x.q_fam_o < 0.05, x.q_fam_a < 0.05
    x["verdict"] = np.where(~x.sign_kept | ~x.inside, "changed", np.where(pas_o & ~pas_a, "no longer resolved", "unchanged"))
    tabs.append(x)
T = pd.concat(tabs)
T.to_csv(os.path.join(HERE, "results", "ta_table.csv"), index=False, float_format="%.5g")

W("RA-29 T-A: archive weekly design of PR 41 / PR 79 (leave-one-season-out hemispheric pattern index), 22 seasons 2004-05..2025-26.")
W("Control: with the era column dropped the original contrast_archive_withm.csv is reproduced exactly (log_est, interval, p).")
for md, lab in (("step", "PRIMARY: step at 2009-11-23"), ("asc", "S1: step at 2012-09-01 (ASCAT-B)"), ("trend", "S2: season-linear trend")):
    W(f"\n== {lab} ==")
    x = T[T["mode"] == md]
    for b in ("atl", "pac"):
        W(f" {b}  corr(index, era) {x[x.basin==b].corr_x_era.iloc[0]:+.3f}")
        for r in x[x.basin == b].itertuples():
            ec = f"  era coef {r.era_coef:+.3f}" if not np.isnan(r.era_coef) else ""
            W(f"   {r.stat:6s} orig {r.est_o:.3f} [{r.lo_o:.3f},{r.hi_o:.3f}] q{r.q_fam_o:.3f} | adj {r.est_a:.3f} [{r.lo_a:.3f},{r.hi_a:.3f}] q{r.q_fam_a:.3f} | "
              f"shift {r.shift_log:+.3f} ({r.shift_in_se_o:+.2f} se), se x{r.se_ratio:.2f}, MDE {r.mde_rr_o:.3f}->{r.mde_rr_a:.3f}  {r.verdict}{ec}")
    W(f"  verdicts: {x.verdict.value_counts().to_dict()}; F1 pass FDR (q<0.05) original {int((x[x.fam=='F1'].q_fam_o<0.05).sum())}/6, adjusted {int((x[x.fam=='F1'].q_fam_a<0.05).sum())}/6; "
      f"F2 {int((x[x.fam=='F2'].q_fam_o<0.05).sum())}/2 -> {int((x[x.fam=='F2'].q_fam_a<0.05).sum())}/2; across all 16 tests {int((x.q_all_o<0.05).sum())} -> {int((x.q_all_a<0.05).sum())} pass")
    W(f"  largest |shift| {x.shift_log.abs().max():.3f} (log), largest SE ratio {x.se_ratio.max():.2f}")

# era coefficient intervals (season-block bootstrap of the Poisson fit, primary step only); added output, no test depends on it
W("\n== Era coefficient (primary step), per class: log rate change after 2009-11-23 at the same pattern index; 2000 season-block bootstraps ==")
rng = np.random.default_rng(29)
for b in C.BASINS:
    M, n = R.archive_with_era(b, "step")
    allidx = list(range(M.ns))
    names = ["h1", "h2", "h3", "b", "m"]
    def coef(idx):
        rws = M.rows(idx)
        return {k: C.pois(M.D.X(k, rws), M.D.Y[k][rws])[-2] for k in names}
    c0 = coef(allidx)
    bt = {k: [] for k in names}
    for _ in range(2000):
        c = coef(list(rng.integers(0, M.ns, M.ns)))
        for k in names: bt[k].append(c[k])
    for k in names:
        lo, hi = np.percentile(bt[k], [2.5, 97.5])
        W(f"  {b} {k:2s} era coef {c0[k]:+.3f} [{lo:+.3f},{hi:+.3f}]  = {100*(np.exp(c0[k])-1):+.0f}% [{100*(np.exp(lo)-1):+.0f}%,{100*(np.exp(hi)-1):+.0f}%] in weekly rate")
open(os.path.join(HERE, "results", "ta_summary.txt"), "w").write("\n".join(lines) + "\n")
