"""RA-21 analysis of the held-out probabilities from run_loso.py. Pipeline A, ERA5 proxy. Rules are in PREREG.md.

usage: analyse.py REPO_ROOT WORKDIR OUTDIR
Reads WORKDIR/probs.npz and meta.csv (not committed). Writes OUTDIR/loso_results.txt, tests.csv, season_brier.csv, strata.csv.
"""
import os, sys
import numpy as np, pandas as pd
root, work, out = sys.argv[1:4]
sys.path.insert(0, os.path.join(root, "research/era5/intensity"))
import model as M
os.makedirs(out, exist_ok=True)
meta = pd.read_csv(os.path.join(work, "meta.csv"), dtype={"time": str})
Z = np.load(os.path.join(work, "probs.npz"))
P = {k: Z[k].astype(float) for k in Z.files}
y = meta.hf24.values.astype(float); seas = meta.season.values; basin = meta.basin.values
us = np.unique(seas)
rng = np.random.default_rng(20261010)
B = rng.integers(0, len(us), (1000, len(us)))
flip = np.where(np.random.default_rng(20261010).random((100000, len(us))) < .5, -1.0, 1.0)
lines = []
def w(s=""):
    lines.append(s); print(s)

def fmt_p(p):
    return '<1e-5' if p < 1e-5 else f'{p:.4f}'

def bs(p, m=None):
    v = (p - y) ** 2
    return v

def sums(v, mask=None):
    m = np.ones(len(v), bool) if mask is None else mask
    return np.array([v[m & (seas == s)].sum() for s in us])

def gain(a, b, ref, mask=None, seasons_mask=None):
    """Skill gain of model a over model b (a minus b in BSS units against ref) with season bootstrap and sign-flip p."""
    sa, sb, sr = sums(bs(P[a]), mask), sums(bs(P[b]), mask), sums(bs(P[ref]), mask)
    use = np.ones(len(us), bool) if seasons_mask is None else seasons_mask
    sa, sb, sr = sa * use, sb * use, sr * use
    d = (sb - sa) / sr.sum()
    est = d.sum()
    bt = (sb[B].sum(1) - sa[B].sum(1)) / sr[B].sum(1) if seasons_mask is None else None
    if bt is None:
        # resample only the seasons in use
        idx = np.where(use)[0]
        Bu = np.random.default_rng(20261010).integers(0, len(idx), (1000, len(idx)))
        bt = (sb[idx][Bu].sum(1) - sa[idx][Bu].sum(1)) / sr[idx][Bu].sum(1)
        dd = d[idx]
        p = np.mean(np.abs((np.where(np.random.default_rng(1).random((100000, len(idx))) < .5, -1.0, 1.0) * dd).sum(1)) >= abs(dd.sum()) - 1e-15)
    else:
        p = np.mean(np.abs((flip * d).sum(1)) >= abs(est) - 1e-15)
    npos = int((d[use] > 0).sum())
    return dict(est=est, lo=np.percentile(bt, 5), hi=np.percentile(bt, 95), se=bt.std(), p=p, npos=npos, nseas=int(use.sum()))

def skill(name, ref, mask=None):
    sm, sr = sums(bs(P[name]), mask), sums(bs(P[ref]), mask)
    return 1 - sm.sum() / sr.sum()

w("RA-21: LOSO check of the boosted-tree gain on P(HF within 24 h). Pipeline A, ERA5 proxy, 159,430 fixes, 22 seasons 2004-05 to 2025-26.")
w(f"n fixes {len(y)}, seasons {len(us)}, hf24 rate {y.mean():.4f}")
b0 = skill("M0", "clim")
w(f"Reproduction check: M0 logistic BSS {b0:+.4f} (PR 12 reports 0.423; allowed band 0.421-0.425) -> {'OK' if abs(b0 - 0.423) <= 0.002 else 'FAILED'}")
if abs(b0 - 0.423) > 0.002:
    w("STOP: the loader or fold design does not reproduce PR 12; no tree score is reported.")
    open(os.path.join(out, "loso_results.txt"), "w").write("\n".join(lines) + "\n"); sys.exit(1)
w("")
w("BSS against basin-month climatology, LOSO (pooled over 22 seasons):")
names = ["M0", "M1", "M2", "M3"]
for n in names:
    w(f"  {n}: {skill(n, 'clim'):+.4f}")
tests = []
def add(tid, desc, g):
    tests.append(dict(id=tid, desc=desc, **g))

# P1
g = gain("M1", "M0", "clim"); add("P1", "BSS(M1) - BSS(M0), LOSO, all 22 seasons", g)
w("")
w(f"P1  M1 (boosted) minus M0 (logistic): {g['est']:+.4f}  90% CI [{g['lo']:+.4f}, {g['hi']:+.4f}]  SE {g['se']:.4f}  sign-flip p {fmt_p(g['p'])}  seasons M1 better: {g['npos']}/22")
est, lo = g["est"], g["lo"]
if est >= 0.02 and lo > 0: verdict = "REAL, as predicted (>= +0.020, lower bound > 0)"
elif est >= 0.01 and lo > 0: verdict = "REAL, smaller than predicted (+0.010 to +0.020, lower bound > 0)"
elif est < 0.01: verdict = "NOT REAL: gain under +0.010" + ("" if g["hi"] < 0.02 else " (but the upper bound is above +0.020, so inconclusive)")
else: verdict = "CAN'T TELL (>= +0.010 but lower bound <= 0)"
w(f"    Pre-registered verdict: {verdict}")
w(f"    80% minimum detectable gain (2.8 x SE): {2.8 * g['se']:.4f}")
# P1b strata
cur = meta.g800.values
strata = {"<55 kt": cur < 55, "55-71.7 kt": (cur >= 55) & (cur < 71.7), ">=71.7 kt": cur >= 71.7}
tot_d = sums(bs(P["M0"])) - sums(bs(P["M1"]))
rows = []
w("")
w("P1b  Where the gain sits, by current g800 (Brier-sum reduction M0 -> M1):")
for k, m in strata.items():
    d = sums(bs(P["M0"]), m) - sums(bs(P["M1"]), m)
    share = d.sum() / tot_d.sum()
    bsh = d[B].sum(1) / tot_d[B].sum(1)
    dm0 = (bs(P["M0"])[m]).mean()
    dpf = (bs(P["M0"])[m] - bs(P["M1"])[m]).mean()
    sm = skill("M1", "clim", m) - skill("M0", "clim", m)
    rows.append(dict(stratum=k, share_fixes=m.mean(), n=int(m.sum()), hf24_rate=y[m].mean(), share_of_gain=share, share_lo=np.percentile(bsh, 5), share_hi=np.percentile(bsh, 95),
                     gain_per_fix=dpf, dBSS_in_stratum=sm, brier_M0=dm0, brier_M1=(bs(P["M1"])[m]).mean()))
    w(f"  {k:11s} fixes {m.mean() * 100:5.1f}%  hf24 rate {y[m].mean():.3f}  share of the gain {share * 100:5.1f}% [{np.percentile(bsh, 5) * 100:.1f}, {np.percentile(bsh, 95) * 100:.1f}]  "
      f"gain per fix {dpf:+.5f}  stratum BSS gain {sm:+.4f}")
pd.DataFrame(rows).to_csv(os.path.join(out, "strata.csv"), index=False)
mid = rows[1]
met = mid["share_of_gain"] >= 0.60 and mid["gain_per_fix"] > max(rows[0]["gain_per_fix"], rows[2]["gain_per_fix"])
w(f"    Stratum prediction (55-71.7 kt carries >= 60% of the gain and the largest per-fix gain): {'MET' if met else 'NOT MET'}")

# secondaries
add("S1", "onset only (not HF at t): M1 - M0", gain("M1", "M0", "clim", meta.hf_now.values == 0))
add("S2a", "Atlantic only: M1 - M0", gain("M1", "M0", "clim", basin == "atl"))
add("S2b", "Pacific only: M1 - M0", gain("M1", "M0", "clim", basin == "pac"))
# S3 temporal split
late = np.isin(seas, us[us > 2014])
sa, sb, sr = [((P["L2_" + n] - y) ** 2)[late] for n in ("M1", "M0", "clim")]
ul = us[us > 2014]
ssa = np.array([sa[seas[late] == s].sum() for s in ul]); ssb = np.array([sb[seas[late] == s].sum() for s in ul]); ssr = np.array([sr[seas[late] == s].sum() for s in ul])
Bl = np.random.default_rng(20261010).integers(0, len(ul), (1000, len(ul)))
d = (ssb - ssa) / ssr.sum(); btl = (ssb[Bl].sum(1) - ssa[Bl].sum(1)) / ssr[Bl].sum(1)
fl = np.where(np.random.default_rng(2).random((100000, len(ul))) < .5, -1.0, 1.0)
add("S3", "temporal split, train 2004-14, test 2015-25: M1 - M0", dict(est=d.sum(), lo=np.percentile(btl, 5), hi=np.percentile(btl, 95), se=btl.std(),
    p=np.mean(np.abs((fl * d).sum(1)) >= abs(d.sum()) - 1e-15), npos=int((d > 0).sum()), nseas=len(ul)))
w("")
w(f"S3 temporal split: M0 BSS {1 - ssb.sum() / ssr.sum():+.4f}, M1 BSS {1 - ssa.sum() / ssr.sum():+.4f} (PR 12 reported M0 0.426)")
# S4 HSS
def hss_counts(fc):
    ob = y.astype(bool)
    return np.array([[np.sum(fc & ob & (seas == s)), np.sum(fc & ~ob & (seas == s)), np.sum(~fc & ob & (seas == s)), np.sum(~fc & ~ob & (seas == s))] for s in us])
def hss_from(c):
    a, b_, c_, d_ = c.sum(0); n = a + b_ + c_ + d_
    ex = ((a + b_) * (a + c_) + (c_ + d_) * (b_ + d_)) / n
    return (a + d_ - ex) / (n - ex)
cut0 = np.quantile(P["M0"], 1 - y.mean()); cut1 = np.quantile(P["M1"], 1 - y.mean())
c0, c1 = hss_counts(P["M0"] >= cut0), hss_counts(P["M1"] >= cut1)
h0, h1 = hss_from(c0), hss_from(c1)
bh = np.array([hss_from(c1[b]) - hss_from(c0[b]) for b in B])
dseas = np.array([hss_from(c1[[i]]) - hss_from(c0[[i]]) for i in range(len(us))])
fs = np.where(np.random.default_rng(3).random((20000, len(us))) < .5, -1.0, 1.0)
# sign-flip on pooled HSS is not defined from per-season terms; use season-level HSS differences' mean
pS4 = np.mean(np.abs((fs * (dseas - 0)).mean(1)) >= abs(dseas.mean()) - 1e-15)
add("S4", "HSS(M1) - HSS(M0) at the pooled count-matched cut (p: sign flip on season-level HSS differences)", dict(est=h1 - h0, lo=np.percentile(bh, 5), hi=np.percentile(bh, 95), se=bh.std(), p=pS4, npos=int((dseas > 0).sum()), nseas=len(us)))
w(f"S4 HSS at the pooled count-matched cut: M0 {h0:.3f}, M1 {h1:.3f} (difference {h1 - h0:+.3f})")
tq0 = P["M0"] >= Z["cut_M0"]; tq1 = P["M1"] >= Z["cut_M1"]
w(f"   contrast, training-quantile cut (PR 12 rule): forecast/observed count M0 {tq0.sum() / y.sum():.2f}, M1 {tq1.sum() / y.sum():.2f}; HSS M0 {hss_from(hss_counts(tq0)):.3f}, M1 {hss_from(hss_counts(tq1)):.3f}")
# decomposition
add("S5", "additive trees M2 - M0 (single-ingredient nonlinearity)", gain("M2", "M0", "clim"))
add("S6", "M1 - M2 (interaction)", gain("M1", "M2", "clim"))
add("S7", "smooth additive splines M3 - M0", gain("M3", "M0", "clim"))
add("S8", "M1 - M3 (sharp thresholds and interactions beyond smooth curves)", gain("M1", "M3", "clim"))
for c in ["msl", "dp12", "lat", "speed", "g800", "logage", "B", "VTL", "VTU", "jet250", "div300", "vadv500", "eady", "sst", "sstgrad", "sst_t500", "flux", "tcwv"]:
    add("S9_" + c, f"M0 + spline of {c} minus M0", gain("S9_" + c, "M0", "clim"))
for k in P:
    if k.startswith("S10_"):
        add(k, "M0 + product " + k[4:].replace("_x_", " x ") + " minus M0", gain(k, "M0", "clim"))

T = pd.DataFrame(tests)
def bh_q(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q
sec = T.id != "P1"
T["q_secondary"] = np.nan; T.loc[sec, "q_secondary"] = bh_q(T.loc[sec, "p"])
T["q_all"] = bh_q(T.p)
T.round(5).to_csv(os.path.join(out, "tests.csv"), index=False)
w("")
w(f"All {len(T)} pre-registered tests ({sec.sum()} secondary). est = BSS-unit gain (HSS units for S4), 90% season-bootstrap CI, sign-flip p, BH q within the secondary family and across all.")
for _, r in T.iterrows():
    qs = "" if np.isnan(r.q_secondary) else f"{r.q_secondary:.3f}"
    w(f"  {r.id:22s} {r.est:+.4f} [{r.lo:+.4f}, {r.hi:+.4f}]  p {fmt_p(r.p)}  q_sec {qs:>6s}  q_all {r.q_all:.3f}  seasons+ {r.npos}/{r.nseas}  {r.desc}")
w(f"Tests with q_all < 0.05: {(T.q_all < .05).sum()} of {len(T)}; secondary q < 0.05: {(T.q_secondary < .05).sum()} of {sec.sum()}")

# per season table
rows = []
for s in us:
    m = seas == s
    r = dict(season=int(s), n=int(m.sum()), n_hf24=int(y[m].sum()))
    for n in P:
        if n.startswith("cut_") or n.startswith("L2_"): continue
        r["brier_" + n] = float(((P[n][m] - y[m]) ** 2).sum())
    rows.append(r)
pd.DataFrame(rows).round(4).to_csv(os.path.join(out, "season_brier.csv"), index=False)
open(os.path.join(out, "loso_results.txt"), "w").write("\n".join(lines) + "\n")
pr = meta[["track", "time", "season", "basin", "hf24", "hf_now", "g800"]].copy()
for n in ("clim", "M0", "M1", "M2", "M3"):
    pr[n] = P[n].round(5)
pr.to_csv(os.path.join(out, "loso_probs.csv.gz"), index=False)
