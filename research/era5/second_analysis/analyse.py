"""RA-27 analysis, exactly as PREREGISTRATION.md specifies. Reads the per-time CSVs written by extract.py.

Run: python3 -I research/era5/second_analysis/analyse.py <repo_root> <floor_kt> <work_dir> [--synthetic]
Writes results/analysis_<floor>.txt and results/track_index_<floor>.csv.gz.
"""
import sys, os, glob, numpy as np, pandas as pd
root, floor, work = sys.argv[1], int(sys.argv[2]), sys.argv[3]
HERE = os.path.join(root, "research/era5/second_analysis")
HF_CUT = 71.7
NB, SEED = 5000, 20261008
S = pd.read_csv(os.path.join(HERE, f"results/sample_fixes_{floor}.csv"))
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(work, "*.csv")))], ignore_index=True)
F = S.merge(D, on=["track", "time"], how="left")
L = []
def out(s=""): L.append(s); print(s)

out(f"RA-27 second analysis system, floor {floor} kt. ERA5 and IFS HRES t0 analysis 10 m sustained wind; pipeline A tracks, seasons 2016-2021. Proxy.")
nmiss = int(F.wE_max.isna().sum())
out(f"fixes {len(F)}; fixes with no owned ocean cell or no match {nmiss}; fixes with ERA5 match > 100 km: {int((F.match_km > 100).sum())}")
F = F.dropna(subset=["wE_max", "wH_max"])

# track level index: max over the track's 00/12 fixes in the sample
agg = F.groupby("track").agg(basin=("basin", "first"), season=("season", "first"), gust=("g800", "max"), idx1200=("idx1200", "first"),
                              E=("wE_max", "max"), H=("wH_max", "max"), E99=("wE_p99", "max"), H99=("wH_p99", "max"), nfix=("time", "size")).reset_index()
agg.to_csv(os.path.join(HERE, f"results/track_index_{floor}.csv.gz"), index=False)
out(f"tracks {len(agg)} (Atlantic {int((agg.basin=='atl').sum())}, Pacific {int((agg.basin=='pac').sum())}); ERA5 gust index >= {HF_CUT}: {int((agg.idx1200>=HF_CUT).sum())}")

def match_cut(x, n):
    """Cut so that exactly n tracks (ties aside) have x >= cut."""
    n = int(max(1, min(n, len(x))))
    return np.sort(x)[::-1][n - 1]

def hss(a, b):
    a = np.asarray(a, bool); b = np.asarray(b, bool)
    tp = (a & b).sum(); fp = (~a & b).sum(); fn = (a & ~b).sum(); tn = (~a & ~b).sum()
    den = (tp + fn) * (fn + tn) + (tp + fp) * (fp + tn)
    return 2.0 * (tp * tn - fp * fn) / den if den else np.nan, (tp, fp, fn, tn)

def agreement(df, col_a, col_b, n_hf):
    ca, cb = match_cut(df[col_a].values, n_hf), match_cut(df[col_b].values, n_hf)
    return hss(df[col_a].values >= ca, df[col_b].values >= cb), ca, cb

def prim(df, col_a="E", col_b="H"):
    n_hf = int((df.idx1200 >= HF_CUT).sum())
    (h, tab), ca, cb = agreement(df, col_a, col_b, n_hf)
    return h, tab, ca, cb, n_hf

rng = np.random.default_rng(SEED)
def boot(df, fn, by=None):
    n = len(df); v = np.empty(NB); arr = df.reset_index(drop=True)
    if by is None:
        idx = rng.integers(0, n, (NB, n))
        for i in range(NB): v[i] = fn(arr.iloc[idx[i]])
    else:
        keys = arr[by].unique(); grp = {k: arr[arr[by] == k] for k in keys}
        for i in range(NB):
            pick = rng.choice(keys, len(keys), replace=True); v[i] = fn(pd.concat([grp[k] for k in pick]))
    return v

def hs_fn(col_a="E", col_b="H"):
    return lambda d: prim(d, col_a, col_b)[0]

tests = []
out("\n== Primary: HSS of crossing, ERA5 vs HRES t0, sustained 10 m wind, count-matched cuts (matched to the number of tracks with ERA5 gust index >= 71.7) ==")
for name, d in (("pooled", agg), ("Atlantic", agg[agg.basin == "atl"]), ("Pacific", agg[agg.basin == "pac"])):
    h, tab, ca, cb, n_hf = prim(d)
    b = boot(d, hs_fn()); lo, hi = np.nanpercentile(b, [2.5, 97.5])
    p = float((b <= 0.6).mean()); se = float(np.nanstd(b))
    tests.append((name, h, p, se, b))
    out(f"{name}: tracks {len(d)}, crossing {n_hf}, cuts ERA5 {ca:.1f} kt / HRES {cb:.1f} kt, 2x2 (both, HRES only, ERA5 only, neither) = {tab}; HSS {h:.3f}, 95% CI [{lo:.3f}, {hi:.3f}] (track bootstrap, cuts re-matched), SE {se:.3f}; one-sided p(H0 HSS <= 0.6) = {p:.4f}")
    if name == "pooled":
        bs = boot(d, hs_fn(), by="season"); out(f"  season-block bootstrap (6 seasons) 95% CI [{np.nanpercentile(bs,2.5):.3f}, {np.nanpercentile(bs,97.5):.3f}]")
        h1, se1 = h, se
# basin difference
hA = prim(agg[agg.basin == "atl"])[0]; hP = prim(agg[agg.basin == "pac"])[0]
bA = tests[1][4]; bP = tests[2][4]; diff = hA - hP; bd = bA - bP
p4 = float(2 * min((bd <= 0).mean(), (bd >= 0).mean()))
out(f"basin difference (Atlantic - Pacific) {diff:+.3f}, 95% CI [{np.nanpercentile(bd,2.5):+.3f}, {np.nanpercentile(bd,97.5):+.3f}], two-sided p {p4:.4f}")
ps = np.array([tests[0][2], tests[1][2], tests[2][2], p4])
order = np.argsort(ps); m = len(ps); q = np.empty(m); prev = 1.0
for r, i in list(enumerate(order))[::-1]:
    prev = min(prev, ps[i] * m / (r + 1)); q[i] = prev
out("\nBH-FDR over the four primary tests (T1 pooled, T2 Atlantic, T3 Pacific: H0 HSS <= 0.6; T4 basin difference): " + ", ".join(f"T{i+1} p {ps[i]:.4f} q {q[i]:.4f}" for i in range(4)) + f"; passing q < 0.05: {int((q<0.05).sum())} of 4")
mde = 0.6 + (1.645 + 0.84) * se1
out(f"power: pooled bootstrap SE {se1:.3f}; one-sided 5% test of HSS <= 0.6 has 80% power against a true HSS of {mde:.3f} (normal approximation)")
out(f"decision statistic: pooled HSS {h1:.3f}, 95% CI [{np.nanpercentile(tests[0][4],2.5):.3f}, {np.nanpercentile(tests[0][4],97.5):.3f}]")

out("\n== Secondary (descriptive; intervals are track bootstrap) ==")
# S2 agreement curve
h, tab, ca, cb, n_hf = prim(agg)
agg["rE"] = agg.E / ca; agg["Hx"] = agg.H >= cb
bins = [0, 0.85, 0.95, 1.05, 1.15, 10]
agg["bin"] = pd.cut(agg.rE, bins, right=False)
cur = agg.groupby("bin", observed=True).agg(n=("Hx", "size"), frac_HRES_crosses=("Hx", "mean")).reset_index()
out("S2 P(HRES crosses | ERA5 index relative to its cut):\n" + cur.to_string(index=False))
# S3 references
agg["Gx"] = agg.gust
for lab, a, b in (("S3a ERA5 sustained vs ERA5 gust (same system, different index)", "E", "gust"), ("S3b HRES sustained vs ERA5 gust (cross system, cross index)", "H", "gust")):
    n_hf = int((agg.idx1200 >= HF_CUT).sum()); (hh, tab), ca2, cb2 = agreement(agg, a, b, n_hf)
    bb = boot(agg, hs_fn(a, b)); out(f"{lab}: HSS {hh:.3f} CI [{np.nanpercentile(bb,2.5):.3f}, {np.nanpercentile(bb,97.5):.3f}], 2x2 {tab}")
# S4 band around the gust cut
bd_ = agg[(agg.idx1200 >= 65) & (agg.idx1200 <= 80)]
out(f"S4 tracks with ERA5 gust index 65-80 kt (n {len(bd_)}), cuts as in the primary ({ca:.1f}/{cb:.1f} kt): HSS {hss(bd_.E.values>=ca, bd_.H.values>=cb)[0]:.3f}, agreement share {((bd_.E.values>=ca)==(bd_.H.values>=cb)).mean():.3f}")
# S5 fix level
n_fx_hf = int((F.g800 >= HF_CUT).sum()); ce, ch = match_cut(F.wE_max.values, n_fx_hf), match_cut(F.wH_max.values, n_fx_hf)
hh, tab = hss(F.wE_max.values >= ce, F.wH_max.values >= ch)
out(f"S5 fix level (fixes with ERA5 gust >= 71.7: {n_fx_hf}; cuts {ce:.1f}/{ch:.1f} kt): HSS {hh:.3f}, 2x2 {tab}")
# S6 relation
sl, ic = np.polyfit(agg.E, agg.H, 1); res = agg.H - (sl * agg.E + ic)
near = agg[(agg.rE >= 0.85) & (agg.rE <= 1.15)]
out(f"S6 HRES index = {sl:.3f} x ERA5 index + {ic:.2f} kt; correlation {np.corrcoef(agg.E, agg.H)[0,1]:.3f}; residual SD {res.std():.2f} kt; SD of ERA5 index among tracks within 15% of its cut {near.E.std():.2f} kt (n {len(near)}); mean HRES - ERA5 {np.mean(agg.H-agg.E):+.2f} kt; ratio of 99th percentiles of the track indices {np.percentile(agg.H,99)/np.percentile(agg.E,99):.3f}")
# S7 p99 variant
n_hf = int((agg.idx1200 >= HF_CUT).sum()); (hh, tab), ca7, cb7 = agreement(agg, "E99", "H99", n_hf)
bb = boot(agg, hs_fn("E99", "H99")); out(f"S7 99th-percentile index instead of the maximum: HSS {hh:.3f}, CI [{np.nanpercentile(bb,2.5):.3f}, {np.nanpercentile(bb,97.5):.3f}]")
# S8 sample floor 60 within this sample
s60 = agg[agg.gust >= 60]; h60 = prim(s60)[0]; out(f"S8 sample restricted to ERA5 gust index >= 60 kt (n {len(s60)}): HSS {h60:.3f}")
open(os.path.join(HERE, f"results/analysis_{floor}.txt"), "w").write("\n".join(L) + "\n")
