"""Post hoc leakage / artifact audit of the storm-patch result (requested after PR 199). Existing patches and the local gust cache only; no pull.
usage: ERA5_WORK=DIR audit.py   -> results/audit.txt
 A fold hygiene: raw z500/msl channels (no climatology) with PCA refit per fold
 B neighbour lows: exclude fixes with another tracked low within 800 km at t or t+24 h (own position at that time)
 C tight-radius label: next-24 h peak ocean gust within 300 km (and 800 km as a replication of the stored label) from the gust cache
 D per-threshold decomposition of the RPSS gain; E per-season gains and interval methods
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import score_patch as SP
from score_patch import *
import numcodecs
from scipy import stats

S, X = SP.build()
y, season = S.gcls.values, S.season.values
Z = np.load(os.path.join(WORK, "oof_probs.npz"), allow_pickle=True)
L = {k: BS.rps(Z[k], y, 4) for k in ("clim", "base", "pc15")}
rng = np.random.default_rng(SP.SEED + 2)
lines = []
def out(t=""):
    print(t, flush=True); lines.append(t)
def gain(a, b, c, m=None):
    m = np.ones(len(S), bool) if m is None else m
    return BS.boot_gain(season[m], a[m], b[m], c[m], rng)
fm = np.ones((NP, NP), bool); cm = PR >= 300.0
SEASONS = sorted(np.unique(season))

# ---------- A fold hygiene
out("== A. Fold hygiene")
out("  Fitted inside each fold on the 7 training seasons only (checked in score_patch.py): channel scaling, patch mean, PCA, composite, analog library, M.Prep (median/clip/standardise), basin-month climatology baseline.")
out("  NOT fold-internal: the WeatherBench2 1990-2019 climatology used to form anomalies (an external, smooth 30-year mean that includes 2014-2019, six of the eight test seasons).")
raw = X[:, [4, 5]].astype(np.float32)          # raw z500 (m), raw msl (hPa): no climatology involved
def fold_raw(s):
    te = season == s; B = S[STATE].values.astype(float)
    sd = raw[~te].std((0, 2, 3), keepdims=True)
    Zf = (raw / sd).reshape(len(S), -1)
    mu = Zf[~te].mean(0)
    pca = PCA(NPC, svd_solver="full", random_state=0).fit(Zf[~te] - mu)
    pc = pca.transform(Zf - mu)
    return s, te, lr(np.c_[B, pc][~te], y[~te], np.c_[B, pc][te])
P = np.zeros((len(S), 4))
for s, te, p in Parallel(n_jobs=4)(delayed(fold_raw)(s) for s in SEASONS):
    P[te] = p
Lraw = BS.rps(P, y, 4)
pt, lo, hi, _ = gain(Lraw, L["base"], L["clim"])
out(f"  Raw z500 + msl patches only (no climatology, PCA per fold): gain over BASE {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]  (anomaly version: {gain(L['pc15'], L['base'], L['clim'])[0]:+.4f})")

# ---------- B neighbour lows
out("\n== B. Neighbour-low contamination")
F = pd.read_csv(os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz"), dtype={"time": str})
F["t"] = pd.to_datetime(F.time, format="%Y%m%d%H")
pos = {(r.track, r.t): (r.lat, r.lon) for r in F.itertuples()}
bytime = {t: g[["track", "lat", "lon"]].values for t, g in F.groupby("t")}
def hav(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2); d = p2 - p1; dl = np.radians(lo2 - lo1)
    a = np.sin(d / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))
def has_neighbour(track, t, lat, lon, lag):
    T = t + pd.Timedelta(hours=lag)
    la, lo = pos.get((track, T), (lat, lon))
    arr = bytime.get(T)
    if arr is None: return False
    m = arr[:, 0] != track
    return bool((hav(la, lo, arr[m, 1].astype(float), arr[m, 2].astype(float)) < 800).any())
tt = pd.to_datetime(S.time, format="%Y%m%d%H")
nb = {lag: np.array([has_neighbour(r.track, t, r.lat, r.lon, lag) for r, t in zip(S.itertuples(), tt)]) for lag in (0, 12, 24)}
excl = nb[0] | nb[24]
out(f"  Fixes with another tracked low within 800 km at t: {nb[0].mean():.3f}; at t+24 h: {nb[24].mean():.3f}; either (asked-for rule): {excl.mean():.3f} ({excl.sum()} fixes); HF outcomes among kept: {int((y[~excl]==3).sum())} of {int((y==3).sum())}")
keep = ~excl
pt, lo, hi, _ = gain(L["pc15"], L["base"], L["clim"], keep)
out(f"  Gain of PC15 over BASE, evaluated on the kept fixes only (models trained on all): {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
ix = np.flatnonzero(keep)
S2 = S.iloc[ix].reset_index(drop=True); X2 = X[ix]
res = Parallel(n_jobs=4)(delayed(fold)(s, S2, X2, X2, fm, cm) for s in SEASONS)
P2 = {k: np.zeros((len(S2), 4)) for k in res[0][2] if k != "analog_alone"}
for s, te, R in res:
    for k in P2: P2[k][te] = R[k]
L2 = {k: BS.rps(P2[k], S2.gcls.values, 4) for k in P2}
pt, lo, hi, _ = BS.boot_gain(S2.season.values, L2["pc15"], L2["base"], L2["clim"], rng)
out(f"  Retrained and scored on the kept fixes only (n {len(S2)}, HF outcomes {int((S2.gcls==3).sum())}): BASE {BS.skill(L2['base'], L2['clim']):+.3f}, PC15 {BS.skill(L2['pc15'], L2['clim']):+.3f}, gain {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
pt, lo, hi, _ = BS.boot_gain(S2.season.values, L2["s1_env"], L2["base"], L2["clim"], rng)
out(f"  Same, outer-pattern-only (r >= 300 km): gain {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
m = ~keep
pt, lo, hi, _ = gain(L["pc15"], L["base"], L["clim"], m)
out(f"  For reference, gain on the excluded fixes (neighbour present): {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
for lag in (12,):
    k3 = ~(nb[0] | nb[12] | nb[24])
    out(f"  Also excluding a neighbour at t+12: kept {k3.sum()} fixes; gain (evaluated on kept) {gain(L['pc15'], L['base'], L['clim'], k3)[0]:+.4f}")

# ---------- C label radius from the gust cache
out("\n== C. Label radius (gust cache, ERA5 proxy)")
CD = "/mnt/project-files/era5-cache"
lsm = np.load(os.path.join(CD, "lsm_nh_0p25.npy")); ocean = lsm < 0.5
need = {}
for t in tt:
    for lag in (12, 24): need[t + pd.Timedelta(hours=lag)] = None
H0 = np.datetime64("1900-01-01T00")
files = sorted(glob.glob(os.path.join(CD, "gust_nh_12h", "gust_*.npz")))
index = {}
for f in files:
    z = np.load(f); hs = z["hours_since_1900"]
    for i, h in enumerate(hs): index[pd.Timestamp(H0 + np.timedelta64(int(h), "h"))] = (f, i)
dec = numcodecs.Blosc()
cache = {}
def field(T):
    f, i = index[T]
    if f not in cache:
        cache.clear()
        z = np.load(f); cache[f] = np.frombuffer(dec.decode(z["data"].tobytes()), np.uint8).reshape(z["shape"])
    return cache[f][i].astype(np.float32) / 2.0
def gmax(g, la, lo, rkm):
    dl = rkm / 111.0 + 0.5; dlo = dl / max(np.cos(np.radians(min(abs(la), 80))), 0.2)
    r0, r1 = int(max(0, (90 - (la + dl)) / 0.25)), int(min(360, np.ceil((90 - (la - dl)) / 0.25)) + 1)
    cols = (np.arange(int(np.floor((lo - dlo) / 0.25)), int(np.ceil((lo + dlo) / 0.25)) + 1)) % 1440
    sub = g[r0:r1][:, cols]; oc = ocean[r0:r1][:, cols]
    la_g = 90 - 0.25 * np.arange(r0, r1)[:, None]; lo_g = 0.25 * cols[None, :]
    d = hav(la, lo, la_g, lo_g)
    mk = (d < rkm) & oc
    return float(sub[mk].max()) if mk.any() else 0.0
Tlist = sorted(need)
G = {}
byT = {}
for r, t in zip(S.itertuples(), tt):
    for lag in (12, 24):
        T = t + pd.Timedelta(hours=lag)
        if (r.track, T) in pos: byT.setdefault(T, []).append((r.Index, lag, *pos[(r.track, T)]))
lab = {800: np.zeros((len(S), 2)), 300: np.zeros((len(S), 2))}
for T in sorted(byT):
    if T not in index: continue
    g = field(T)
    for (i, lag, la, lo) in byT[T]:
        for R_ in (800, 300): lab[R_][i, lag // 12 - 1] = gmax(g, la, lo, R_)
for R_ in (800, 300): S[f"gm{R_}"] = np.maximum(lab[R_][:, 0], lab[R_][:, 1])
a = S.gm800.values; b = S.gmax24.values
out(f"  Replication of the stored label from the cache at 800 km: agreement of the 4-class outcome {np.mean(np.digitize(a, BS.CUTS) == y):.3f}; correlation of the peak gust {np.corrcoef(a, b)[0,1]:.3f} (stored label from the 0.25 degree tracker table)")
for R_ in (300,):
    yr = np.digitize(S[f"gm{R_}"].values, BS.CUTS)
    out(f"  {R_} km label class frequencies: " + ", ".join(f"{n} {f:.3f}" for n, f in zip(BS.NAMES, np.bincount(yr, minlength=4) / len(yr))) + f"; HF outcomes {int((yr==3).sum())}")
    S3 = S.copy(); S3["gcls"] = yr
    res = Parallel(n_jobs=4)(delayed(fold)(s, S3, X, X, fm, cm) for s in SEASONS)
    P3 = {k: np.zeros((len(S3), 4)) for k in res[0][2] if k != "analog_alone"}
    for s, te, R in res:
        for k in P3: P3[k][te] = R[k]
    L3 = {k: BS.rps(P3[k], yr, 4) for k in P3}
    pt, lo, hi, _ = BS.boot_gain(season, L3["pc15"], L3["base"], L3["clim"], rng)
    out(f"  Label = peak gust within {R_} km (retrained, same fixes): BASE {BS.skill(L3['base'], L3['clim']):+.3f}, PC15 {BS.skill(L3['pc15'], L3['clim']):+.3f}, gain {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    pt, lo, hi, _ = BS.boot_gain(season, L3["s1_env"], L3["base"], L3["clim"], rng)
    out(f"    outer pattern only (r >= 300 km): gain {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    k_ = ~excl
    pt, lo, hi, _ = BS.boot_gain(season[k_], L3["pc15"][k_], L3["base"][k_], L3["clim"][k_], rng)
    out(f"    300 km label and no neighbour low (evaluated on kept): gain {pt:+.4f} [{lo:+.4f}, {hi:+.4f}]")

# ---------- D per-threshold
out("\n== D. Where the RPSS gain sits (RPS = sum of Brier scores of the three cumulative events; gain / sum RPS_clim)")
cum = lambda Pm: np.cumsum(Pm, 1)[:, :3]
Yc = np.cumsum(np.eye(4)[y], 1)[:, :3]
tot = L["clim"].sum()
names = ["below gale | at least gale (P(class>=1))", "gale-or-less | storm or more (class>=2)", "storm-or-less | HF (class>=3)"]
for j, n in enumerate(names):
    bb = (cum(Z["base"])[:, j] - Yc[:, j]) ** 2; pp = (cum(Z["pc15"])[:, j] - Yc[:, j]) ** 2
    out(f"  threshold {j+1} {n:40s} share of RPSS gain {(bb.sum()-pp.sum())/tot:+.4f}  (Brier skill vs clim: base {1-bb.sum()/((cum(Z['clim'])[:, j]-Yc[:, j])**2).sum():+.3f} -> pc15 {1-pp.sum()/((cum(Z['clim'])[:, j]-Yc[:, j])**2).sum():+.3f})")
out("  (the three shares add to the total RPSS gain; the RPSS gain is relative to the whole climatology RPS, so thresholds with more variance weigh more)")

# ---------- E seasons and intervals
out("\n== E. Per-season gains and interval methods")
g_s = np.array([(L["base"][season == s].sum() - L["pc15"][season == s].sum()) / L["clim"][season == s].sum() for g_ in [0] for s in SEASONS])
tot_s = np.array([L["clim"][season == s].sum() for s in SEASONS]); num = np.array([L["base"][season == s].sum() - L["pc15"][season == s].sum() for s in SEASONS])
out("  per-season gain: " + ", ".join(f"{s}-{str(s+1)[2:]} {g:+.4f}" for s, g in zip(SEASONS, g_s)))
out(f"  pooled {num.sum()/tot_s.sum():+.4f}; across seasons mean {g_s.mean():+.4f}, SD {g_s.std(ddof=1):.4f}, range {g_s.min():+.4f}..{g_s.max():+.4f}")
loo = [(num.sum() - num[i]) / (tot_s.sum() - tot_s[i]) for i in range(8)]
out(f"  leave-one-season-out pooled gain: {min(loo):+.4f} .. {max(loo):+.4f}")
se = g_s.std(ddof=1) / np.sqrt(8); tc = stats.t.ppf(0.95, 7)
out(f"  90% interval, method 1 (as reported): resample 8 seasons with replacement 2000 times, pool the per-season loss sums, 5-95 percentile: {gain(L['pc15'], L['base'], L['clim'])[1]:+.4f} .. {gain(L['pc15'], L['base'], L['clim'])[2]:+.4f}")
out(f"  90% interval, method 2: t interval on the 8 per-season gains (t7 = {tc:.2f}): {g_s.mean()-tc*se:+.4f} .. {g_s.mean()+tc*se:+.4f}")
out("  Both treat seasons as independent and fixed the model (15 PCs, k = 30, C = 1) before the data were read; neither includes uncertainty from the ERA5 proxy or from choices not tried. Within a season the 4,000 fixes are strongly autocorrelated, which is why seasons, not fixes, are resampled.")
open(os.path.join(OUT, "audit.txt"), "w").write("\n".join(lines) + "\n")
