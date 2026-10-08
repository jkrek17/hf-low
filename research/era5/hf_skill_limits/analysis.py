"""Error analysis of the PR 12 P(HF within 24 h) model (pipeline A, ERA5 proxy). Plan: PREREG.md.

usage: analysis.py REPO_ROOT OUTDIR
Reads committed tables only (intensity/results/fixes_2004, env_2004, hf_history/results/all_tracks.csv.gz,
data/hf_lows/HF_Data_-_*.csv). No ERA5 pull. Seasons 2004-05..2025-26, leave-one-season-out.
"""
import os, sys, json, warnings
import numpy as np, pandas as pd
from scipy.stats import norm
warnings.filterwarnings("ignore")
root, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(root, "research/era5/intensity"))
import model as M
from sklearn.linear_model import LinearRegression
os.makedirs(out, exist_ok=True)
W = os.path.join(root, "research/era5/hf_skill_limits/work"); os.makedirs(W, exist_ok=True)
HF = 71.7
NB = 1000
rng = np.random.default_rng(20261010)
lines = []
def w(s=""):
    lines.append(s); print(s, flush=True)

# ---------------------------------------------------------------- data
res = os.path.join(root, "research/era5/intensity/results")
D = M.load(os.path.join(res, "fixes_2004.csv.gz"), os.path.join(res, "env_2004.csv.gz"))
D["hf24"] = D.hf24.astype(int); D["hf_now"] = D.hf_now.astype(int)
D = D.sort_values(["track", "time"]).reset_index(drop=True)
# 12-hourly future values on the same track
D["dt"] = pd.to_datetime(D.time, format="%Y%m%d%H")
key = D.set_index(["track", "dt"])
def future(col, h):
    idx = pd.MultiIndex.from_arrays([D.track, D.dt + pd.Timedelta(hours=h)])
    return key[col].reindex(idx).values
for h in (12, 24, 36, 48):
    D[f"g{h}"] = future("g800", h); D[f"m{h}"] = future("msl", h)
D["dmsl12"] = D.m12 - D.msl
# peak gust over the two 12-hourly fixes in (t, t+24]; 06/18 UTC fixes are not in the committed table
D["gpk24"] = np.fmax(D.g12, D.g24)
D["gpk24"] = D.gpk24.fillna(0.0)
D["pk_consistent"] = ((D.gpk24 >= HF).astype(int) == D.hf24)
tr_all = pd.read_csv(os.path.join(root, "research/era5/hf_history/results/all_tracks.csv.gz"))
trk_hf = dict(zip(tr_all.track, tr_all.gust800_kt >= HF))
D["track_hf"] = D.track.map(trk_hf).fillna(False)
w(f"Pipeline A, ERA5 proxy. {len(D)} fixes, {D.season.nunique()} seasons {D.season.min()}-{D.season.max()}, hf24 base rate {D.hf24.mean():.4f}, "
  f"hf_now {D.hf_now.mean():.4f}, onset (not HF at t) {(1 - D.hf_now).sum()} fixes")
pos = D.hf24 == 1
w(f"hf24 positives {pos.sum()}; peak of the two 12-hourly fixes in the window is >= {HF} kt for {(pos & (D.gpk24 >= HF)).sum()} "
  f"({(pos & (D.gpk24 >= HF)).sum() / pos.sum():.3f}); the rest reach HF only at a 06/18 UTC point (not in the committed table)")

# ---------------------------------------------------------------- PR 12 model, leave-one-season-out
M.SETS = {k: M.SETS[k] for k in ("state", "full")}
S, P, CUT = M.loso(D, "hf24", 2)
assert len(S) == len(D)
D["p"] = P["full"][:, 1]; D["p_state"] = P["state"][:, 1]; D["cut"] = CUT["full"][:, 1]
D["fc"] = D.p >= D.cut
D["y"] = D.hf24
D["brier"] = (D.p - D.y) ** 2
D["pclim"] = P["clim"][:, 1]
D["brier_clim"] = (D.pclim - D.y) ** 2
bss = 1 - D.brier.sum() / D.brier_clim.sum()
w(f"Refit check: full-model BSS vs basin-month climatology {bss:+.3f} (PR 12 reported +0.423)")

def tab(fc, ob):
    return M.table(np.asarray(fc, bool), np.asarray(ob, bool))

# season-bootstrap helpers: per-season contingency counts for a mask
seasons = np.sort(D.season.unique())
Bidx = rng.integers(0, len(seasons), (NB, len(seasons)))
def season_counts(mask):
    m = mask.values if hasattr(mask, "values") else mask
    a = np.zeros((len(seasons), 4))
    for i, s in enumerate(seasons):
        g = D[(D.season == s).values & m]
        fc, ob = g.fc.values, g.y.values.astype(bool)
        a[i] = [np.sum(fc & ob), np.sum(fc & ~ob), np.sum(~fc & ob), np.sum(~fc & ~ob)]
    return a
def hss_from(a):
    A, Bq, C, Dd = a.sum(-2).T if a.ndim == 2 else a.T
    n = A + Bq + C + Dd
    ex = ((A + Bq) * (A + C) + (C + Dd) * (Bq + Dd)) / np.maximum(n, 1)
    return (A + Dd - ex) / np.maximum(n - ex, 1e-9)
def hss_boot(a):
    return np.array([hss_from(a[b]) for b in Bidx])

# ---------------------------------------------------------------- Q1 strata
D["onset"] = D.hf_now == 0
D["dist"] = pd.cut(D.g800, [-1, 55, 65, HF - 1e-9, 1e9], labels=["g<55", "55-65", "65-71.7", ">=71.7"])
D["monthgrp"] = np.select([D.month.isin([10, 11]), D.month.isin([12, 1, 2]), D.month.isin([3, 4])], ["Oct-Nov", "Dec-Feb", "Mar-Apr"], "May-Sep")
D["agegrp"] = pd.cut(D.age, [-1, 23.9, 72, 1e9], labels=["<24h", "24-72h", ">72h"])
D["ndrcls"] = np.where(D.cls >= 0, np.array(M.CLASSES, dtype=object)[np.clip(D.cls.values.astype(int), 0, 4)], "track ends / no class")
D["atl60"] = (D.basin == "atl") & (D.lat >= 60)
strata = [("hf_now", D.hf_now == 1), ("onset", D.onset),
          ("Atlantic", D.basin == "atl"), ("Pacific", D.basin == "pac"),
          ("Atlantic >= 60N", D.atl60), ("Atlantic < 60N", (D.basin == "atl") & ~D.atl60)]
strata += [(f"month {m}", D.monthgrp == m) for m in ["Oct-Nov", "Dec-Feb", "Mar-Apr", "May-Sep"]]
strata += [(f"g800 {d}", D.dist == d) for d in ["g<55", "55-65", "65-71.7", ">=71.7"]]
strata += [(f"24h change: {c}", D.ndrcls == c) for c in M.CLASSES + ["track ends / no class"]]
strata += [(f"age {a}", D.agegrp == a) for a in ["<24h", "24-72h", ">72h"]]
tot_loss = D.brier.sum(); nmiss = ((~D.fc) & (D.y == 1)).sum(); nfa = (D.fc & (D.y == 0)).sum()
w("")
w(f"== Q1 error by stratum (count-matched cut; {nmiss} misses, {nfa} false alarms, total Brier loss {tot_loss:.0f})")
w("stratum | share of fixes | share of Brier loss | share of misses | share of false alarms | base rate | POD | FAR | HSS [90% CI over seasons] | HSS vs rest, p")
rows = []
for name, m in strata:
    m = np.asarray(m)
    g = D[m]
    t = tab(g.fc, g.y)
    a = season_counts(m); bh = hss_boot(a)
    rest = ~m
    ar = season_counts(rest); bhr = hss_boot(ar)
    d = bh - bhr
    pv = min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()) + 1 / NB)
    rows.append(dict(stratum=name, n=int(m.sum()), share_fixes=m.mean(), share_brier=g.brier.sum() / tot_loss,
                     share_miss=((~g.fc) & (g.y == 1)).sum() / nmiss, share_fa=(g.fc & (g.y == 0)).sum() / nfa,
                     base=g.y.mean(), pod=t["pod"], far=t["far"], hss=t["hss"],
                     hss_lo=np.percentile(bh, 5), hss_hi=np.percentile(bh, 95), hss_rest=hss_from(ar), diff_p=pv))
R = pd.DataFrame(rows)
# Benjamini-Hochberg over the contrast family (strata that are not complements of each other are all kept; reported as a family)
o = np.argsort(R.diff_p.values); q = np.empty(len(R)); prev = 1.0
for rank, i in reversed(list(enumerate(o, 1))):
    prev = min(prev, R.diff_p.values[i] * len(R) / rank); q[i] = prev
R["diff_q"] = q
R.round(4).to_csv(os.path.join(out, "strata.csv"), index=False)
for r in R.itertuples():
    w(f"{r.stratum:28s} | {r.share_fixes:6.3f} | {r.share_brier:6.3f} | {r.share_miss:6.3f} | {r.share_fa:6.3f} | {r.base:6.3f} | "
      f"{r.pod:4.2f} | {r.far:4.2f} | {r.hss:5.2f} [{r.hss_lo:5.2f}, {r.hss_hi:5.2f}] | rest {r.hss_rest:5.2f} p {r.diff_p:.3f} q {r.diff_q:.3f}" if not np.isnan(r.pod) else f"{r.stratum}: no positives")

# H1: knife-edge share. misses and false alarms; peak within 6 kt of the cut
err = D[(D.fc & (D.y == 0)) | ((~D.fc) & (D.y == 1))].copy()
def knife(df, imputed):
    pk = df.gpk24.values.copy()
    if imputed:   # positives whose peak sits at a 06/18 point: the true peak is >= cut, put it on the cut
        pk = np.where((df.y.values == 1) & (pk < HF), HF, pk)
    return np.abs(pk - HF) <= 6.0
w("")
w("== H1 knife-edge: share of misses + false alarms whose 24 h peak gust is within 6 kt of 71.7")
for label, imputed in (("a) 12-hourly peak as is (positives reached only at a 06/18 point count by their 12-hourly peak)", False),
                       ("b) positives reached only at a 06/18 point placed on the cut", True)):
    sh = []
    for b in Bidx:
        num = den = 0
        for i in b:
            e = err[err.season == seasons[i]]
            k = knife(e, imputed); num += k.sum(); den += len(k)
        sh.append(num / den)
    k = knife(err, imputed)
    w(f"  {label}: {k.mean():.3f} [{np.percentile(sh, 5):.3f}, {np.percentile(sh, 95):.3f}]  (n errors {len(err)})")
    if not imputed: h1a = (k.mean(), np.percentile(sh, 5))
    else: h1b = (k.mean(), np.percentile(sh, 5))
for kind, m in (("misses", (~err.fc) & (err.y == 1)), ("false alarms", err.fc & (err.y == 0))):
    e = err[m.values]
    w(f"  {kind} n={len(e)}: within 6 kt {knife(e, True).mean():.3f} (b); peak gust quartiles (12-hourly) {np.percentile(e.gpk24, [25, 50, 75]).round(1).tolist()}")
fa = err[err.fc & (err.y == 0)]
w(f"  false alarms: HF reached 24-48 h later (hf48) {fa.hf48.mean():.3f}; storm reaches HF at some time {fa.track_hf.mean():.3f}; never HF {1 - fa.track_hf.mean():.3f}")
mi = err[(~err.fc) & (err.y == 1)]
w(f"  misses: already HF at t {mi.hf_now.mean():.3f}; storm-onset misses (not HF at t) {(1 - mi.hf_now.mean()):.3f}")
# H2
on = D.onset
w("")
w(f"== H2 onset share of Brier loss: {D.brier[on].sum() / tot_loss:.3f} ({on.mean():.3f} of fixes)")
sh = []
for b in Bidx:
    n_ = d_ = 0
    for i in b:
        g = D[D.season == seasons[i]]; n_ += g.brier[g.onset].sum(); d_ += g.brier.sum()
    sh.append(n_ / d_)
w(f"   90% CI over seasons [{np.percentile(sh, 5):.3f}, {np.percentile(sh, 95):.3f}]")
h2 = (D.brier[on].sum() / tot_loss, np.percentile(sh, 5))
# BSS split
o_ = D[on]; i_ = D[~on]
w(f"   BSS vs climatology: all {bss:+.3f}; onset fixes {1 - o_.brier.sum() / o_.brier_clim.sum():+.3f}; already-HF fixes {1 - i_.brier.sum() / i_.brier_clim.sum():+.3f}")
D[["track", "time", "basin", "season", "y", "p", "cut", "g800", "gpk24"]].round(4).to_csv(os.path.join(W, "loso_probs.csv.gz"), index=False)

# ---------------------------------------------------------------- Q2 oracle (uses the future)
w("")
w("== Q2 oracle fits (use the realised future; diagnostic, not a forecast). Rows: track survives 24 h and |NDR| <= 3 (class available)")
full = M.SETS["full"]
sub = D[(D.cls >= 0) & D.dmsl12.notna()].reset_index(drop=True)
sub["ndr2"] = sub.ndr24 ** 2
sub["dmsl12_2"] = sub.dmsl12 ** 2
def loso_bss(df, cols, label):
    pp = np.zeros(len(df)); 
    for s in seasons:
        te = (df.season == s).values; trn = ~te
        prep = M.Prep(df.loc[trn, cols].values.astype(float))
        m = M.fit(prep(df.loc[trn, cols].values.astype(float)), df.loc[trn, "hf24"].values)
        pp[te] = m.predict_proba(prep(df.loc[te, cols].values.astype(float)))[:, 1]
    return pp
pc = np.zeros(len(sub))
for s in seasons:
    te = (sub.season == s).values
    trn = sub[~te]
    for (b, mo), g in sub[te].groupby(["basin", "month"]):
        t_ = trn[(trn.basin == b) & (trn.month == mo)].hf24
        pc[g.index] = (t_.sum() + .5) / (len(t_) + 1)
yy = sub.hf24.values
bclim = ((pc - yy) ** 2).sum()
out2 = {}
for label, cols in (("base (PR 12 full)", full), ("+ realised 12 h pressure change", full + ["dmsl12", "dmsl12_2"]),
                    ("+ realised 24 h deepening rate", full + ["ndr24", "ndr2"]),
                    ("+ both", full + ["dmsl12", "dmsl12_2", "ndr24", "ndr2"])):
    pp = loso_bss(sub, cols, label)
    bs = ((pp - yy) ** 2)
    # bootstrap over seasons
    sums = np.array([bs[(sub.season == s).values].sum() for s in seasons]); cs = np.array([((pc - yy) ** 2)[(sub.season == s).values].sum() for s in seasons])
    bb = 1 - sums[Bidx].sum(1) / cs[Bidx].sum(1)
    out2[label] = (1 - bs.sum() / bclim, np.percentile(bb, 5), np.percentile(bb, 95), pp)
    w(f"  {label:36s} BSS {out2[label][0]:+.3f} [{out2[label][1]:+.3f}, {out2[label][2]:+.3f}]  AUC {M.roc_auc_score(yy, pp):.3f}")
base_bss = out2["base (PR 12 full)"][0]
sums_b = np.array([((out2['base (PR 12 full)'][3] - yy) ** 2)[(sub.season == s).values].sum() for s in seasons])
sums_o = np.array([((out2['+ realised 24 h deepening rate'][3] - yy) ** 2)[(sub.season == s).values].sum() for s in seasons])
cs = np.array([((pc - yy) ** 2)[(sub.season == s).values].sum() for s in seasons])
gain = (sums_b.sum() - sums_o.sum()) / cs.sum(); gb = (sums_b[Bidx].sum(1) - sums_o[Bidx].sum(1)) / cs[Bidx].sum(1)
w(f"  H3 oracle (b) gain over base: {gain:+.3f} [{np.percentile(gb, 5):+.3f}, {np.percentile(gb, 95):+.3f}] (rule: >= 0.10)")
h3 = (gain, np.percentile(gb, 5))
# how much of the 24 h deepening the environment explains (linear R2, LOSO)
y_ = sub.ndr24.values; pr = np.zeros(len(sub))
for s in seasons:
    te = (sub.season == s).values
    prep = M.Prep(sub.loc[~te, full].values.astype(float))
    lr = LinearRegression().fit(prep(sub.loc[~te, full].values.astype(float)), y_[~te])
    pr[te] = lr.predict(prep(sub.loc[te, full].values.astype(float)))
r2 = 1 - ((y_ - pr) ** 2).sum() / ((y_ - y_.mean()) ** 2).sum()
w(f"  24 h deepening rate (Bergerons) explained by the full predictor set: LOSO R2 {r2:.3f} (all fixes; PR 12 class RPSS 0.305)")
hn = sub.hf_now == 0
y2 = y_[hn.values]; pr2 = pr[hn.values]
w(f"  same, fixes not yet HF: R2 {1 - ((y2 - pr2) ** 2).sum() / ((y2 - y2.mean()) ** 2).sum():.3f}")

# ---------------------------------------------------------------- Q3 knife-edge ceiling
w("")
w("== Q3 knife-edge: forecasting the 24 h peak gust (not-yet-HF fixes whose track survives 24 h and has both 12-hourly fixes)")
g3 = D[(D.hf_now == 0) & D.g12.notna() & D.g24.notna()].reset_index(drop=True)
ycont = g3.gpk24.values
pm = np.zeros(len(g3))
for s in seasons:
    te = (g3.season == s).values
    prep = M.Prep(g3.loc[~te, full].values.astype(float))
    lr = LinearRegression().fit(prep(g3.loc[~te, full].values.astype(float)), ycont[~te])
    pm[te] = lr.predict(prep(g3.loc[te, full].values.astype(float)))
resid = ycont - pm
sig_now = resid.std()
w(f"  n {len(g3)}; peak gust mean {ycont.mean():.1f} kt, SD {ycont.std():.1f}; LOSO linear R2 {1 - (resid ** 2).sum() / ((ycont - ycont.mean()) ** 2).sum():.3f}; residual SD {sig_now:.1f} kt")
dens = ((np.abs(ycont - HF) <= 3).mean()) / 6
w(f"  fraction of fixes with peak gust within 3 kt of the cut: {(np.abs(ycont - HF) <= 3).mean():.3f}; HF base rate in this set {(ycont >= HF).mean():.4f}")
def bss_hss(pm_, truth, sigma, labels):
    p_ = norm.cdf((pm_ - HF) / sigma)
    base = labels.mean()
    bs = ((p_ - labels) ** 2).sum(); bc = ((base - labels) ** 2).sum()
    cut = np.quantile(p_, 1 - base)
    return 1 - bs / bc, M.table(p_ >= cut, labels.astype(bool))["hss"]
lab = (ycont >= HF).astype(int)
w("  BSS and HSS if the gust forecast error were smaller at fixed truth (forecast pulled toward the truth by k; Gaussian-error forecaster; climatology = not-yet-HF base rate)")
w("  k = 1 is the linear gust model, a consistency check on the logistic model's onset BSS 0.29 / HSS about 0.5")
q3 = []
for k in (1.0, 0.75, 0.5, 0.35, 0.25):
    pk = ycont + k * (pm - ycont)
    sg = k * sig_now
    b_, h_ = bss_hss(pk, ycont, sg, lab)
    q3.append((k, sg, b_, h_)); w(f"    k {k:4.2f}  error SD {sg:4.1f} kt  BSS {b_:+.3f}  HSS {h_:.2f}")
for sg in (2, 4, 6, 8, 10):
    k = sg / sig_now
    if k > 1.0: continue
    b_, h_ = bss_hss(ycont + k * (pm - ycont), ycont, sg, lab)
    w(f"    error SD {sg:2d} kt (k {k:4.2f})  BSS {b_:+.3f}  HSS {h_:.2f}")


# POST HOC (logged in PREREG.md): the homoscedastic Gaussian version above failed its own consistency check at k = 1
# (BSS +0.087 against the logistic model's onset +0.291), so the residual SD is also estimated by predicted-gust band.
w("")
w("== Q3b (post hoc) the same with a residual SD that depends on the predicted gust (bands of the linear prediction)")
edges = [-1e9, 30, 40, 50, 60, 70, 80, 1e9]
band = np.digitize(pm, edges) - 1
sig_b = np.array([resid[band == b].std() if (band == b).sum() > 50 else sig_now for b in range(len(edges) - 1)])
for b in range(len(edges) - 1):
    w(f"    predicted {edges[b]:>6.0f}..{edges[b+1]:<6.0f} n {(band == b).sum():6d}  residual SD {sig_b[b]:5.1f} kt  observed HF rate {lab[band == b].mean():.4f}")
sg_arr = sig_b[band]
for k in (1.0, 0.75, 0.5, 0.35, 0.25):
    pk = ycont + k * (pm - ycont); sg = k * sg_arr
    p_ = norm.cdf((pk - HF) / sg); base = lab.mean()
    bsv = 1 - ((p_ - lab) ** 2).sum() / ((base - lab) ** 2).sum()
    cutv = np.quantile(p_, 1 - base); hs = M.table(p_ >= cutv, lab.astype(bool))["hss"]
    w(f"    k {k:4.2f}  BSS {bsv:+.3f}  HSS {hs:.2f}   (k=1 should sit near the logistic onset BSS 0.29 / HSS 0.46)")
    q3.append((k, "band", bsv, hs))

# concurrent diagnosis: can the same fields tell whether the gust index is HF at this very time without the gust itself?
w("")
w("== Q2c concurrent diagnosis: HF now (g800 >= 71.7) from the full predictors WITHOUT g800, same fix (no forecasting; LOSO)")
nog = [c for c in full if c != "g800"]
D["now"] = D.hf_now
pp = np.zeros(len(D)); pcl_ = np.zeros(len(D))
for s in seasons:
    te = (D.season == s).values
    prep = M.Prep(D.loc[~te, nog].values.astype(float))
    m = M.fit(prep(D.loc[~te, nog].values.astype(float)), D.loc[~te, "now"].values)
    pp[te] = m.predict_proba(prep(D.loc[te, nog].values.astype(float)))[:, 1]
    trn = D[~te]
    for (b, mo), g in D[te].groupby(["basin", "month"]):
        t_ = trn[(trn.basin == b) & (trn.month == mo)].now
        pcl_[g.index] = (t_.sum() + .5) / (len(t_) + 1)
yn = D.now.values
bsn = (pp - yn) ** 2; bcn = (pcl_ - yn) ** 2
sums = np.array([bsn[(D.season == s).values].sum() for s in seasons]); cs_ = np.array([bcn[(D.season == s).values].sum() for s in seasons])
bb = 1 - sums[Bidx].sum(1) / cs_[Bidx].sum(1)
cut = np.quantile(pp, 1 - yn.mean()); t_ = M.table(pp >= cut, yn.astype(bool))
w(f"  base rate {yn.mean():.4f}  BSS {1 - bsn.sum() / bcn.sum():+.3f} [{np.percentile(bb, 5):+.3f}, {np.percentile(bb, 95):+.3f}]  AUC {M.roc_auc_score(yn, pp):.3f}  POD {t_['pod']:.2f} FAR {t_['far']:.2f} HSS {t_['hss']:.2f}")
q2c = (1 - bsn.sum() / bcn.sum(), t_["hss"])
# the same with the realised pressure at the target time as well (pressure history is the strongest single proxy for gust)

# ---------------------------------------------------------------- Q4 lead
w("")
w("== Q4 skill by lead (HF at exactly t+h among fixes not HF at t; full predictors; LOSO; BSS vs basin-month climatology)")
q4 = []
for h in (12, 24, 36, 48):
    d4 = D[(D.hf_now == 0)].copy().reset_index(drop=True)
    d4["yh"] = (d4[f"g{h}"].fillna(0) >= HF).astype(int)
    pp = np.zeros(len(d4)); pcl = np.zeros(len(d4))
    for s in seasons:
        te = (d4.season == s).values
        prep = M.Prep(d4.loc[~te, full].values.astype(float))
        m = M.fit(prep(d4.loc[~te, full].values.astype(float)), d4.loc[~te, "yh"].values)
        pp[te] = m.predict_proba(prep(d4.loc[te, full].values.astype(float)))[:, 1]
        trn = d4[~te]
        for (b, mo), g in d4[te].groupby(["basin", "month"]):
            t_ = trn[(trn.basin == b) & (trn.month == mo)].yh
            pcl[g.index] = (t_.sum() + .5) / (len(t_) + 1)
    yh = d4.yh.values
    bs = (pp - yh) ** 2; bcl = (pcl - yh) ** 2
    sums = np.array([bs[(d4.season == s).values].sum() for s in seasons]); cs_ = np.array([bcl[(d4.season == s).values].sum() for s in seasons])
    bb = 1 - sums[Bidx].sum(1) / cs_[Bidx].sum(1)
    cut = np.quantile(pp, 1 - yh.mean())
    t_ = M.table(pp >= cut, yh.astype(bool))
    q4.append((h, 1 - bs.sum() / bcl.sum(), t_["hss"]))
    w(f"  +{h:2d} h: base rate {yh.mean():.4f}  BSS {1 - bs.sum() / bcl.sum():+.3f} [{np.percentile(bb, 5):+.3f}, {np.percentile(bb, 95):+.3f}]  AUC {M.roc_auc_score(yh, pp):.3f}  POD {t_['pod']:.2f} FAR {t_['far']:.2f} HSS {t_['hss']:.2f}")

# ---------------------------------------------------------------- Q5, Q6 archive
w("")
w("== Q5 archive label (category HF only; archive fixes matched to the nearest same-time ERA5 fix of any track within 600 km; 06/18 UTC fixes matched to the nearest 00/12 UTC fix within 6 h and 800 km)")
A = []
for b, f in (("atl", "HF_Data_-_Atl.csv"), ("pac", "HF_Data_-_Pac.csv")):
    a = pd.read_csv(os.path.join(root, "data/hf_lows", f)); a["basin"] = b; A.append(a)
A = pd.concat(A)
A["dt"] = pd.to_datetime(A.date.astype(str), format="%Y%m%d%H", errors="coerce")
A = A[(A.Category == "HF") & A.dt.notna()]
A["season"] = np.where(A.dt.dt.month >= 6, A.dt.dt.year, A.dt.dt.year - 1)
A = A[(A.season >= 2004) & (A.season <= 2025)].reset_index(drop=True)
A["lon360"] = A.Longitude % 360
def hav(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2); dl = np.radians(lon2 - lon1)
    a_ = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * 3440.065 * 1.852 * np.arcsin(np.sqrt(a_))
byt = {t: g for t, g in D.groupby("time")}
match_tr, match_d, match_m = [], [], []
for r in A.itertuples():
    best = (None, 1e9, np.nan)
    for dh, lim in ((0, 600),) if r.dt.hour in (0, 12) else ((-6, 800), (6, 800)):
        t = (r.dt + pd.Timedelta(hours=dh)).strftime("%Y%m%d%H")
        if t in byt:
            g = byt[t]; dd = hav(r.Latitude, r.lon360, g.lat.values, g.lon.values)
            j = dd.argmin()
            if dd[j] < lim and dd[j] < best[1]:
                best = (g.track.values[j], dd[j], g.msl.values[j])
    match_tr.append(best[0]); match_d.append(best[1]); match_m.append(best[2])
A["track"], A["dist"], A["era5_msl"] = match_tr, match_d, match_m
nm = A.track.isna()
w(f"  archive HF fixes 2004-05..2025-26: {len(A)}; matched to a pipeline A track point {(~nm).sum()} ({(~nm).mean():.3f}); unmatched {nm.sum()}")
w(f"  matched tracks {A.track.nunique()}; archive events (IDs) {A.ID.nunique()}; events with at least one matched fix {A[~nm].ID.nunique()}")
Am = A[~nm].copy(); Am["track"] = Am.track.astype(int)
arch_by_track = {t: np.sort(g.dt.values) for t, g in Am.groupby("track")}
def arch24(r):
    ts = arch_by_track.get(r.track)
    if ts is None: return 0
    return int(((ts > np.datetime64(r.dt)) & (ts <= np.datetime64(r.dt + pd.Timedelta(hours=24)))).any())
D["arch24"] = [arch24(r) for r in D[["track", "dt"]].itertuples()]
D["arch_ever"] = D.track.isin(arch_by_track)
w(f"  archive-label positives {D.arch24.sum()} of {len(D)} fixes (ERA5 label {D.hf24.sum()}; both {((D.arch24 == 1) & (D.hf24 == 1)).sum()})")
def row(name, fc, ob):
    t = tab(fc, ob)
    w(f"  {name:62s} POD {t['pod']:.2f} FAR {t['far']:.2f} CSI {t['csi']:.2f} HSS {t['hss']:.2f} bias {t['bias']:.2f}")
    return t
c1 = row("ERA5 hf24 label as a forecast of the archive label (ceiling)", D.hf24 == 1, D.arch24 == 1)
row("   ... fixes with hf_now == 0 only", (D.hf24 == 1)[D.hf_now == 0], (D.arch24 == 1)[D.hf_now == 0])
row("PR 12 model at its ERA5-count cut, vs the archive label", D.fc, D.arch24 == 1)
cut_a = np.quantile(D.p, 1 - D.arch24.mean())
row("PR 12 model at the archive-count cut, vs the archive label", D.p >= cut_a, D.arch24 == 1)
row("PR 12 model vs the ERA5 label (as reported)", D.fc, D.hf24 == 1)
# bootstrap of ceiling HSS and model HSS vs archive
def boot_hss(fc, ob):
    a = np.array([[np.sum(fc[(D.season == s).values] & ob[(D.season == s).values]), np.sum(fc[(D.season == s).values] & ~ob[(D.season == s).values]),
                   np.sum(~fc[(D.season == s).values] & ob[(D.season == s).values]), np.sum(~fc[(D.season == s).values] & ~ob[(D.season == s).values])] for s in seasons])
    bh = hss_boot(a); return hss_from(a), np.percentile(bh, 5), np.percentile(bh, 95)
for name, fc, ob in (("ceiling (ERA5 label vs archive label)", (D.hf24 == 1).values, (D.arch24 == 1).values),
                     ("PR 12 model vs archive label", D.fc.values, (D.arch24 == 1).values),
                     ("PR 12 model vs ERA5 label", D.fc.values, (D.hf24 == 1).values)):
    h_, lo, hi = boot_hss(fc, ob); w(f"  HSS {name}: {h_:.3f} [{lo:.3f}, {hi:.3f}]")
q5 = c1
# Q6 pressure bias
w("")
w("== Q6 ERA5 central pressure minus archive pressure at matched archive HF fixes (same-time matches only, 00/12 UTC)")
Am["Pressure"] = pd.to_numeric(Am.Pressure, errors="coerce")
Ab = Am[Am.dt.dt.hour.isin([0, 12])].copy(); Ab["bias"] = Ab.era5_msl - Ab.Pressure
Ab = Ab[Ab.Pressure.between(900, 1030)]
Ab["band"] = pd.cut(Ab.Pressure, [0, 950, 960, 970, 980, 1040], labels=["<950", "950-960", "960-970", "970-980", ">=980"])
for bnd, g in Ab.groupby("band", observed=True):
    w(f"  archive {bnd:8s} n {len(g):5d}  ERA5 - archive: mean {g.bias.mean():+5.1f} hPa, median {g.bias.median():+5.1f}")
w(f"  all n {len(Ab)}: mean {Ab.bias.mean():+.1f}, median {Ab.bias.median():+.1f}; matched-distance median {Ab.dist.median():.0f} km")
# ---------------------------------------------------------------- summary of confirmatory tests
w("")
w("== Confirmatory tests (Holm over 3; bounds are the lower 5% season-bootstrap limits)")
w(f"  H1 knife-edge >= 60%: (a) {h1a[0]:.3f}, lower {h1a[1]:.3f}; (b) {h1b[0]:.3f}, lower {h1b[1]:.3f}  -> rule needs lower >= 0.50")
w(f"  H2 onset share of Brier loss >= 80%: {h2[0]:.3f}, lower {h2[1]:.3f}")
w(f"  H3 oracle deepening gain >= 0.10: {h3[0]:+.3f}, lower {h3[1]:+.3f}")
open(os.path.join(out, "skill_limits.txt"), "w").write("\n".join(lines) + "\n")
