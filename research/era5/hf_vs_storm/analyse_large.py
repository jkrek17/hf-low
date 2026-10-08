"""Large-scale composites: HF onset / HF peak against storm-force-only peak (hf_vs_storm, PREREGISTRATION.md Part 1).

ERA5 proxy, pipeline A. 5.625 degree, 12 UTC daily fields from fields_large.py (season_<s>.npz in $ERA5_WORK/hf_vs_storm_large).
usage: analyse_large.py REPO_ROOT FIELDS_DIR [B]
Outputs results/large/: boxtests.csv, series.csv, maps_<basin>_<cmp>.npz, summary.txt
"""
import sys, os, glob, datetime as dt, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from strat import cells, boot, bh
root, fdir = sys.argv[1], sys.argv[2]
B = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
OUT = os.path.join(root, "research/era5/hf_vs_storm/results/large"); os.makedirs(OUT, exist_ok=True)
VARS = ["z500", "u250", "mslp", "sst", "tcwv"]
LAT = -87.1875 + 5.625 * np.arange(32); LAT = LAT[20:32]; LON = 5.625 * np.arange(64)
BOX = {"atl": (53.5, 318.9), "pac": (42.25, 171.5)}
SECT = {"atl": lambda lon: (lon >= 260) | (lon <= 40), "pac": lambda lon: (lon >= 110) & (lon <= 250)}
rng = np.random.default_rng(20261012)

# ---- fields and anomalies
dates, X = [], {v: [] for v in VARS}
for f in sorted(glob.glob(os.path.join(fdir, "season_*.npz"))):
    z = np.load(f)
    dates += [dt.date.fromisoformat(d) for d in z["dates"]]
    for v in VARS:
        X[v].append(z[v])
X = {v: np.concatenate(X[v]) for v in VARS}
didx = {d: i for i, d in enumerate(dates)}
mmdd = np.array([d.strftime("%m-%d") for d in dates])
order = [d.strftime("%m-%d") for d in (dt.date(2003, 8, 20) + dt.timedelta(days=i) for i in range((dt.date(2004, 5, 31) - dt.date(2003, 8, 20)).days + 1))]
A = {}
for v in VARS:
    clim = np.stack([np.nanmean(X[v][mmdd == k], 0) for k in order])
    sm = np.empty_like(clim)
    for i in range(len(order)):
        sm[i] = np.nanmean(clim[max(0, i - 15):i + 16], 0)
    pos = {k: i for i, k in enumerate(order)}
    A[v] = X[v] - sm[[pos[k] for k in mmdd]]
print("days", len(dates), "fields", {v: X[v].shape for v in VARS}, flush=True)

# ---- storms
S = pd.read_csv(os.path.join(root, "research/era5/hf_vs_storm/results/storms.csv"))
def day_for(anchor_dt, k):
    ref = anchor_dt - pd.Timedelta(days=int(k))
    return (ref if ref.hour >= 12 else ref - pd.Timedelta(days=1)).date()
KS = list(range(11))
def storm_table(kind):
    if kind == "HF_onset":
        d = S[S.grp == "HF"].copy(); d["t"] = d.onset_time
    elif kind == "HF_peak":
        d = S[S.grp == "HF"].copy(); d["t"] = d.peak_time
    else:
        d = S[S.grp == "SF"].copy(); d["t"] = d.peak_time
    d["t"] = pd.to_datetime(d.t.astype("int64").astype(str), format="%Y%m%d%H")
    idx = np.full((len(d), len(KS)), -1)
    for i, t in enumerate(d.t):
        for j, k in enumerate(KS):
            idx[i, j] = didx.get(day_for(t, k), -1)
    d["ok"] = (idx >= 0).all(1)
    return d, idx
cfgs = ["L", "k0", "k1", "k2", "k4", "k7"]
cfg_k = {"L": list(range(4, 11)), "k0": [0], "k1": [1], "k2": [2], "k4": [4], "k7": [7]}

def box_weights(b):
    la0, lo0 = BOX[b]
    dl = np.abs((LON - lo0 + 180) % 360 - 180)
    m = (np.abs(LAT[:, None] - la0) <= 15) & (dl[None, :] <= 30)
    w = np.cos(np.radians(LAT))[:, None] * m
    return w / w.sum(), m

rows_box, rows_ser, summ = [], [], []
for b in ("atl", "pac"):
    wbox, mbox = box_weights(b)
    sect = SECT[b](LON)
    for cmpname, hfkind in (("C1", "HF_onset"), ("C2", "HF_peak")):
        th, ih = storm_table(hfkind); ts, is_ = storm_table("SF_peak")
        # restrict to basin & complete lags
        mh = (th.basin.values == b); ms = (ts.basin.values == b)
        th, ih = th[mh], ih[mh]; ts, is_ = ts[ms], is_[ms]
        okh, oks = th.ok.values, ts.ok.values
        th, ih, ts, is_ = th[okh], ih[okh], ts[oks], is_[oks]
        summ.append(f"{b} {cmpname}: HF {len(th)} of {int(mh.sum())} (complete lags), SF {len(ts)} of {int(ms.sum())}")
        grp = np.r_[np.zeros(len(th), int), np.ones(len(ts), int)]
        season = np.r_[th.season.values, ts.season.values]
        stratum = np.r_[th.mon.values, ts.mon.values] - 1
        seasons = np.arange(2004, 2026)
        idx = np.vstack([ih, is_])
        # time series of box means, all lags
        def bmean(a):      # box mean over non-missing cells (SST is NaN over land), area weights renormalised
            m = ~np.isnan(a)
            return (np.where(m, a, 0) * wbox).sum((1, 2)) / (m * wbox).sum((1, 2))
        ser = np.stack([np.stack([bmean(A[v][idx[:, j]]) for v in VARS], 1) for j in range(len(KS))], 1)   # [n, 11, 5]
        sm, cn = cells(ser, grp, season, stratum, seasons, 12)
        obs, reps, p = boot(sm, cn, B, rng)
        lo, hi = np.nanpercentile(reps, [2.5, 97.5], axis=0)
        for j, k in enumerate(KS):
            for iv, v in enumerate(VARS):
                rows_ser.append(dict(basin=b, cmp=cmpname, var=v, k=k, diff=obs[0][j, iv], hf=obs[1][j, iv], sf=obs[2][j, iv],
                                     lo=lo[j, iv], hi=hi[j, iv], se=np.nanstd(reps[:, j, iv])))
        for cfg in cfgs:
            ks = cfg_k[cfg]
            d = reps[:, ks, :].mean(1)
            o = obs[0][ks, :].mean(0)
            for iv, v in enumerate(VARS):
                dd = d[:, iv]
                pp = min(1.0, 2 * min((np.sum(dd > 0) + 1) / (B + 1), (np.sum(dd < 0) + 1) / (B + 1)))
                rows_box.append(dict(basin=b, cmp=cmpname, cfg=cfg, var=v, diff=o[iv], hf=obs[1][ks, iv].mean(), sf=obs[2][ks, iv].mean(),
                                     lo=np.percentile(dd, 2.5), hi=np.percentile(dd, 97.5), se=dd.std(), p=pp,
                                     mde80=2.8 * dd.std(), n_hf=len(th), n_sf=len(ts)))
        # maps: composites for each config, anomalies for all vars plus raw z500/u250/mslp
        fields = {}
        for cfg in cfgs:
            ks = cfg_k[cfg]
            fields[cfg] = np.stack([A[v][idx[:, ks]].mean(1) for v in VARS] + [X[v][idx[:, ks]].mean(1) for v in ("z500", "u250", "mslp")], 1)
        mp = np.stack([fields[c] for c in cfgs], 1)       # [n, 6, 8, 12, 64]
        sm2, cn2 = cells(mp.astype(np.float32), grp, season, stratum, seasons, 12)
        obs2, _, p2 = boot(sm2, cn2, max(B // 2, 500), rng, keep=False)
        q2 = np.full(p2.shape, np.nan); rej = np.zeros(p2.shape, bool)
        for ic in range(len(cfgs)):
            for iv in range(8):
                pm = p2[ic, iv].copy(); pm[:, ~sect] = np.nan
                qq, rr = bh(pm)
                q2[ic, iv], rej[ic, iv] = qq, rr
        np.savez_compressed(os.path.join(OUT, f"maps_{b}_{cmpname}.npz"), diff=obs2[0].astype(np.float32), hf=obs2[1].astype(np.float32),
                            sf=obs2[2].astype(np.float32), p=p2.astype(np.float32), q=q2.astype(np.float32), rej=rej,
                            cfgs=np.array(cfgs), vars=np.array(VARS + ["z500_raw", "u250_raw", "mslp_raw"]), lat=LAT, lon=LON,
                            box=mbox, n_hf=len(th), n_sf=len(ts))
        for ic, cfg in enumerate(cfgs):
            for iv in range(5):
                summ.append(f"  map {b} {cmpname} {cfg} {VARS[iv]}: {int(rej[ic, iv].sum())} of {int(sect.sum() * len(LAT))} sector cells pass BH q<0.05")
    print(b, "done", flush=True)

BX = pd.DataFrame(rows_box)
# primary family: L and k0, 2 basins x 5 fields x 2 anchors comparisons? pre-registered: 20 tests = 2 basins x 5 fields x {L, k0}, comparison C1
prim = BX[(BX.cmp == "C1") & BX.cfg.isin(["L", "k0"])].copy()
prim["q_primary"] = bh(prim.p.values)[0]
sec = BX[(BX.cmp == "C1") & BX.cfg.isin(["k7", "k4", "k2", "k1"])].copy()
sec["q_secondary"] = bh(sec.p.values)[0]
c2 = BX[(BX.cmp == "C2") & BX.cfg.isin(["L", "k0"])].copy()
c2["q_c2"] = bh(c2.p.values)[0]
BX = BX.merge(prim[["basin", "cmp", "cfg", "var", "q_primary"]], how="left").merge(sec[["basin", "cmp", "cfg", "var", "q_secondary"]], how="left").merge(c2[["basin", "cmp", "cfg", "var", "q_c2"]], how="left")
BX.to_csv(os.path.join(OUT, "boxtests.csv"), index=False)
pd.DataFrame(rows_ser).to_csv(os.path.join(OUT, "series.csv"), index=False)
open(os.path.join(OUT, "summary.txt"), "w").write("\n".join(summ) + f"\nprimary family: {int((prim.q_primary < 0.05).sum())} of {len(prim)} pass q<0.05\nsecondary (k=7,4,2,1): {int((sec.q_secondary < 0.05).sum())} of {len(sec)}\nC2 (HF peak vs SF peak; L,k0): {int((c2.q_c2 < 0.05).sum())} of {len(c2)}\n")
print(open(os.path.join(OUT, "summary.txt")).read()[-400:])
