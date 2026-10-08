"""RA-20 analysis (Arms A and B), exactly as PREREGISTRATION.md. ERA5 proxy, pipeline A.

usage: python3 -I analyse.py <repo_root> --power      power by group (needs features, uses case counts only, no outcome-feature join)
       python3 -I analyse.py <repo_root> --run [B|A|BA]   all tests; writes results/analysis_*.txt/csv
Feature CSVs are read from $ERA5_WORK/explosive_onset/{armB,armA}/ (made by extract_armB.py / extract_armA.py).
"""
import sys, os, glob, json, numpy as np, pandas as pd
from scipy import stats
root = sys.argv[1]; MODE = sys.argv[2]; ARMS = sys.argv[3] if len(sys.argv) > 3 else "BA"
HERE = os.path.join(root, "research/era5/explosive_onset"); RES = os.path.join(HERE, "results")
WORK = os.environ.get("ERA5_WORK", os.path.join(HERE, "work"))
rng = np.random.RandomState(20261008)
GROUPS = {"B": {"B1 gust concentration": ["A45", "PK", "RG"], "B2 pressure structure": ["SHARP", "CURV", "GRAD"],
                "B3 wind vs gust": ["WS99", "GF", "A34"], "B4 shear and convergence": ["SHR", "CONV"]},
          "A": {"A_stab": ["stab"], "A_trough": ["trough"], "A_omega700": ["omega700"]}}
ALLF = {a: sum(GROUPS[a].values(), []) for a in GROUPS}
expit = lambda x: 1 / (1 + np.exp(-x))

# ---------------------------------------------------------------- data
def load(arm):
    f = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"))
    p = pd.read_csv(os.path.join(root, "research/era5/hf_boosted/results/loso_probs.csv.gz"), usecols=["track", "time", "M0", "M1"])
    f = f.merge(p, on=["track", "time"])
    S = f[(f.g800 < 55) & (~f.hf_now)].copy()
    T = pd.read_csv(os.path.join(RES, "times_armB.csv"))[["time", "kind", "weight"]]
    d = S.merge(T, on="time")
    ft = pd.concat([pd.read_csv(x) for x in glob.glob(os.path.join(WORK, "armB" if arm == "B" else "armA", "*.csv"))])
    ft = ft.drop(columns=[c for c in ["g800"] if c in ft.columns and arm == "B"], errors="ignore")
    d = d.merge(ft.drop(columns=[c for c in ["match_km", "g800_re"] if c in ft.columns and arm == "A"]), on=["track", "time"])
    d["case"] = d.hf24.astype(float); d["grp"] = pd.factorize(d.time)[0]
    d["lb1"] = np.log(d.M1.clip(1e-5, 1 - 1e-5) / (1 - d.M1.clip(1e-5, 1 - 1e-5)))
    d["lb0"] = np.log(d.M0.clip(1e-5, 1 - 1e-5) / (1 - d.M0.clip(1e-5, 1 - 1e-5)))
    return d.reset_index(drop=True)

def standardise(Xtr, wtr, Xte_list):
    mu = np.array([np.nansum(wtr * Xtr[:, k]) / wtr[~np.isnan(Xtr[:, k])].sum() for k in range(Xtr.shape[1])])
    sd = np.array([np.sqrt(np.nansum(wtr * (Xtr[:, k] - mu[k]) ** 2) / wtr[~np.isnan(Xtr[:, k])].sum()) for k in range(Xtr.shape[1])])
    sd[sd == 0] = 1
    return [np.nan_to_num((X - mu) / sd, nan=0.0) for X in Xte_list]

# ---------------------------------------------------------------- conditional (Breslow-type) logit on risk sets
def cl_parts(off, Z, y, g, ng, b):
    q = np.exp(off + Z @ b); S0 = np.bincount(g, q, ng); k = Z.shape[1]
    S1 = np.stack([np.bincount(g, q * Z[:, j], ng) for j in range(k)], 1)
    S2 = np.zeros((ng, k, k))
    for i in range(k):
        for j in range(i, k):
            S2[:, i, j] = S2[:, j, i] = np.bincount(g, q * Z[:, i] * Z[:, j], ng)
    nc = np.bincount(g, y, ng); ok = (S0 > 0) & (nc > 0)
    m = S1[ok] / S0[ok, None]
    U = (y[:, None] * Z).sum(0) - (nc[ok, None] * m).sum(0)
    I = (nc[ok, None, None] * (S2[ok] / S0[ok, None, None] - m[:, :, None] * m[:, None, :])).sum(0)
    ll = (y * (off + Z @ b)).sum() - (nc[ok] * np.log(S0[ok])).sum()
    return ll, U, I

def cl_fit(off, Z, y, g, ridge=0.0, iters=25):
    ng = g.max() + 1; b = np.zeros(Z.shape[1])
    for _ in range(iters):
        ll, U, I = cl_parts(off, Z, y, g, ng, b)
        step = np.linalg.solve(I + ridge * np.eye(len(b)), U - ridge * b)
        b = b + step
        if np.abs(step).max() < 1e-7: break
    ll, U, I = cl_parts(off, Z, y, g, ng, b)
    return b, ll, I

def cl_ll(off, Z, y, g, b):
    return cl_parts(off, Z, y, g, g.max() + 1, b)[0]

def relabel(g):
    return pd.factorize(g)[0]

def lr_test(off, Z, y, g, season, nboot):
    g = relabel(g); b, ll, I = cl_fit(off, Z, y, g); ll0 = cl_ll(off, Z, y, g, np.zeros(Z.shape[1]))
    lr = 2 * (ll - ll0); k = Z.shape[1]; p_lr = stats.chi2.sf(max(lr, 0), k)
    seas = np.unique(season); idx = {s: np.where(season == s)[0] for s in seas}
    bs = []
    for _ in range(nboot):
        pick = rng.choice(seas, len(seas), replace=True)
        rows = np.concatenate([idx[s] for s in pick]); gg = np.concatenate([g[idx[s]] + c * (g.max() + 1) for c, s in enumerate(pick)])
        try: bs.append(cl_fit(off[rows], Z[rows], y[rows], relabel(gg), iters=12)[0])
        except Exception: pass
    bs = np.array(bs); C = np.atleast_2d(np.cov(bs.T))
    wald = float(b @ np.linalg.pinv(C) @ b); p_bt = stats.chi2.sf(wald, k)
    return dict(k=k, b=b, se_boot=np.sqrt(np.diag(C)), lr=lr, p_lr=p_lr, p_boot=p_bt, p=max(p_lr, p_bt))

# ---------------------------------------------------------------- held-out within-time gain and BSS gain (LOSO)
def loso_models(d, cols, off_col, tree=False):
    off = d[off_col].to_numpy(); X = d[cols].to_numpy(float); y = d.case.to_numpy(); w = d.weight.to_numpy(); g = d.grp.to_numpy(); se = d.season.to_numpy()
    inrisk = (d.kind == "case_time").to_numpy()
    seasons = np.unique(se); ll_gain, ncs = [], []; qa = np.zeros(len(d)); qt = np.zeros(len(d))
    for s in seasons:
        tr, te = se != s, se == s
        Ztr, Zte = standardise(X[tr], w[tr], [X[tr], X[te]])
        # (1) within-time: case-time risk sets, ridge 1
        mtr, mte = tr & inrisk, te & inrisk
        Zm_tr, Zm_te = standardise(X[mtr], np.ones(mtr.sum()), [X[mtr], X[mte]])
        b, _, _ = cl_fit(off[mtr], Zm_tr, y[mtr], relabel(g[mtr]), ridge=1.0)
        gte = relabel(g[mte]); yte = y[mte]
        ll_gain.append(cl_ll(off[mte], Zm_te, yte, gte, b) - cl_ll(off[mte], Zm_te, yte, gte, np.zeros(len(b)))); ncs.append(yte.sum())
        # (2) weighted offset logistic with intercept, ridge 1 on slopes
        Xd = np.c_[np.ones(tr.sum()), Ztr]; bb = np.zeros(Xd.shape[1]); R = np.diag([0.0] + [1.0] * (Xd.shape[1] - 1))
        for _ in range(25):
            q = expit(off[tr] + Xd @ bb); H = (Xd * (w[tr] * q * (1 - q))[:, None]).T @ Xd + R
            st = np.linalg.solve(H, Xd.T @ (w[tr] * (y[tr] - q)) - R @ bb); bb += st
            if np.abs(st).max() < 1e-8: break
        qa[te] = expit(off[te] + np.c_[np.ones(te.sum()), Zte] @ bb)
        if tree:
            from sklearn.ensemble import HistGradientBoostingClassifier
            m = HistGradientBoostingClassifier(max_iter=100, learning_rate=0.06, max_leaf_nodes=8, l2_regularization=1.0, early_stopping=False, random_state=0)
            m.fit(np.c_[off[tr], Ztr], y[tr], sample_weight=w[tr]); qt[te] = m.predict_proba(np.c_[off[te], Zte])[:, 1]
    pb = expit(off)
    def bs(q, mask): return (w * (q - y) ** 2)[mask].sum()
    pi = (w * y).sum() / w.sum()
    res = dict(ll_gain_season=np.array(ll_gain), ncase=np.array(ncs))
    sb = np.array([bs(pb, se == s) for s in seasons]); sa = np.array([bs(qa, se == s) for s in seasons]); st_ = np.array([bs(qt, se == s) for s in seasons]) if tree else None
    clim = pi * (1 - pi) * w.sum()
    res.update(season_brier_base=sb, season_brier_aug=sa, season_brier_tree=st_, clim=clim, seasons=seasons)
    return res

def signflip(x, n=100000):
    f = rng.choice([-1.0, 1.0], size=(n, len(x))); return float((np.abs((f * x).sum(1)) >= abs(x.sum())).mean())

def boot_ci(sb, sa, clim, seas_w, n=2000):
    k = len(sb); out = []
    for _ in range(n):
        i = rng.randint(0, k, k); out.append((sb[i].sum() - sa[i].sum()) / (clim * (np.sum(seas_w[i]) / np.sum(seas_w))))
    return np.percentile(out, [5, 95])

def bh(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); q = np.empty(m); r = p[o] * m / (np.arange(m) + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1]; return np.minimum(q, 1)

# ---------------------------------------------------------------- power by group (no outcome joined to features)
def power(arm):
    d = load(arm); cols = ALLF[arm]
    ct = d[d.kind == "case_time"]; out = [f"Arm {arm}: case-time risk sets {ct.grp.nunique()}, cases {int(ct.case.sum())}, controls {int((1-ct.case).sum())}, fixes with features {len(ct)}"]
    for gname, gc in {**GROUPS[arm], f"{arm}_all": cols}.items():
        X = ct[gc].to_numpy(float); Z, = standardise(X, np.ones(len(X)), [X]); y = ct.case.to_numpy(); g = relabel(ct.grp.to_numpy())
        _, _, I = cl_parts(ct.lb1.to_numpy(), Z, y, g, g.max() + 1, np.zeros(len(gc)))
        se = np.sqrt(np.diag(np.linalg.inv(I))); mde = (stats.norm.isf(0.05 / 2 / (5 if arm == "B" else 4)) + 0.8416) * se; mde05 = (1.96 + 0.8416) * se
        out.append(f"  {gname}: features {gc}; 80% MDE per feature, log-odds per SD (adjusted for the others in the group): " + ", ".join(f"{c} {m:.3f} (alpha 0.05) / {m2:.3f} (corrected)" for c, m, m2 in zip(gc, mde05, mde)))
    open(os.path.join(RES, f"power_groups_{arm}.txt"), "w").write("\n".join(out) + "\n"); print("\n".join(out))

# ---------------------------------------------------------------- the run
def run(arm):
    d = load(arm); cols = ALLF[arm]; L = []; rows = []
    # reproduction check (arm B only): re-detected g800 vs the table where the table has an in-domain value
    if arm == "B":
        ok = d.g800 > 0; diff = (d.g800_re - d.g800)[ok].abs()
        L.append(f"g800 reproduction (fixes with table g800 > 0, n={int(ok.sum())}): {100*(diff<=0.5).mean():.2f}% within 0.5 kt, max difference {diff.max():.2f} kt; fixes with table g800 = 0 (outside pipeline A's index domain): {int((~ok).sum())}")
        if (diff <= 0.5).mean() < 0.99: L.append("REPRODUCTION CHECK FAILED (<99%): STOP"); open(os.path.join(RES, f"analysis_{arm}.txt"), "w").write("\n".join(L)); print("\n".join(L)); return
        L.append(f"match distance to re-detected low: {100*(d.match_km<=50).mean():.1f}% within 50 km; fixes with no owned ocean cell: {int(d[cols[0]].isna().sum())}")
    d = d[d[cols].notna().any(axis=1) | True].reset_index(drop=True)
    if arm == "A":
        L.append("Arm A restricted to times up to 2023-01-09 (WeatherBench2); baseline scored on the same fixes.")
    ct = d[d.kind == "case_time"]
    L.append(f"Arm {arm}: sampled fixes {len(d)}; case-time risk sets {ct.grp.nunique()}, cases {int(ct.case.sum())}, controls {int((1-ct.case).sum())}; sampled non-case-time fixes {int((d.kind!='case_time').sum())}; seasons with cases {ct[ct.case==1].season.nunique()}")
    tests = []
    glist = {**GROUPS[arm], f"{arm}_all": cols}
    for base, off_col in (("M1", "lb1"), ("M0", "lb0")):
        for gname, gc in glist.items():
            X = ct[gc].to_numpy(float); Z, = standardise(X, np.ones(len(X)), [X])
            r = lr_test(ct[off_col].to_numpy(), Z, ct.case.to_numpy(), ct.grp.to_numpy(), ct.season.to_numpy(), 2000 if base == "M1" else 500)
            m = loso_models(d, gc, off_col, tree=(base == "M1" and gname.endswith("_all")))
            gain = m["ll_gain_season"]; p_sf = signflip(gain)
            dbss = (m["season_brier_base"].sum() - m["season_brier_aug"].sum()) / m["clim"]
            sw = np.array([d.weight[d.season == s].sum() for s in m["seasons"]]); lo, hi = boot_ci(m["season_brier_base"], m["season_brier_aug"], m["clim"], sw)
            p_bss = signflip(m["season_brier_base"] - m["season_brier_aug"])
            row = dict(arm=arm, baseline=base, group=gname, k=len(gc), lr=r["lr"], p_lr=r["p_lr"], p_boot=r["p_boot"], p_decision=r["p"], coef_per_sd=";".join(f"{c}:{x:+.3f}±{s:.3f}" for c, x, s in zip(gc, r["b"], r["se_boot"])),
                       heldout_ll_gain_per_case=gain.sum() / m["ncase"].sum(), heldout_signflip_p=p_sf, dBSS_S=dbss, dBSS_lo5=lo, dBSS_hi95=hi, dBSS_signflip_p=p_bss)
            if m["season_brier_tree"] is not None:
                row["tree_dBSS_S"] = (m["season_brier_base"].sum() - m["season_brier_tree"].sum()) / m["clim"]
            rows.append(row); print(row["baseline"], gname, f"p={r['p']:.4g} dBSS={dbss:+.4f}", flush=True)
    R = pd.DataFrame(rows)
    prim = R.baseline == "M1"
    R.loc[prim, "q_within_arm"] = bh(R.loc[prim, "p_decision"].to_numpy())
    sec = R.baseline == "M0"
    R.loc[sec, "q_secondary"] = bh(R.loc[sec, "p_decision"].to_numpy())
    # subsets (descriptive): basin, dp12 present/missing, matched cases; 'all' model vs M1
    sub = []
    allc = cols
    for name, mask in (("Atlantic", ct.basin == "atl"), ("Pacific", ct.basin == "pac"), ("dp12 available", ct.dp12.notna()), ("dp12 missing (new storms)", ct.dp12.isna())):
        c2 = ct[mask.to_numpy()]
        if c2.case.sum() < 30: continue
        Z, = standardise(c2[allc].to_numpy(float), np.ones(len(c2)), [c2[allc].to_numpy(float)])
        r = lr_test(c2.lb1.to_numpy(), Z, c2.case.to_numpy(), c2.grp.to_numpy(), c2.season.to_numpy(), 300)
        sub.append(dict(arm=arm, subset=name, cases=int(c2.case.sum()), lr=r["lr"], p_decision=r["p"]))
    # matched sensitivity: same time, |g800 diff| <= 3, |dp12 diff| <= 2 (or both missing)
    keep = np.zeros(len(ct), bool); ctr = ct.reset_index(drop=True)
    for t, gg in ctr.groupby("time"):
        cs, ks = gg[gg.case == 1], gg[gg.case == 0]
        for i, r_ in cs.iterrows():
            okk = (abs(ks.g800 - r_.g800) <= 3) & ((abs(ks.dp12 - r_.dp12) <= 2) if not pd.isna(r_.dp12) else ks.dp12.isna())
            if okk.any(): keep[i] = True; keep[ks.index[okk.to_numpy()]] = True
    cm = ctr[keep]
    if cm.case.sum() >= 30:
        Z, = standardise(cm[allc].to_numpy(float), np.ones(len(cm)), [cm[allc].to_numpy(float)])
        r = lr_test(cm.lb1.to_numpy(), Z, cm.case.to_numpy(), cm.grp.to_numpy(), cm.season.to_numpy(), 300)
        sub.append(dict(arm=arm, subset="matched controls only", cases=int(cm.case.sum()), lr=r["lr"], p_decision=r["p"]))
    R.to_csv(os.path.join(RES, f"analysis_{arm}_tests.csv"), index=False); pd.DataFrame(sub).to_csv(os.path.join(RES, f"analysis_{arm}_subsets.csv"), index=False)
    L.append(R.drop(columns=["coef_per_sd"]).to_string(float_format=lambda x: f"{x:.4g}")); L.append(pd.DataFrame(sub).to_string(float_format=lambda x: f"{x:.4g}"))
    open(os.path.join(RES, f"analysis_{arm}.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))

if __name__ == "__main__":
    for a in [c for c in ARMS if c in "AB"]:
        power(a) if MODE == "--power" else run(a)
    if MODE == "--run" and len(ARMS) == 2:
        R = pd.concat([pd.read_csv(os.path.join(RES, f"analysis_{a}_tests.csv")) for a in "BA"]); pm = R.baseline == "M1"
        R.loc[pm, "q_all_primary"] = bh(R.loc[pm, "p_decision"].to_numpy()); R.to_csv(os.path.join(RES, "analysis_all_tests.csv"), index=False)
