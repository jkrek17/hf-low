"""Pre-registered tests (PREREGISTRATION.md): ENSO flavor vs Kuroshio genesis, ERA5 proxy.

usage: analysis.py TRACKS_CSV_GZ PIPELINE_A_MATCH_CSV WINTER_INDICES_CSV OUT_DIR [NPERM]
Writes results.csv (every test), summary.txt, genesis_by_winter.csv.
"""
import sys, os, json
import numpy as np, pandas as pd
from scipy import stats

NPERM = int(sys.argv[5]) if len(sys.argv) > 5 else 10000
RNG = np.random.default_rng(20261008)
T = pd.read_csv(sys.argv[1], dtype={"t0": str, "t_end": str})
M = pd.read_csv(sys.argv[2])
IDX = pd.read_csv(sys.argv[3]).set_index("winter")
OUT = sys.argv[4]; os.makedirs(OUT, exist_ok=True)
HFKT = 71.7

# ---- HF label from pipeline A (any matched Pacific pipeline A track with gust index >= 71.7 kt)
hf = M[M.gust800_kt >= HFKT].tid.dropna().astype(int).unique()
T["hf"] = T.tid.isin(hf)

# ---- windows (genesis time strings YYYYMMDDHH)
def inwin(row_t0, s, name):
    if name == "DJF":
        lo, hi = f"{s}120300", f"{s+1}030100"
    elif name == "JFM":
        lo, hi = f"{s+1}010100", f"{s+1}032900"
    else:
        lo, hi = f"{s}120300", f"{s+1}032900"
    return (row_t0 >= lo) & (row_t0 < hi)
def days(s, name):
    leap = (s + 1) % 4 == 0 and ((s + 1) % 100 != 0 or (s + 1) % 400 == 0)
    return {"DJF": 88, "JFM": 87, "DM": 116}[name] + (1 if leap else 0)
WIN = {"DJF": "DJF", "JFM": "JFM", "DM": "DM"}

# ---- boxes
def box(df, name):
    la, lo = df.lat0, df.lon0
    def rect(a0, a1, o0, o1): return (la >= a0) & (la <= a1) & (lo >= o0) & (lo <= o1)
    if name == "K": m = rect(25, 35, 120, 142) | rect(30, 42, 142, 165)
    elif name == "K_N": m = rect(30, 40, 120, 142) | rect(35, 47, 142, 165)
    elif name == "K_E": m = rect(25, 35, 130, 152) | rect(30, 42, 152, 175)
    elif name == "K_ECS": m = rect(25, 35, 120, 142)
    elif name == "K_KE": m = rect(30, 42, 142, 165)
    return m & df.ocean0
BASIN = (T.lat0 >= 25) & (T.lat0 <= 67) & (T.lon0 >= 120) & (T.lon0 <= 240) & T.ocean0

def table(boxname, win, hfonly=False):
    rows = []
    for s in range(1979, 2026):
        t = T[(T.winter == s) & inwin(T.t0, s, win)]
        b = box(t, boxname); bas = BASIN.loc[t.index]
        rows.append(dict(winter=s, days=days(s, win), n=int(b.sum()), basin=int(bas.sum()),
                         hf=int((b & t.hf).sum())))
    return pd.DataFrame(rows).set_index("winter")

# ---- GLM (quasi-Poisson / quasi-binomial) by IRLS
def irls(y, X, off=None, n=None, fam="pois", it=60):
    off = np.zeros(len(y)) if off is None else off
    if fam == "pois":
        mu = y + 0.5; eta = np.log(mu)
    else:
        mu = (y + 0.5) / (n + 1); eta = np.log(mu / (1 - mu))
    b = np.zeros(X.shape[1])
    for _ in range(it):
        if fam == "pois":
            w = mu; z = eta - off + (y - mu) / mu
        else:
            w = n * mu * (1 - mu); z = eta + (y / n - mu) / (mu * (1 - mu))
        XtW = X.T * w
        bn = np.linalg.solve(XtW @ X + 1e-10 * np.eye(X.shape[1]), XtW @ z)
        eta = X @ bn + (off if fam == "pois" else 0)
        mu = np.exp(eta) if fam == "pois" else 1 / (1 + np.exp(-eta))
        if np.max(np.abs(bn - b)) < 1e-9:
            b = bn; break
        b = bn
    if fam == "pois":
        w = mu; pr = (y - mu) / np.sqrt(mu)
    else:
        w = n * mu * (1 - mu); pr = (y - n * mu) / np.sqrt(n * mu * (1 - mu))
    phi = max((pr ** 2).sum() / (len(y) - X.shape[1]), 1.0)
    cov = phi * np.linalg.inv((X.T * w) @ X)
    return b, np.sqrt(np.diag(cov)), phi, mu

def zs(v): v = np.asarray(v, float); return (v - v.mean()) / v.std(ddof=1)

def design(winters, pred, controls, ind=IDX):
    cols = [np.ones(len(winters))]
    names = ["const"]
    for nm in [pred] + controls:
        v = np.asarray(winters, float) if nm == "year" else (ind.loc[winters, nm].values if nm != "era" else (np.asarray(winters) >= 2001).astype(float))
        cols.append(zs(v) if nm != "era" else v); names.append(nm)
    return np.column_stack(cols), names

def fit_test(tab, y_col, pred, controls, fam="pois", winters=None, perm=True, label=""):
    w = np.array(tab.index if winters is None else winters)
    tb = tab.loc[w]
    X, names = design(w, pred, controls)
    if fam == "pois":
        y = tb[y_col].values.astype(float); off = np.log(tb.days.values.astype(float)); n = None
    else:
        y = tb[y_col].values.astype(float); n = tb["basin" if y_col == "n" else "n"].values.astype(float); off = None
    b, se, phi, mu = irls(y, X, off, n, fam)
    out = dict(label=label, n_winters=len(w), b=b[1], se=se[1], rr=np.exp(b[1]),
               lo=np.exp(b[1] - 1.96 * se[1]), hi=np.exp(b[1] + 1.96 * se[1]), detect80=np.exp(2.8 * se[1]),
               phi=phi, mean_y=float(np.mean(y)) if fam == "pois" else float(np.mean(y / n)))
    if perm:
        # Freedman-Lane (counts), as registered; predictor-residual permutation for binomial outcomes
        Xr = np.delete(X, 1, axis=1)
        if fam == "pois":
            br, _, _, mur = irls(y, Xr, off, None, "pois")
            r = (y - mur) / np.sqrt(mur)
            cnt = 0; ge = 0
            for _ in range(NPERM):
                ys = np.maximum(mur + RNG.permutation(r) * np.sqrt(mur), 1e-6)
                bs = irls(ys, X, off, None, "pois", it=30)[0][1]
                ge += abs(bs) >= abs(b[1]) - 1e-12
            out["p_perm"] = (ge + 1) / (NPERM + 1)
        else:
            g = np.linalg.lstsq(Xr, X[:, 1], rcond=None)[0]; fit = Xr @ g; e = X[:, 1] - fit
            ge = 0
            for _ in range(NPERM):
                Xs = X.copy(); Xs[:, 1] = fit + RNG.permutation(e)
                bs = irls(y, Xs, off, n, fam, it=30)[0][1]
                ge += abs(bs) >= abs(b[1]) - 1e-12
            out["p_perm"] = (ge + 1) / (NPERM + 1)
        out["p_wald"] = 2 * stats.norm.sf(abs(b[1] / se[1]))
    out["coefs"] = dict(zip(names, b))
    return out

def predperm_p(tab, y_col, pred, controls, winters=None):   # predictor-residual p for the counts, as a cross-check
    w = np.array(tab.index if winters is None else winters); tb = tab.loc[w]
    X, _ = design(w, pred, controls); y = tb[y_col].values.astype(float); off = np.log(tb.days.values.astype(float))
    b = irls(y, X, off)[0][1]
    Xr = np.delete(X, 1, axis=1); g = np.linalg.lstsq(Xr, X[:, 1], rcond=None)[0]; fit = Xr @ g; e = X[:, 1] - fit
    ge = 0
    for _ in range(NPERM):
        Xs = X.copy(); Xs[:, 1] = fit + RNG.permutation(e)
        ge += abs(irls(y, Xs, off, None, "pois", it=30)[0][1]) >= abs(b) - 1e-12
    return (ge + 1) / (NPERM + 1)

R = []
def add(id_, desc, res, **kw):
    r = dict(id=id_, desc=desc, **{k: v for k, v in res.items() if k != "coefs"}, **kw); R.append(r); return r

CTRL = ["n34", "year"]
tabK = table("K", "DJF"); tabK.to_csv(os.path.join(OUT, "genesis_by_winter.csv"))
P = add("P1", "DJF K genesis count ~ EMI + N34 + year (n=47)", fit_test(tabK, "n", "emi", CTRL, label="P1"))
fit_P1_coefs = fit_test(tabK, "n", "emi", CTRL, perm=False)["coefs"]
add("S1", "DJF share of basin genesis in K ~ EMI + N34 + year", fit_test(tabK, "n", "emi", CTRL, fam="binom", label="S1"))
w22 = list(range(2004, 2026))
add("S2", "DJF K-genesis HF-reaching count ~ EMI + N34 (22 winters 2004-05..)", fit_test(tabK, "hf", "emi", ["n34"], winters=w22, label="S2"))
# S3 binomial with trials = K genesis count
def fit_share_hf():
    tb = tabK.loc[w22]; X, names = design(w22, "emi", ["n34"])
    y = tb.hf.values.astype(float); n = tb.n.values.astype(float)
    b, se, phi, mu = irls(y, X, None, n, "binom")
    Xr = np.delete(X, 1, axis=1); g = np.linalg.lstsq(Xr, X[:, 1], rcond=None)[0]; fit = Xr @ g; e = X[:, 1] - fit
    ge = sum(abs(irls(y, np.column_stack([X[:, 0], fit + RNG.permutation(e), X[:, 2]]), None, n, "binom", it=30)[0][1]) >= abs(b[1]) - 1e-12 for _ in range(NPERM))
    return dict(label="S3", n_winters=22, b=b[1], se=se[1], rr=np.exp(b[1]), lo=np.exp(b[1] - 1.96 * se[1]), hi=np.exp(b[1] + 1.96 * se[1]),
                detect80=np.exp(2.8 * se[1]), phi=phi, mean_y=float((y.sum()) / n.sum()), p_perm=(ge + 1) / (NPERM + 1),
                p_wald=2 * stats.norm.sf(abs(b[1] / se[1])), coefs={})
add("S3", "HF share among K-genesis storms ~ EMI + N34 (22 winters; effect = odds ratio per SD)", fit_share_hf())
add("S4", "P1 with window JFM", fit_test(table("K", "JFM"), "n", "emi", CTRL, label="S4"))
add("S5", "P1 with window Dec-Mar", fit_test(table("K", "DM"), "n", "emi", CTRL, label="S5"))
add("S6", "P1 with N4-N3 in place of EMI", fit_test(tabK, "n", "n4mn3", CTRL, label="S6"))
for sid, bn in (("S7", "K_N"), ("S8", "K_E"), ("S9", "K_ECS"), ("S10", "K_KE")):
    add(sid, f"P1 with box {bn}", fit_test(table(bn, "DJF"), "n", "emi", CTRL, label=sid))
add("S11", "HF-reaching K-genesis count 1979-2025 ~ EMI + N34 + era dummy (>=2001)", fit_test(tabK, "hf", "emi", ["n34", "era"], label="S11"))
# S12 CP vs EP El Nino winters (permutation of class labels)
en = IDX[IDX.n34 >= 0.5].copy(); en["cls"] = np.where(en.n4 > en.n3, "CP", "EP")
rate = (tabK.n / tabK.days)[en.index].values * 88
cp = rate[en.cls.values == "CP"]; ep = rate[en.cls.values == "EP"]
d0 = cp.mean() - ep.mean(); allv = rate.copy(); lab = (en.cls.values == "CP")
ge = sum(abs(allv[(pl := RNG.permutation(lab))].mean() - allv[~pl].mean()) >= abs(d0) - 1e-12 for _ in range(NPERM))
R.append(dict(id="S12", desc=f"El Nino winters (N34>=0.5): CP (n={len(cp)}) minus EP (n={len(ep)}) mean DJF K genesis per 88 days",
              n_winters=len(rate), b=d0, se=np.sqrt(cp.var(ddof=1) / len(cp) + ep.var(ddof=1) / len(ep)), rr=cp.mean() / ep.mean(),
              lo=np.nan, hi=np.nan, detect80=np.nan, phi=np.nan, mean_y=rate.mean(), p_perm=(ge + 1) / (NPERM + 1), p_wald=np.nan,
              extra=f"CP {sorted(en.index[en.cls=='CP'].tolist())}; EP {sorted(en.index[en.cls=='EP'].tolist())}; CP mean {cp.mean():.2f} EP mean {ep.mean():.2f}"))
add("S13", "P1 without the year term", fit_test(tabK, "n", "emi", ["n34"], label="S13"))
add("S14", "P1 without the N34 term", fit_test(tabK, "n", "emi", ["year"], label="S14"))

# ---- BH over the 15 tests
D = pd.DataFrame(R)
pv = D.p_perm.values; o = np.argsort(pv); m = len(pv)
q = np.empty(m); run = 1.0
for rank, i in reversed(list(enumerate(o, 1))):
    run = min(run, pv[i] * m / rank); q[i] = run
D["q_bh"] = q
# cross-check of P1 with the predictor-residual permutation and a Wald interval
p1_alt = predperm_p(tabK, "n", "emi", CTRL)
# ---- mechanism (outside the FDR family)
mech = {}
tabP = tabK.copy()
with_pna = fit_test(tabP, "n", "emi", ["n34", "year", "pna"], label="M1")
mech["M1_b_emi_with_pna"] = with_pna["b"]; mech["M1_rr_with_pna"] = with_pna["rr"]; mech["M1_rr_without_pna"] = P["rr"]
mech["M1_ci_with_pna"] = [with_pna["lo"], with_pna["hi"]]
mech["M2_pna_coef"] = with_pna["coefs"]["pna"]
Xp, _ = design(list(tabP.index), "emi", ["n34", "year", "pna"])
pn = fit_test(tabP, "n", "pna", ["emi", "n34", "year"], label="M2")
mech["M2_rr_pna"] = pn["rr"]; mech["M2_ci"] = [pn["lo"], pn["hi"]]; mech["M2_p_perm"] = pn["p_perm"]
mech["corr_emi_pna"] = float(IDX.emi.corr(IDX.pna)); mech["corr_emi_n34"] = float(IDX.emi.corr(IDX.n34))
mech["corr_n34_official_oni"] = float(IDX.n34.corr(IDX.oni_official))
mech["P1_crosscheck_p_predictor_residual_perm"] = p1_alt
mech["vif_emi"] = float(1 / (1 - np.corrcoef(zs(IDX.emi), zs(IDX.n34))[0, 1] ** 2))
# descriptives
desc = dict(mean_K_DJF=float(tabK.n.mean()), sd_K_DJF=float(tabK.n.std()), var_to_mean=float(tabK.n.var() / tabK.n.mean()),
            mean_basin_DJF=float(tabK.basin.mean()), mean_hf_K_DJF_all=float(tabK.hf.mean()), mean_hf_K_DJF_2004on=float(tabK.loc[w22].hf.mean()),
            mean_K_DJF_2004on=float(tabK.loc[w22].n.mean()))
D.to_csv(os.path.join(OUT, "results.csv"), index=False)
json.dump(dict(mechanism=mech, descriptives=desc, P1_coefs=fit_P1_coefs, nperm=NPERM), open(os.path.join(OUT, "mechanism_descriptives.json"), "w"), indent=1)
with open(os.path.join(OUT, "summary.txt"), "w") as f:
    f.write("ENSO flavor vs Kuroshio genesis. ERA5 proxy; new MSLP tracker on pipeline A's detector; HF labels from pipeline A.\n")
    f.write(f"permutations: {NPERM}; BH over 15 tests (P1, S1-S14).\n\n")
    for _, r in D.iterrows():
        f.write(f"{r.id:4s} n={int(r.n_winters):2d} effect={r.rr:.3f} ({r.lo:.3f}-{r.hi:.3f}) detect80={r.detect80:.3f} p_perm={r.p_perm:.4f} q_BH={r.q_bh:.3f}  {r.desc}\n")
        if isinstance(r.get('extra'), str): f.write(f"      {r.extra}\n")
    f.write("\nmechanism:\n" + json.dumps(mech, indent=1) + "\n\ndescriptives:\n" + json.dumps(desc, indent=1) + "\n")
print(open(os.path.join(OUT, "summary.txt")).read())
