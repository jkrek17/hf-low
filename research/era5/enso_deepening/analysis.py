"""Pre-registered tests of PREREGISTRATION.md: ENSO flavor vs where Pacific storms deepen fastest.
ERA5 proxy; PR 38 tracker output (pipeline A detector/linker); HF label from pipeline A.

usage: analysis.py POINTS_CSV_GZ PIPELINE_A_MATCH_CSV WINTER_INDICES_CSV SST_RAW_CSV OUT_DIR [NPERM] [NBOOT]
Writes results.csv (every test), summary.txt, winter_table.csv, power.csv.
"""
import sys, os
import numpy as np, pandas as pd

P = pd.read_csv(sys.argv[1], dtype={"tm": str, "t_min": str})
M = pd.read_csv(sys.argv[2])
IDX = pd.read_csv(sys.argv[3]).set_index("winter")
SST = pd.read_csv(sys.argv[4])
OUT = sys.argv[5]; os.makedirs(OUT, exist_ok=True)
NPERM = int(sys.argv[6]) if len(sys.argv) > 6 else 10000
NBOOT = int(sys.argv[7]) if len(sys.argv) > 7 else 5000
RNG = np.random.default_rng(20261008)
WINT = np.arange(1979, 2026)

# ---- SON indices (same anomaly construction as enso_kuroshio/indices.py)
d = SST.copy()
for k in ["nino3", "nino4", "nino34", "emi_a", "emi_b", "emi_c"]:
    clim = d[(d.year >= 1991) & (d.year <= 2020)].groupby("month")[k].mean()
    d[k + "_a"] = d[k] - d.month.map(clim)
d["emi"] = d.emi_a_a - 0.5 * d.emi_b_a - 0.5 * d.emi_c_a
son = d[d.month.isin([9, 10, 11])].groupby("year")[["emi", "nino34_a"]].mean()
IDX["emi_son"] = son.emi.reindex(IDX.index).values
IDX["n34_son"] = son.nino34_a.reindex(IDX.index).values

def zs(v): v = np.asarray(v, float); return (v - v.mean()) / v.std(ddof=1)
Z = {c: pd.Series(zs(IDX.loc[WINT, c].values), index=WINT) for c in ["emi", "n34", "n4mn3", "pna", "emi_son", "n34_son"]}
Z["year"] = pd.Series(zs(WINT), index=WINT)

# ---- populations
P["win_ok"] = P.winter.notna()
def inwin(tm, s, win):
    if win == "DJF": lo, hi = f"{s}120300", f"{s+1}030100"
    else: lo, hi = f"{s+1}010100", f"{s+1}032900"
    return (tm >= lo) & (tm < hi)
BOX = (P.latm >= 25) & (P.latm <= 67) & (P.lonm >= 120) & (P.lonm <= 240)
BOXMIN = (P.lat_min >= 25) & (P.lat_min <= 67) & (P.lon_min >= 120) & (P.lon_min <= 240)
HFT = M[M.gust800_kt >= 71.7].tid.dropna().astype(int).unique()

def pop(kind, win="DJF"):
    m = BOX.copy()
    if kind == "bomb": m &= P.B >= 24
    elif kind == "rapid": m &= P.B >= 12
    elif kind == "hf": m &= P.tid.isin(HFT)
    ok = np.zeros(len(P), bool)
    for s in WINT:
        ok |= ((P.winter == s) & inwin(P.tm, s, win)).values
    return P[m.values & ok]

def wtab(df, col, how="mean", wcol=None):
    rows = []
    for s in WINT:
        x = df[df.winter == s]
        if len(x) == 0: rows.append((s, np.nan, 0)); continue
        v = x[col].values
        if how == "mean": y = v.mean()
        elif how == "p75": y = np.percentile(v, 75)
        elif how == "p25": y = np.percentile(v, 25)
        elif how == "wmean": y = np.average(v, weights=x[wcol].values)
        rows.append((s, y, len(x)))
    return pd.DataFrame(rows, columns=["winter", "y", "n"]).set_index("winter")

# ---- weighted least squares + Freedman-Lane permutation + winter pairs bootstrap
def wls(y, X, w):
    sw = np.sqrt(w); A = X * sw[:, None]
    b = np.linalg.lstsq(A, y * sw, rcond=None)[0]
    return b

def fit(tab, pred, controls, winters=None, nperm=None, nboot=None):
    nperm = NPERM if nperm is None else nperm; nboot = NBOOT if nboot is None else nboot
    w_ = np.array(WINT if winters is None else winters)
    tb = tab.loc[w_].dropna(subset=["y"]); w_ = tb.index.values
    cols = [np.ones(len(w_))] + [Z[c].loc[w_].values for c in [pred] + controls]
    X = np.column_stack(cols); y = tb.y.values; w = tb.n.values.astype(float)
    sw = np.sqrt(w); A = X * sw[:, None]; ys = y * sw
    G = np.linalg.pinv(A)                       # b = G @ ys
    b = G @ ys; res = ys - A @ b
    se_wald = np.sqrt(np.diag(np.linalg.inv(A.T @ A)) * (res @ res) / (len(y) - X.shape[1]))
    # Freedman-Lane: reduced model without the predictor
    Ar = np.delete(A, 1, axis=1); br = np.linalg.lstsq(Ar, ys, rcond=None)[0]
    fr = Ar @ br; rr = ys - fr
    perm = np.array([RNG.permutation(rr) for _ in range(nperm)])
    bp = (fr[None, :] + perm) @ G[1]
    p = (np.sum(np.abs(bp) >= abs(b[1]) - 1e-12) + 1) / (nperm + 1)
    # pairs bootstrap over winters
    bs = np.empty(nboot); k = len(y)
    for i in range(nboot):
        ix = RNG.integers(0, k, k)
        try: bs[i] = np.linalg.lstsq(A[ix], ys[ix], rcond=None)[0][1]
        except Exception: bs[i] = np.nan
    lo, hi = np.nanpercentile(bs, [2.5, 97.5]); se = np.nanstd(bs)
    return dict(n_winters=len(y), mean_storms=float(w.mean()), b=b[1], se_boot=se, se_wald=se_wald[1],
                lo=lo, hi=hi, detect80=2.8 * se, p_perm=p, sd_y=float(np.sqrt(np.average((y - np.average(y, weights=w)) ** 2, weights=w))))

R = []
def add(id_, desc, res, **kw): R.append(dict(id=id_, desc=desc, **res, **kw))

bombs = pop("bomb"); rapid = pop("rapid")
T_lon = wtab(bombs, "lonm"); T_lat = wtab(bombs, "latm")
T_lon.join(T_lat, rsuffix="_lat").to_csv(os.path.join(OUT, "winter_table.csv"))
C = ["n34", "year"]
add("P1", "bomb deepening longitude ~ EMI + N34 + year (47 winters)", fit(T_lon, "emi", C))
add("P2", "bomb deepening latitude ~ EMI + N34 + year", fit(T_lat, "emi", C))
add("S1", "rapid (B>=12) deepening longitude", fit(wtab(rapid, "lonm"), "emi", C))
add("S2", "rapid (B>=12) deepening latitude", fit(wtab(rapid, "latm"), "emi", C))
add("S3", "bomb longitude with SON EMI + SON N34 + year", fit(T_lon, "emi_son", ["n34_son", "year"]))
add("S4", "bomb latitude with SON EMI + SON N34 + year", fit(T_lat, "emi_son", ["n34_son", "year"]))
# S5 position of lowest pressure, same bomb tracks
bm = bombs[BOXMIN.loc[bombs.index].values]
add("S5", "longitude of minimum-pressure fix of bomb tracks (min fix in box)", fit(wtab(bm, "lon_min"), "emi", C))
add("S6", "bomb longitude with N4-N3 instead of EMI", fit(T_lon, "n4mn3", C))
add("S7", "bomb longitude without year term", fit(T_lon, "emi", ["n34"]))
add("S8", "bomb longitude without N34 term", fit(T_lon, "emi", ["year"]))
add("S9", "B-weighted centroid longitude of rapid storms", fit(wtab(rapid, "lonm", "wmean", "B"), "emi", C))
add("S10", "bomb longitude, JFM window", fit(wtab(pop("bomb", "JFM"), "lonm"), "emi", C))
hf = pop("hf")
add("S11", "longitude of the deepening point of HF-reaching storms (pipeline A, 22 winters 2004-05..), no year term",
    fit(wtab(hf, "lonm"), "emi", ["n34"], winters=range(2004, 2026)))
# S12 CP vs EP El Nino winters
en = IDX[IDX.n34 >= 0.5].copy(); en["cls"] = np.where(en.n4 > en.n3, "CP", "EP")
yv = T_lon.loc[en.index, "y"].values; lab = (en.cls.values == "CP")
d0 = yv[lab].mean() - yv[~lab].mean()
ge = 0
for _ in range(NPERM):
    pl = RNG.permutation(lab); ge += abs(yv[pl].mean() - yv[~pl].mean()) >= abs(d0) - 1e-12
se12 = np.sqrt(yv[lab].var(ddof=1) / lab.sum() + yv[~lab].var(ddof=1) / (~lab).sum())
R.append(dict(id="S12", desc=f"El Nino winters: CP (n={lab.sum()}) minus EP (n={(~lab).sum()}) mean bomb longitude (deg, not per SD)",
              n_winters=len(yv), mean_storms=np.nan, b=d0, se_boot=se12, se_wald=se12, lo=d0 - 1.96 * se12, hi=d0 + 1.96 * se12,
              detect80=2.8 * se12, p_perm=(ge + 1) / (NPERM + 1), sd_y=np.nan,
              extra=f"CP {sorted(en.index[lab].tolist())}; EP {sorted(en.index[~lab].tolist())}; CP mean {yv[lab].mean():.2f} EP mean {yv[~lab].mean():.2f}"))
add("S13", "winter 75th percentile of bomb longitude", fit(wtab(bombs, "lonm", "p75"), "emi", C))
add("S14", "winter 25th percentile of bomb longitude", fit(wtab(bombs, "lonm", "p25"), "emi", C))
add("S15", "bomb longitude, 1979-80..2003-04 (25 winters)", fit(T_lon, "emi", C, winters=range(1979, 2004)))
add("S16", "bomb longitude, 2004-05..2025-26 (22 winters)", fit(T_lon, "emi", C, winters=range(2004, 2026)))

D = pd.DataFrame(R)
pv = D.p_perm.values; o = np.argsort(pv); m = len(pv); q = np.empty(m); run = 1.0
for rank, i in reversed(list(enumerate(o, 1))):
    run = min(run, pv[i] * m / rank); q[i] = run
D["q_bh"] = q
# ---- mechanism M1 (outside the family)
m1 = fit(T_lon, "emi", ["n34", "year", "pna"], nperm=2000, nboot=2000)
m2 = fit(T_lon, "pna", ["emi", "n34", "year"], nperm=2000, nboot=2000)
mech = dict(M1_emi_with_pna=m1["b"], M1_lo=m1["lo"], M1_hi=m1["hi"], M1_p=m1["p_perm"], M2_pna_b=m2["b"], M2_lo=m2["lo"], M2_hi=m2["hi"],
            M2_p=m2["p_perm"], corr_emi_pna=float(np.corrcoef(Z["emi"], Z["pna"])[0, 1]), corr_emi_n34=float(np.corrcoef(Z["emi"], Z["n34"])[0, 1]),
            corr_emi_emison=float(np.corrcoef(Z["emi"], Z["emi_son"])[0, 1]))

# ---- planted-effect power (null surrogate: permute reduced-model residuals, add delta * EMI_z)
def power(tab, pred, controls, deltas, nplant=300, nperm=499):
    tb = tab.dropna(subset=["y"]); w_ = tb.index.values
    X = np.column_stack([np.ones(len(w_))] + [Z[c].loc[w_].values for c in [pred] + controls])
    w = tb.n.values.astype(float); sw = np.sqrt(w); A = X * sw[:, None]; G = np.linalg.pinv(A)
    ys = tb.y.values * sw
    Ar = np.delete(A, 1, axis=1); br = np.linalg.lstsq(Ar, ys, rcond=None)[0]; fr = Ar @ br; rr = ys - fr
    zc = Z[pred].loc[w_].values * sw
    out = []
    for dl in deltas:
        hit = 0
        for _ in range(nplant):
            yp = fr + RNG.permutation(rr) + dl * zc          # null surrogate + planted effect
            b = (G @ yp)[1]
            rp = yp - Ar @ np.linalg.lstsq(Ar, yp, rcond=None)[0]; fp = yp - rp
            perm = np.array([RNG.permutation(rp) for _ in range(nperm)])
            bp = (fp[None, :] + perm) @ G[1]
            hit += (np.sum(np.abs(bp) >= abs(b) - 1e-12) + 1) / (nperm + 1) < 0.05
        out.append(hit / nplant)
    return out
grid = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6]
pw = []
for name, tab, pred, ctl in (("P1", T_lon, "emi", C), ("P2", T_lat, "emi", C), ("S3", T_lon, "emi_son", ["n34_son", "year"])):
    pw.append(pd.DataFrame(dict(test=name, delta_deg_per_sd=grid, power=power(tab, pred, ctl, grid))))
PW = pd.concat(pw); PW.to_csv(os.path.join(OUT, "power.csv"), index=False)

D.to_csv(os.path.join(OUT, "results.csv"), index=False)
with open(os.path.join(OUT, "summary.txt"), "w") as fh:
    fh.write("ENSO flavor vs where Pacific storms deepen fastest. ERA5 proxy; PR 38 tracker output (pipeline A detector/linker); HF label from pipeline A.\n")
    fh.write(f"permutations {NPERM}; bootstrap {NBOOT}; BH over 18 tests. Effects are degrees per SD of the index (S12: degrees).\n\n")
    for _, r in D.iterrows():
        fh.write(f"{r.id:4s} n={int(r.n_winters):2d} b={r.b:+.3f} ({r.lo:+.3f},{r.hi:+.3f}) detect80={r.detect80:.2f} p={r.p_perm:.4f} q={r.q_bh:.3f}  {r.desc}\n")
        if isinstance(r.get('extra', None), str): fh.write(f"       {r.extra}\n")
    fh.write("\nbombs per winter: mean %.1f, min %d, max %d; rapid per winter mean %.1f; HF-reaching per winter (2004-25) mean %.1f\n" % (
        T_lon.n.mean(), T_lon.n.min(), T_lon.n.max(), wtab(rapid, "lonm").n.mean(), wtab(hf, "lonm").loc[2004:2025, "n"].mean()))
    fh.write("SD of winter-mean bomb longitude %.2f deg; latitude %.2f deg\n" % (D.loc[D.id == "P1", "sd_y"].iloc[0], D.loc[D.id == "P2", "sd_y"].iloc[0]))
    fh.write("\nMechanism: " + ", ".join(f"{k}={v:.3f}" for k, v in mech.items()) + "\n")
    fh.write("\nPlanted-effect power (share of 300 plants with permutation p<0.05)\n" + PW.pivot(index="delta_deg_per_sd", columns="test", values="power").to_string() + "\n")
print(open(os.path.join(OUT, "summary.txt")).read())
