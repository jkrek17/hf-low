"""The pre-registered RA-24 test: upstream trough depth against bombs, 1979-80..2003-04 (ERA5 proxy, pipeline A).

Runs the analysis declared in PREREGISTRATION.md, once. Both trough definitions (D1 climatological anomaly,
D2 zonal eddy) are always run and reported side by side. Scores the new seasons once; 2015-26 is never read.

usage: analyse.py WORK_DIR     -> results/*.csv, results/trough_bombs-result.txt, results/or_forest.png
"""
import sys, os, datetime
import numpy as np, pandas as pd
import tblib as T

work = sys.argv[1]
os.makedirs("results", exist_ok=True)
N, R = T.load_tables(work)
std = T.Std(R)
OUT = []


def w(s=""):
    OUT.append(s)
    print(s, flush=True)


NB = int(os.environ.get("NBOOT", 2000))
rowsP, rowsS, rowsX = [], [], []


def one(S, col, label, family, outcome="bomb", **kw):
    X = T.design(S, std, col, **kw)
    y = S[outcome].values
    r = T.summarize(S, y, X)
    r.update(label=label, family=family, definition=[d for d, c in T.DEFS.items() if c == col][0])
    return r, X, y


def new_rows(df=N, outcome_ok="elig", col=None):
    m = df[outcome_ok] & df[col].notna()
    return df[m].reset_index(drop=True)


# ------------------------------------------------------------------ record of what is in the table
w("RA-24: upstream trough depth and bombs, 1979-80..2003-04 (ERA5 proxy, pipeline A re-tracked on pressure only)")
w(f"generated {datetime.datetime.utcnow():%Y-%m-%d %H:%MZ}")
w()
for nm, D in (("new seasons 1979-80..2003-04", N), ("reference seasons 2004-05..2014-15", R)):
    E = D[D.elig]
    w(f"{nm}: {len(D)} fixes, {len(E)} class-eligible, {int(E.bomb.sum())} bomb fixes; "
      f"no sector cells (trough missing): {int(E.trough_up.isna().sum())} eligible fixes dropped")
w()

# ------------------------------------------------------------------ diagnostic 8: bomb frequency by basin and decade
Dg = pd.concat([N.assign(blk="new"), R.assign(blk="ref")])
Dg = Dg[Dg.elig]
Dg["decade"] = (Dg.season // 5) * 5
diag = Dg.groupby(["basin", "decade"]).agg(fixes=("bomb", "size"), bombs=("bomb", "sum"), rate=("bomb", "mean"),
                                           msl=("msl", "mean"), lat=("lat", "mean")).round(4).reset_index()
diag.to_csv("results/diag_bomb_rate.csv", index=False)
w("Diagnostic (selects nothing): bomb frequency of eligible fixes by basin and 5-season block")
w(diag.to_string(index=False))
w()

# ------------------------------------------------------------------ primary family: pooled, new seasons
FOREST = []
prim = {}
for d, col in T.DEFS.items():
    S = new_rows(col=col)
    r, X, y = one(S, col, "pooled 1979-2003", "primary")
    r["boot_lo"], r["boot_hi"] = T.season_boot(S, y, X, B=NB, seed=1)
    rowsP.append(r)
    prim[d] = (S, X, y, r)
P = pd.DataFrame(rowsP)
P["q_family"] = T.bh(P.p.values)

# ------------------------------------------------------------------ secondary family
for d, col in T.DEFS.items():
    for b in ("atl", "pac"):
        S = new_rows(col=col)
        S = S[S.basin == b].reset_index(drop=True)
        r, X, y = one(S, col, f"{b} 1979-2003", "secondary", pooled=False)
        rowsS.append(r)
    S = new_rows(col=col)
    S = S[S.season <= 2000].reset_index(drop=True)
    r, X, y = one(S, col, "pooled 1979-80..2000-01", "secondary")
    rowsS.append(r)
Sx = pd.DataFrame(rowsS)
Sx["q_family"] = T.bh(Sx.p.values)
A = pd.concat([P, Sx], ignore_index=True)
A["q_all"] = T.bh(A.p.values)
A = A[["family", "definition", "label", "OR", "lo", "hi", "p", "q_family", "q_all", "boot_lo", "boot_hi",
       "n", "events", "seasons"]]
A.to_csv("results/tests.csv", index=False)

# ------------------------------------------------------------------ sensitivities (not in a family)
for d, col in T.DEFS.items():
    S0 = new_rows(col=col)
    base = dict(pooled=True, year=True)
    runs = [("jet_up added", S0, dict(extra=("jet_up",))),
            ("no year term", S0, dict(year=False)),
            ("basin x covariate interactions", S0, dict(inter=True)),
            ("Oct-Apr fixes only", S0[S0.month.isin([10, 11, 12, 1, 2, 3, 4])].reset_index(drop=True), {})]
    for lab, S, kw in runs:
        S = S.dropna(subset=["jet_up"]) if "extra" in kw else S
        r, X, y = one(S, col, lab, "sensitivity", **kw)
        rowsX.append(r)
    S = N[N.elig_l24 & N[col].notna()].reset_index(drop=True)
    r, X, y = one(S, col, "lead 24 (outcome in t+24..t+48 h)", "sensitivity", outcome="bomb_l24")
    rowsX.append(r)
    Rr = R[R.elig & R[col].notna()].reset_index(drop=True)
    r, X, y = one(Rr, col, "reference seasons 2004-2014 alone (not fresh)", "sensitivity")
    rowsX.append(r)
    Np = pd.concat([N, R], ignore_index=True)
    Sp = Np[Np.elig & Np[col].notna()].reset_index(drop=True)
    r, X, y = one(Sp, col, "pooled 1979-2014 (36 seasons; 2004-14 not fresh)", "sensitivity")
    rowsX.append(r)
Xs = pd.DataFrame(rowsX)[["definition", "label", "OR", "lo", "hi", "p", "n", "events", "seasons"]]
Xs.to_csv("results/sensitivities.csv", index=False)

# ------------------------------------------------------------------ sensitivity 7: spline vs linear, leave-one-season-out
def rcs(z, knots):
    k = len(knots)
    sc = (knots[-1] - knots[0]) ** 2
    out = [z]
    for j in range(k - 2):
        b = (np.maximum(z - knots[j], 0) ** 3
             - np.maximum(z - knots[k - 2], 0) ** 3 * (knots[k - 1] - knots[j]) / (knots[k - 1] - knots[k - 2])
             + np.maximum(z - knots[k - 1], 0) ** 3 * (knots[k - 2] - knots[j]) / (knots[k - 1] - knots[k - 2]))
        out.append(b / sc)
    return np.column_stack(out)


nl = []
for d, col in T.DEFS.items():
    S, X, y, _ = prim[d]
    z = X[:, -1]
    kn = np.percentile(z, (1, 25, 50, 75, 99))
    Xs_ = np.column_stack([X[:, :-1], rcs(z, kn)])
    sea = S.season.values
    ll_lin, ll_spl = np.zeros(len(y)), np.zeros(len(y))
    for s in np.unique(sea):
        tr, te = sea != s, sea == s
        for XX, out in ((X, ll_lin), (Xs_, ll_spl)):
            b, _ = T.fit_w(XX[tr], y[tr])
            p = 1 / (1 + np.exp(-np.clip(XX[te] @ b, -30, 30)))
            p = np.clip(p, 1e-9, 1 - 1e-9)
            out[te] = -(y[te] * np.log(p) + (1 - y[te]) * np.log(1 - p))
    us = np.unique(sea)
    g = np.array([np.sum(ll_lin[sea == s] - ll_spl[sea == s]) for s in us])     # positive favours the spline
    rng = np.random.default_rng(3)
    bs = g[rng.integers(0, len(us), (NB, len(us)))].sum(1)
    nl.append(dict(definition=d, spline_gain_logloss=float(g.sum()), lo90=float(np.percentile(bs, 5)),
                   hi90=float(np.percentile(bs, 95)), p_nonlinear=float((np.sum(bs <= 0) + 1) / (NB + 1))))
NL = pd.DataFrame(nl)
NL.to_csv("results/nonlinear_loso.csv", index=False)

# ------------------------------------------------------------------ verdict (rule fixed in the pre-registration)
pw = pd.read_csv("results/power.csv") if os.path.exists("results/power.csv") else None
p1 = pd.DataFrame(rowsP).set_index("definition")
q1 = dict(zip(P.definition, P.q_family))
yes = all((p1.loc[d, "OR"] > 1) and (q1[d] < 0.05) for d in T.DEFS)
tight = all(p1.loc[d, "hi"] < 1.10 for d in T.DEFS)
pw115 = None
if pw is not None:
    pw115 = {d: float(pw[(pw.definition == d) & (pw.planted_OR == 1.15)].power.iloc[0]) for d in T.DEFS}
no = tight and pw115 is not None and all(v >= 0.8 for v in pw115.values())
sig = {d: (p1.loc[d, "OR"] > 1) and (q1[d] < 0.05) for d in T.DEFS}
if yes:
    verdict = "YES"
elif no:
    verdict = "NO (well-powered)"
elif sig["D1"] != sig["D2"] or (np.sign(np.log(p1.loc["D1", "OR"])) != np.sign(np.log(p1.loc["D2", "OR"]))):
    verdict = "CAN'T TELL (definitions disagree)"
else:
    verdict = "CAN'T TELL"

# ------------------------------------------------------------------ report
w("Primary family (pooled basins, new seasons; OR per reference SD; CR1 season-clustered 95% interval; 2-sided p)")
w(P[["definition", "OR", "lo", "hi", "p", "q_family", "boot_lo", "boot_hi", "n", "events", "seasons"]].round(4).to_string(index=False))
w()
w("Secondary family")
w(Sx[["definition", "label", "OR", "lo", "hi", "p", "q_family", "n", "events", "seasons"]].round(4).to_string(index=False))
w()
w("All 8 tests with q across all:")
w(A[["family", "definition", "label", "OR", "lo", "hi", "p", "q_family", "q_all"]].round(4).to_string(index=False))
npass_f = int(((A.q_all < 0.05) & (A.OR > 1)).sum())
w(f"Tests passing FDR across all 8 (q < 0.05, OR > 1): {npass_f} of {len(A)}")
w()
w("Sensitivities (not in a family)")
w(Xs.round(4).to_string(index=False))
w()
w("Spline against linear, leave-one-season-out log-loss gain (positive favours the spline), season bootstrap")
w(NL.round(4).to_string(index=False))
w()
if pw is not None:
    w(open("results/power.txt").read())
w(f"Power at planted OR 1.15: {pw115}")
w(f"Agenda prediction: pooled OR 1.10-1.25 per SD under D1 with an interval excluding 1.")
w()
w(f"VERDICT by the pre-registered rule: {verdict}")

with open("results/trough_bombs-result.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

# ------------------------------------------------------------------ looks
os.makedirs("../looks", exist_ok=True)
with open("../looks/ra24-trough-bombs.log", "a") as f:
    f.write(f"{datetime.datetime.utcnow().isoformat(timespec='seconds')}Z RA-24 trough depth and bombs: scored the "
            f"1979-80..2003-04 block once (25 seasons, 22 of them pre-2001: counts as one pre-2001 look); "
            f"2015-26 not scored (0 looks spent there); reference seasons 2004-05..2014-15 used as baseline and "
            f"as a labelled non-fresh sensitivity.\n")

# ------------------------------------------------------------------ figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7.5, 4.8))
Z = pd.concat([A.assign(kind=A.family), Xs.assign(kind="sensitivity")], ignore_index=True)
Z = Z.reset_index(drop=True)
yy = np.arange(len(Z))[::-1]
for i, r in Z.iterrows():
    c = "#1f77b4" if r.definition == "D1" else "#d95f02"
    ax.plot([r.lo, r.hi], [yy[i]] * 2, color=c, lw=1.6)
    ax.plot(r.OR, yy[i], "o", color=c, ms=4)
ax.axvline(1, color="k", lw=0.8)
ax.axvspan(1.10, 1.25, color="0.9", zorder=0)
ax.set_yticks(yy)
ax.set_yticklabels([f"{r.definition} {r.label}" for _, r in Z.iterrows()], fontsize=6)
ax.set_xlabel("odds ratio of a bomb per SD of upstream trough depth (95% CI, season-clustered)")
ax.set_title("RA-24 (ERA5 proxy, pipeline A): D1 blue, D2 orange; grey band = agenda prediction 1.10-1.25", fontsize=8)
fig.tight_layout()
fig.savefig("results/or_forest.png", dpi=130)
