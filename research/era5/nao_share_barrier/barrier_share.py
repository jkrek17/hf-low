"""Atlantic NAO share effect with barrier-type fixes removed, and the Greenland high as
mediator. ERA5 PROXY, pipeline A. Plan: PREREGISTRATION.md (committed first, d642ff7).

Reuses PR 14's machinery (research/era5/freq_split/split.py): same tracks, window,
lagged CPC NAO (days -10..-4), month effects, season trend, season-block resampling.

usage: barrier_share.py ALL_TRACKS CPC_DIR REPO_ROOT GH_DAILY GUSTLOC_FIXES OUT_DIR [N]
"""
import json, os, sys, warnings
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "freq_split"))
import split as S  # noqa: E402

warnings.filterwarnings("ignore")
SEED = 20261008
THR = 71.7
S0, S1, WIN = 2004, 2025, "octapr"


# ------------------------------------------------------------------ track flags
def flags(T, fixes):
    """Per-track HF-without-barrier flags from the fix-level table (HF-strength fixes only)."""
    F = fixes.copy()
    F["T"] = ~((F.max_dkm > 400) & (F.max_dgl < 300))
    g = lambda c: F[c].fillna(-1).ge(THR)
    keep = {"T": F["T"], "G100": g("g800_gl100"), "G300": g("g800_gl300"), "R400": g("g800_r400")}
    out = {}
    for k, m in keep.items():
        out[k] = set(F.loc[m, "track"])
    T = T.copy()
    has_rows = T.track.isin(set(F.track))
    for k in keep:
        # an HF track with no fix rows cannot be classified: kept as HF (logged)
        T[k] = ((T.hf == 1) & T.track.isin(out[k])) | ((T.hf == 1) & ~has_rows)
        T[k] = T[k].astype(int)
    T["B"] = ((T.hf == 1) & (T["T"] == 0)).astype(int)
    T["unclassified_hf"] = ((T.hf == 1) & ~has_rows).astype(int)
    return T


# ------------------------------------------------------------------ estimates
def b(C, X, col, k=0):
    return S.fast_pois(np.asarray(C[col], float), np.asarray(X, float))[k]


def stats(C, X1, X2, X3):
    """All log-scale statistics. X1: NAO; X2: NAO + GH; X3: GH."""
    c1 = ["all", "hf", "T", "G100", "G300", "R400", "B", "all60", "hf60", "allD", "TD"]
    b1 = {c: b(C, X1, c) for c in c1}
    b2 = {c: b(C, X2, c) for c in ("all", "hf", "T")}
    b3 = {c: b(C, X3, c) for c in ("all", "hf", "T")}
    sh = b1["hf"] - b1["all"]
    o = {}
    o["orig_share"] = sh
    o["P1"] = b1["T"] - b1["all"]
    o["P2"] = o["P1"] - sh
    o["P3_share_NAO|GH"] = b2["hf"] - b2["all"]
    o["P3"] = o["P3_share_NAO|GH"] - sh
    for v in ("G100", "G300", "R400"):
        o[f"F2_{v}"] = b1[v] - b1["all"]
        o[f"F2_{v}_diff"] = o[f"F2_{v}"] - sh
    o["F2_N60"] = b1["hf60"] - b1["all60"]
    o["F2_N60_diff"] = o["F2_N60"] - sh
    o["F2_D"] = b1["TD"] - b1["allD"]
    o["F2_D_diff"] = o["F2_D"] - sh
    o["F3_stormHF_NAO"] = b1["T"]
    o["F3_barrierHF_NAO"] = b1["B"]
    o["F3_diff"] = b1["T"] - b1["B"]
    o["F4_GH_share_orig"] = b3["hf"] - b3["all"]
    o["F4_GH_share_T"] = b3["T"] - b3["all"]
    o["F4_GH_diff"] = o["F4_GH_share_T"] - o["F4_GH_share_orig"]
    o["F4_NAO|GH_share_T"] = b2["T"] - b2["all"]
    o["logRR_HF_NAO"] = b1["hf"]
    o["logRR_all_NAO"] = b1["all"]
    o["q_retained"] = o["P1"] / sh
    o["r_med"] = o["P3_share_NAO|GH"] / sh
    return o


PERM1 = ["P1", "F2_G100", "F2_G300", "F2_R400", "F2_N60", "F2_D", "F3_stormHF_NAO", "F3_barrierHF_NAO"]
PERM3 = ["F4_GH_share_orig", "F4_GH_share_T"]


def main():
    tracks, cpc, repo, ghfile, fixfile, out = sys.argv[1:7]
    nboot = int(sys.argv[7]) if len(sys.argv) > 7 else 2000
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(SEED)
    T = S.load_tracks(tracks)
    T = T[T.basin == "atl"].reset_index(drop=True)
    fixes = pd.read_csv(fixfile)
    T = flags(T, fixes)
    I = S.indices(cpc, repo)
    gh = pd.read_csv(ghfile, parse_dates=["date"]).set_index("date").GH.asfreq("D")
    I["GH_lag"] = gh.rolling(7, min_periods=7).mean().shift(4)

    D = S.window_days(WIN, S0, S1)
    raw = I.reindex(D.date)[["NAO_lag", "GH_lag"]].reset_index(drop=True)
    ok = raw.notna().all(axis=1).values
    D, raw = D[ok].reset_index(drop=True), raw[ok].reset_index(drop=True)
    z = (raw - raw.mean()) / raw.std(ddof=0)
    Tw = S.assign(T, WIN, S0, S1)
    # extra samples
    n60 = Tw.peak_lat < 60
    dropD = (Tw.hf == 1) & (Tw["T"] == 0)
    C = S.daily_counts(Tw, D, "atl", ["all", "hf", "T", "G100", "G300", "R400", "B"])
    C60 = S.daily_counts(Tw[n60], D, "atl", ["all", "hf"]).rename(columns={"all": "all60", "hf": "hf60"})
    CD = Tw[~dropD].copy()
    CD = S.daily_counts(CD, D, "atl", ["all", "T"]).rename(columns={"all": "allD", "T": "TD"})
    for c in ("all60", "hf60"):
        C[c] = C60[c].values
    for c in ("allD", "TD"):
        C[c] = CD[c].values
    X1, X2, X3 = (S.design(C, z, cols) for cols in (["NAO_lag"], ["NAO_lag", "GH_lag"], ["GH_lag"]))
    seasons = np.arange(S0, S1 + 1)
    sidx = {s: np.where(C.season.values == s)[0] for s in seasons}

    est = stats(C, X1, X2, X3)
    info = dict(n_days=len(C), n_seasons=len(seasons),
                n_tracks=int(C["all"].sum()), n_hf=int(C.hf.sum()), n_T=int(C["T"].sum()),
                n_G100=int(C.G100.sum()), n_G300=int(C.G300.sum()), n_R400=int(C.R400.sum()), n_B=int(C.B.sum()),
                n_tracks_N60=int(C.all60.sum()), n_hf_N60=int(C.hf60.sum()),
                n_tracks_D=int(C.allD.sum()), n_hf_D=int(C.TD.sum()),
                unclassified_hf_in_window=int(Tw.unclassified_hf.sum()),
                corr_NAO_GH_lag=float(np.corrcoef(z.NAO_lag, z.GH_lag)[0, 1]),
                sd_NAO=float(raw.NAO_lag.std(ddof=0)), sd_GH_hPa=float(raw.GH_lag.std(ddof=0)))

    # ---- season-block bootstrap (all statistics jointly)
    boots = []
    for _ in range(nboot):
        pick = rng.choice(seasons, len(seasons), replace=True)
        rows = np.concatenate([sidx[s] for s in pick])
        Cb = C.iloc[rows].reset_index(drop=True)
        try:
            boots.append(stats(Cb, X1.iloc[rows].reset_index(drop=True),
                               X2.iloc[rows].reset_index(drop=True), X3.iloc[rows].reset_index(drop=True)))
        except Exception:
            continue
    B = pd.DataFrame(boots)

    # ---- season-block permutation for single coefficients
    def permute(col_names, X, names):
        blocks = {s: z.loc[C.season.values == s, col_names].values for s in seasons}
        lens = {s: len(v) for s, v in blocks.items()}
        ge = {n: 0 for n in names}
        for _ in range(nboot):
            perm = rng.permutation(seasons)
            Zp = np.vstack([np.resize(blocks[p], (lens[s], len(col_names))) if len(blocks[p]) != lens[s] else blocks[p]
                            for s, p in zip(seasons, perm)])
            Xp = X.copy()
            Xp[[c for c in col_names]] = Zp
            return_stats = stats_perm(C, Xp, names, col_names)
            for n in names:
                ge[n] += abs(return_stats[n]) >= abs(est[n]) - 1e-12
        return {n: (ge[n] + 1) / (nboot + 1) for n in names}

    def stats_perm(C_, Xp, names, col_names):
        if col_names == ["NAO_lag"]:
            b1 = {c: b(C_, Xp, c) for c in ("all", "hf", "T", "G100", "G300", "R400", "B", "all60", "hf60", "allD", "TD")}
            o = {"P1": b1["T"] - b1["all"], "F2_G100": b1["G100"] - b1["all"], "F2_G300": b1["G300"] - b1["all"],
                 "F2_R400": b1["R400"] - b1["all"], "F2_N60": b1["hf60"] - b1["all60"], "F2_D": b1["TD"] - b1["allD"],
                 "F3_stormHF_NAO": b1["T"], "F3_barrierHF_NAO": b1["B"]}
        else:
            b3 = {c: b(C_, Xp, c) for c in ("all", "hf", "T")}
            o = {"F4_GH_share_orig": b3["hf"] - b3["all"], "F4_GH_share_T": b3["T"] - b3["all"]}
        return o

    pperm = {}
    pperm.update(permute(["NAO_lag"], X1, PERM1))
    pperm.update(permute(["GH_lag"], X3, PERM3))

    # ---- leave one season out for the headline quantities
    loso = []
    for s in seasons:
        keep = C.season.values != s
        loso.append(stats(C[keep].reset_index(drop=True), X1[keep].reset_index(drop=True),
                          X2[keep].reset_index(drop=True), X3[keep].reset_index(drop=True)))
    L = pd.DataFrame(loso)

    # ---- tests table
    family = {"P1": "F1", "P2": "F1", "P3": "F1",
              "F2_G100": "F2", "F2_G100_diff": "F2", "F2_G300": "F2", "F2_G300_diff": "F2", "F2_R400": "F2",
              "F2_R400_diff": "F2", "F2_N60": "F2", "F2_N60_diff": "F2", "F2_D": "F2", "F2_D_diff": "F2",
              "F3_stormHF_NAO": "F3", "F3_barrierHF_NAO": "F3", "F3_diff": "F3",
              "F4_GH_share_orig": "F4", "F4_GH_share_T": "F4", "F4_GH_diff": "F4", "F4_NAO|GH_share_T": "F4"}
    rows = []
    for n, fam in family.items():
        bd = B[n].dropna()
        lo, hi = np.percentile(bd, [2.5, 97.5])
        if n in pperm:
            p, kind = pperm[n], "permutation"
        else:
            p, kind = 2 * min((bd <= 0).mean(), (bd >= 0).mean()), "bootstrap"
            p = max(p, 1 / (len(bd) + 1))
        rows.append(dict(test=n, family=fam, log_est=est[n], RR=np.exp(est[n]), RR_lo=np.exp(lo), RR_hi=np.exp(hi),
                         p=p, p_kind=kind, se_boot=bd.std(ddof=1), mde_2p8se=2.8 * bd.std(ddof=1),
                         loso_min=L[n].min(), loso_max=L[n].max()))
    R = pd.DataFrame(rows)

    def bh(p):
        p = np.asarray(p); o = np.argsort(p); n = len(p)
        q = np.empty(n); prev = 1.0
        for rank, i in zip(range(n, 0, -1), o[::-1]):
            prev = min(prev, p[i] * n / rank); q[i] = prev
        return q
    R["q_family"] = np.nan
    for fam in R.family.unique():
        m = (R.family == fam).values
        R.loc[m, "q_family"] = bh(R.p[m])
    R["q_all"] = bh(R.p)
    R.to_csv(os.path.join(out, "tests.csv"), index=False)

    # ---- headline retention and mediation, with intervals
    def ci(x):
        return np.percentile(x.dropna(), [2.5, 97.5])
    head = dict(info)
    head.update(orig_RR_share=float(np.exp(est["orig_share"])), orig_RR_share_ci=[float(np.exp(v)) for v in ci(B["orig_share"])],
                orig_RR_HF=float(np.exp(est["logRR_HF_NAO"])), orig_RR_all=float(np.exp(est["logRR_all_NAO"])),
                q_retained=est["q_retained"], q_ci=[float(v) for v in ci(B["q_retained"])],
                r_med=est["r_med"], r_ci=[float(v) for v in ci(B["r_med"])],
                q_loso=[float(L.q_retained.min()), float(L.q_retained.max())],
                r_loso=[float(L.r_med.min()), float(L.r_med.max())],
                RR_share_NAO_given_GH=float(np.exp(est["P3_share_NAO|GH"])),
                RR_share_NAO_given_GH_ci=[float(np.exp(v)) for v in ci(B["P3_share_NAO|GH"])],
                n_boot=len(B), seed=SEED)
    for v in ("G100", "G300", "R400", "N60", "D"):
        qv = est[f"F2_{v}"] / est["orig_share"]
        head[f"q_{v}"] = qv
        head[f"q_{v}_ci"] = [float(x) for x in ci(B[f"F2_{v}"] / B["orig_share"])]
    head["q_ci_T_boot_median"] = float((B["P1"] / B["orig_share"]).median())
    json.dump(head, open(os.path.join(out, "headline.json"), "w"), indent=1)
    pd.set_option("display.width", 250, "display.max_columns", 30, "display.float_format", lambda v: f"{v:.4g}")
    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write("NAO share effect without barrier winds. Atlantic, pipeline A (ERA5 proxy), Oct-Apr 2004-05..2025-26, NAO lagged days -10..-4.\n\n")
        f.write(json.dumps(head, indent=1) + "\n\n")
        f.write(R.to_string(index=False) + "\n")
    print(open(os.path.join(out, "summary.txt")).read())


if __name__ == "__main__":
    main()
