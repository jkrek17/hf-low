"""Why does a strong Greenland high lower the Atlantic HF share? ERA5 PROXY, pipeline A. EXPLORATORY.
Plan: PREREGISTRATION.md. usage: run.py ALL_TRACKS CPC_DIR REPO_ROOT GH_DAILY FIXES ENV OUT_DIR [N]"""
import json, os, sys
import numpy as np, pandas as pd
import sg
from sg import S

SEED = sg.SEED
OUT_B = [  # id, column, kind, label
    ("B1", "lat0", "ols", "first-position latitude (deg N)"), ("B2", "lon0", "ols", "first-position longitude (deg E)"),
    ("B3", "latm", "ols", "min-pressure-fix latitude (deg N)"), ("B4", "lonm", "ols", "min-pressure-fix longitude (deg E)"),
    ("B5", "peak_lat", "ols", "peak-gust latitude (deg N)"), ("B6", "lonp", "ols", "peak-gust longitude (deg E)"),
    ("B7", "minp", "ols", "track minimum pressure (hPa)"), ("B8", "gust800_kt", "ols", "track peak 800 km gust (kt)"),
    ("B9", "ndr_max", "ols", "max 24 h normalised deepening rate"), ("B10", "box", "logit", "min-pressure fix in deep-low box (log odds)"),
    ("B11", "jet0", "ols", "250 hPa jet at first position (m/s)"), ("B12", "eady0", "ols", "Eady growth at first position")]
MED = {"C1": ["lat0", "lon0", "lat0_2", "lon0_2"],
       "C2": ["lat0", "lon0", "lat0_2", "lon0_2", "minp", "ndr_max", "ndr_na"],
       "C3": ["lat0", "lon0", "lat0_2", "lon0_2", "jet0", "eady0", "sstgrad0", "sstgrad_na", "tcwv0"],
       "C4": ["lat0", "lon0", "lat0_2", "lon0_2", "minp", "ndr_max", "ndr_na", "jet0", "eady0", "sstgrad0", "sstgrad_na", "tcwv0"]}


def logit(X, y, iters=40):
    b = np.zeros(X.shape[1])
    b[-1] = np.log(max(y.mean(), 1e-6) / max(1 - y.mean(), 1e-6))
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(X @ b, -30, 30)))
        w = p * (1 - p) + 1e-9
        H = X.T @ (X * w[:, None]) + 1e-8 * np.eye(X.shape[1])
        step = np.linalg.solve(H, X.T @ (y - p))
        b += step
        if np.max(np.abs(step)) < 1e-9:
            break
    return b


def ols(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]


def main():
    tracks, cpc, repo, ghfile, fixes, env, out = sys.argv[1:8]
    nboot = int(sys.argv[8]) if len(sys.argv) > 8 else 2000
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(SEED)
    P = sg.prep(tracks, cpc, repo, ghfile, fixes, env)
    D, z, Tw = P["D"], P["z"], P["Tw"]
    seasons = np.arange(sg.S0, sg.S1 + 1)

    # ---------------- daily counts (Poisson), as PR 63
    C = S.daily_counts(Tw, D, "atl", ["all", "hf"])
    Xg, Xn, Xgn = (S.design(C, z, c) for c in (["GH_lag"], ["NAO_lag"], ["NAO_lag", "GH_lag"]))
    Xg, Xn, Xgn = (np.asarray(x, float) for x in (Xg, Xn, Xgn))  # Xgn columns: NAO, GH
    yd = {k: np.asarray(C[k], float) for k in ("all", "hf")}
    d_season = C.season.values
    d_idx = {s: np.where(d_season == s)[0] for s in seasons}

    # ---------------- storm table
    W = Tw[Tw.lat0.notna()].reset_index(drop=True).copy()  # tracks with fix rows
    W["lat0_2"], W["lon0_2"] = W.lat0 ** 2, W.lon0 ** 2
    W["ndr_na"] = W.ndr_max.isna().astype(float)
    W["ndr_max"] = W.ndr_max.fillna(0.0)
    W["sstgrad_na"] = W.sstgrad0.isna().astype(float)  # 27% of storms: first position over ice/land mask; filled with the mean, flagged
    W["sstgrad0"] = W.sstgrad0.fillna(W.sstgrad0.mean())
    zc = ["lat0", "lon0", "lat0_2", "lon0_2", "minp", "ndr_max", "jet0", "eady0", "sstgrad0", "tcwv0"]
    for c in zc:  # standardise mediators once on the full sample (nuisance scaling only)
        W[c + "_s"] = (W[c] - W[c].mean()) / W[c].std(ddof=0)
    # ndr_max/ndr_na for B9 itself is the raw (un-imputed) value
    W["ndr_raw"] = Tw.set_index("track").loc[W.track, "ndr_max"].values
    mons = sorted(W.mon.unique())[1:]
    Z = np.column_stack([(W.mon.values == m).astype(float) for m in mons] + [(W.season.values - W.season.mean()) / 10.0, np.ones(len(W))])
    gh_, na_ = W.zGH.values, W.zNAO.values
    X_gn, X_g, X_n = (np.column_stack(c + [Z]) for c in ([gh_, na_], [gh_], [na_]))  # GH is col 0, NAO col 1 in X_gn
    s_season = W.season.values
    s_idx = {s: np.where(s_season == s)[0] for s in seasons}
    hf = W.hf.values.astype(float)
    ycol = {i: (W.ndr_raw.values if c == "ndr_max" else W[c].values.astype(float)) for i, c, k, l in OUT_B}
    Msets = {k: np.column_stack([W[c + "_s"].values if c + "_s" in W else W[c].values for c in v]) for k, v in MED.items()}

    def stats(drows, srows):
        o = {}
        # --- A2: Poisson, per +1 SD, rows are daily
        Xa, Xb, Xc = Xg[drows], Xn[drows], Xgn[drows]
        yy = {k: v[drows] for k, v in yd.items()}
        for tag, X, col in (("GH", Xa, 0), ("NAO", Xb, 0), ("GHgNAO", Xc, 1), ("NAOgGH", Xc, 0)):
            ba = S.fast_pois(yy["all"], X)[col]
            bh = S.fast_pois(yy["hf"], X)[col]
            o[f"A2_{tag}_all"], o[f"A2_{tag}_hf"], o[f"A2_{tag}_share"] = ba, bh, bh - ba
        # --- B: storm level
        Zs = Z[srows]
        for i, c, k, l in OUT_B:
            y = ycol[i][srows]
            m = ~np.isnan(y)
            for tag, X, col in (("GHgNAO", X_gn, 0), ("GH", X_g, 0), ("NAOgGH", X_gn, 1)):
                Xm = X[srows][m]
                o[f"{i}_{tag}"] = (ols if k == "ols" else logit)(Xm, y[m])[col]
        # --- C
        yh = hf[srows]
        o["C0_GH_NAOin"] = logit(X_gn[srows], yh)[0]
        o["C0_GH_NAOout"] = logit(X_g[srows], yh)[0]
        for k, M in Msets.items():
            for tag, base in (("NAOin", X_gn), ("NAOout", X_g)):
                Xm = np.column_stack([base[srows][:, :base.shape[1] - Z.shape[1]], M[srows], Zs])
                o[f"{k}_GH_{tag}"] = logit(Xm, yh)[0]
        return o

    alld, alls = np.arange(len(C)), np.arange(len(W))
    est = stats(alld, alls)
    n_info = dict(n_days=len(C), n_tracks=int(C["all"].sum()), n_hf=int(C.hf.sum()), n_storm_rows=len(W),
                  n_storm_hf=int(W.hf.sum()), n_tracks_without_fixes=int(len(Tw) - len(W)),
                  n_ndr_missing=int(W.ndr_na.sum()), seed=SEED, n_boot=nboot)

    boots = []
    for _ in range(nboot):
        pick = rng.choice(seasons, len(seasons), replace=True)
        boots.append(stats(np.concatenate([d_idx[s] for s in pick]), np.concatenate([s_idx[s] for s in pick])))
    B = pd.DataFrame(boots)
    loso = []
    for s in seasons:
        loso.append(stats(np.where(d_season != s)[0], np.where(s_season != s)[0]))
    L = pd.DataFrame(loso)

    # ---------------- A1, A3, A4
    zr = z.corr().loc["NAO_lag", "GH_lag"]
    r_storm = float(np.corrcoef(W.zNAO, W.zGH)[0, 1])
    sm = pd.DataFrame({"s": C.season.values, "n": z.NAO_lag.values, "g": z.GH_lag.values}).groupby("s").mean()
    r_seas = float(np.corrcoef(sm.n, sm.g)[0, 1])
    a1 = dict(r_daily=float(zr), r_storm=r_storm, r_season_means=r_seas, vif=float(1 / (1 - zr ** 2)),
              partial_sd_GH_given_NAO=float(np.sqrt(1 - zr ** 2)), sd_GH_hPa=float(P["raw"].GH_lag.std(ddof=0)), sd_NAO=float(P["raw"].NAO_lag.std(ddof=0)))
    # joint region
    pts = B[["A2_GHgNAO_share", "A2_NAOgGH_share"]].dropna().values
    cov, mean = np.cov(pts.T), pts.mean(0)
    corners = {"GH only (GH|NAO = GH alone, NAO|GH = 0)": (est["A2_GH_share"], 0.0),
               "NAO only (GH|NAO = 0, NAO|GH = NAO alone)": (0.0, est["A2_NAO_share"]),
               "neither (0, 0)": (0.0, 0.0),
               "both at PR 63 NAO-given-GH 1.039 and GH 0.870 pairing": (est["A2_GHgNAO_share"], est["A2_NAOgGH_share"])}
    cinv = np.linalg.inv(cov)
    reg = {}
    for k, (a, b2) in corners.items():
        v = np.array([a, b2]) - mean
        d2 = float(v @ cinv @ v)
        reg[k] = dict(mahalanobis_d2=d2, inside_95=bool(d2 < 5.991))
    a1["joint_boot_corr_GHgNAO_NAOgGH"] = float(np.corrcoef(pts.T)[0, 1])
    a1["region"] = reg

    # ---------------- tests table
    rows = []

    def add(block, tid, desc, key, family, contrasts=()):
        bd = B[key].dropna()
        lo, hi = np.percentile(bd, [2.5, 97.5])
        p = max(2 * min((bd <= 0).mean(), (bd >= 0).mean()), 1 / (len(bd) + 1))
        r = dict(block=block, test=tid, desc=desc, key=key, family=family, est=est[key], lo=lo, hi=hi, p=p,
                 se=bd.std(ddof=1), mde=2.8 * bd.std(ddof=1), loso_min=L[key].min(), loso_max=L[key].max())
        for c in contrasts:
            r[c] = est.get(key.replace("GHgNAO", c).replace("GH", c) if False else c, np.nan)
        rows.append(r)

    for tag, nm in (("GH", "GH alone"), ("GHgNAO", "GH given NAO"), ("NAO", "NAO alone"), ("NAOgGH", "NAO given GH")):
        add("A", f"A2{'abcd'[['GH','GHgNAO','NAO','NAOgGH'].index(tag)]}", f"log share coefficient, {nm}", f"A2_{tag}_share", "A")
    for i, c, k, l in OUT_B:
        add("B", i, f"{l}: GH given NAO", f"{i}_GHgNAO", "B")
    # C: paired differences
    for k in ("C1", "C2", "C3", "C4"):
        B[f"{k}_diff"], B[f"{k}_frac"] = B[f"{k}_GH_NAOin"] - B["C0_GH_NAOin"], B[f"{k}_GH_NAOin"] / B["C0_GH_NAOin"]
        est[f"{k}_diff"], est[f"{k}_frac"] = est[f"{k}_GH_NAOin"] - est["C0_GH_NAOin"], est[f"{k}_GH_NAOin"] / est["C0_GH_NAOin"]
        L[f"{k}_diff"], L[f"{k}_frac"] = L[f"{k}_GH_NAOin"] - L["C0_GH_NAOin"], L[f"{k}_GH_NAOin"] / L["C0_GH_NAOin"]
    for k, d in (("C1", "STEER: first position"), ("C2", "WEAKEN: + min pressure, deepening"),
                 ("C3", "WEAKEN: + first-position environment"), ("C4", "all mediators")):
        add("C", k, f"GH coefficient, change on adding {d}", f"{k}_diff", "C")
    add("C", "C5", "GH coefficient with C2 mediators (sign; WEAKEN predicts >= 0)", "C2_GH_NAOin", "C")
    R = pd.DataFrame(rows)

    def bh(p):
        p = np.asarray(p); o = np.argsort(p); n = len(p); q = np.empty(n); prev = 1.0
        for rank, i in zip(range(n, 0, -1), o[::-1]):
            prev = min(prev, p[i] * n / rank); q[i] = prev
        return q
    R["q_family"] = np.nan
    for f in R.family.unique():
        m = (R.family == f).values
        R.loc[m, "q_family"] = bh(R.p[m])
    R["q_all"] = bh(R.p)
    # contrasts for B: GH alone and NAO given GH (not tested)
    for i, c, k, l in OUT_B:
        for tag in ("GH", "NAOgGH"):
            bd = B[f"{i}_{tag}"].dropna()
            R.loc[R.test == i, f"{tag}_est"] = est[f"{i}_{tag}"]
            R.loc[R.test == i, f"{tag}_lo"] = np.percentile(bd, 2.5)
            R.loc[R.test == i, f"{tag}_hi"] = np.percentile(bd, 97.5)
    R.to_csv(os.path.join(out, "tests.csv"), index=False)

    # ---------------- retained fractions with intervals
    head = dict(n_info, **{"A1": a1})
    def ci(x):
        return [float(v) for v in np.percentile(pd.Series(x).dropna(), [2.5, 97.5])]
    head["share_log_coef"] = {t: dict(est=est[f"A2_{t}_share"], ci=ci(B[f"A2_{t}_share"]), RR=float(np.exp(est[f"A2_{t}_share"])), RR_ci=[float(np.exp(v)) for v in ci(B[f"A2_{t}_share"])]) for t in ("GH", "GHgNAO", "NAO", "NAOgGH")}
    head["count_log_coef"] = {t: dict(hf=est[f"A2_{t}_hf"], hf_ci=ci(B[f"A2_{t}_hf"]), all=est[f"A2_{t}_all"], all_ci=ci(B[f"A2_{t}_all"])) for t in ("GH", "GHgNAO", "NAO", "NAOgGH")}
    head["C"] = {"C0_GH_NAOin": dict(est=est["C0_GH_NAOin"], ci=ci(B["C0_GH_NAOin"])), "C0_GH_NAOout": dict(est=est["C0_GH_NAOout"], ci=ci(B["C0_GH_NAOout"]))}
    for k in MED:
        head["C"][k] = dict(coef_NAOin=est[f"{k}_GH_NAOin"], coef_ci=ci(B[f"{k}_GH_NAOin"]), frac=est[f"{k}_frac"], frac_ci=ci(B[f"{k}_frac"]),
                            frac_loso=[float(L[f"{k}_frac"].min()), float(L[f"{k}_frac"].max())],
                            coef_NAOout=est[f"{k}_GH_NAOout"], coef_NAOout_ci=ci(B[f"{k}_GH_NAOout"]))
    head["C5_frac_boot_ge0"] = float((B["C2_GH_NAOin"] >= 0).mean())
    # NAO retained fractions for the stand-in rule
    head["standin_NAOgGH_over_NAO"] = est["A2_NAOgGH_share"] / est["A2_NAO_share"]
    head["standin_NAOgGH_over_NAO_ci"] = ci(B["A2_NAOgGH_share"] / B["A2_NAO_share"])
    head["standin_GHgNAO_over_GH"] = est["A2_GHgNAO_share"] / est["A2_GH_share"]
    head["standin_GHgNAO_over_GH_ci"] = ci(B["A2_GHgNAO_share"] / B["A2_GH_share"])
    json.dump(head, open(os.path.join(out, "headline.json"), "w"), indent=1, default=float)
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.float_format", lambda v: f"{v:.4g}")
    with open(os.path.join(out, "summary.txt"), "w") as f:
        f.write("Greenland high and the Atlantic HF share. Pipeline A (ERA5 proxy), Oct-Apr 2004-05..2025-26, EXPLORATORY.\n\n")
        f.write(json.dumps(head, indent=1, default=float) + "\n\n")
        f.write(R[["block", "test", "desc", "est", "lo", "hi", "p", "q_family", "q_all", "mde", "loso_min", "loso_max"]].to_string(index=False) + "\n")
    print(open(os.path.join(out, "summary.txt")).read())


if __name__ == "__main__":
    main()
