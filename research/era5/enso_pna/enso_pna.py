"""Does ENSO (monthly ONI) change what PNA does to Pacific HF lows? ERA5 PROXY, pipeline A.

Plan: PREREGISTRATION.md (committed before any fit, c6a085d). Reuses helpers from
research/era5/freq_split/split.py.

Models (daily genesis counts, Poisson, month FE + season trend [+ era term], season-clustered SEs):
  M0 y ~ ONI                  total effect T
  M1 y ~ ONI + PNA            direct effect D, PNA effect b
  M2 y ~ ONI + PNA + ONIxPNA  interaction gamma
  A  PNA ~ ONI (OLS)          path a;  indirect effect = a x b
Each is fitted for the HF count, the all-cyclone count and (stacked, share = HF - all).

usage: enso_pna.py ALL_TRACKS_CSV_GZ CPC_DIR REPO_ROOT OUT_DIR [NPERM]
"""
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "freq_split"))
import split as S  # noqa: E402

warnings.filterwarnings("ignore")
SEED = 20261008
BASIN = "pac"


# ------------------------------------------------------------- estimation
def fast_pois(y, X, iters=40):
    b = np.zeros(X.shape[1])
    b[-1] = np.log(max(y.mean(), 1e-9))
    for _ in range(iters):
        mu = np.exp(X @ b)
        step = np.linalg.solve(X.T @ (X * mu[:, None]), X.T @ (y - mu))
        b += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b


def cluster_cov(X, w, resid, G, n_per):
    """Sandwich covariance, clusters = G equal-length contiguous season blocks, with
    statsmodels' default small-sample factor G/(G-1) * (N-1)/(N-k)."""
    N, k = X.shape
    bread = np.linalg.inv(X.T @ (X * w[:, None]))
    sc = (X * resid[:, None]).reshape(G, n_per, k).sum(1)
    meat = sc.T @ sc
    return bread @ meat @ bread * (G / (G - 1)) * ((N - 1) / (N - k))


def pois_fit(y, X, G, n_per):
    b = fast_pois(y, X)
    mu = np.exp(X @ b)
    return b, cluster_cov(X, mu, y - mu, G, n_per)


def ols_fit(y, X, G, n_per):
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    return b, cluster_cov(X, np.ones(len(y)), y - X @ b, G, n_per)


def stacked_fit(ya, yb, X, G, n_per):
    """Outcome a (denominator) and b stacked, every column interacted with the b indicator.
    Second block of coefficients = b - a, with clustered covariance (clusters span both halves)."""
    k = X.shape[1]
    Z = np.zeros_like(X)
    XX = np.vstack([np.hstack([X, Z]), np.hstack([X, X])])
    y = np.concatenate([ya, yb])
    b = fast_pois(y, XX)
    mu = np.exp(XX @ b)
    sc1 = (XX[: len(ya)] * (y[: len(ya)] - mu[: len(ya)])[:, None]).reshape(G, n_per, 2 * k).sum(1)
    sc2 = (XX[len(ya):] * (y[len(ya):] - mu[len(ya):])[:, None]).reshape(G, n_per, 2 * k).sum(1)
    sc = sc1 + sc2
    N = len(y)
    bread = np.linalg.inv(XX.T @ (XX * mu[:, None]))
    cov = bread @ (sc.T @ sc) @ bread * (G / (G - 1)) * ((N - 1) / (N - 2 * k))
    return b, cov


# ------------------------------------------------------------- data
class Data:
    """Standardised index columns and daily counts for one analysis cell."""

    def __init__(self, T, I, window, s0, s1, outcomes, now=False, era=False):
        self.window, self.s0, self.s1, self.now = window, s0, s1, now
        D = S.window_days(window, s0, s1)
        Tw = S.assign(T, window, s0, s1)
        Tw = Tw[Tw.basin == BASIN]
        pc = "PNA_now" if now else "PNA_lag"
        raw = I.reindex(D.date)[[pc, "ONI_lag"]].reset_index(drop=True)
        if now:  # ONI of the genesis month (day 0)
            oni = S.read_oni(CPC, REPO)
            raw["ONI_lag"] = [oni.get((d.year, d.month), np.nan) for d in D.date]
        self.n_missing = int(raw.isna().any(axis=1).sum())
        self.pna_raw, self.oni_raw = raw[pc].values, raw["ONI_lag"].values
        mu, sd = raw.mean(), raw.std(ddof=0)
        z = (raw - mu) / sd
        z = z.fillna(0.0)  # missing index day -> standardised mean, so every season keeps its full calendar
        self.pna, self.oni = z[pc].values, z["ONI_lag"].values
        self.sd_pna, self.sd_oni = float(sd[pc]), float(sd["ONI_lag"])
        self.mu_oni = float(mu["ONI_lag"])
        self.C = S.daily_counts(Tw, D, BASIN, outcomes)
        self.season = self.C.season.values
        self.G = s1 - s0 + 1
        self.n_per = len(self.C) // self.G
        B = pd.DataFrame({f"m{mo}": (self.C.month.values == mo).astype(float)
                          for mo in sorted(self.C.month.unique())[1:]})
        B["trend"] = (self.season - self.season.mean()) / 10.0
        if era:
            B["era"] = (self.season <= 2000).astype(float)
        B["const"] = 1.0
        self.B = B.values
        self.outcomes = outcomes
        # N+ / N-: ONI split at zero, standardised after the split
        o = self.oni_raw.copy()
        o = np.where(np.isnan(o), 0.0, o)
        self.nplus = _zs(np.maximum(o, 0))
        self.nminus = _zs(np.minimum(o, 0))


def _zs(x):
    return (x - x.mean()) / x.std()


# ------------------------------------------------------------- one set of estimates
def point(d, ysets, oni, pna, G=None):
    """Point estimates for every outcome: T, D, b, gamma, a, indirect. ysets: dict name -> y."""
    n = len(oni)
    B = d.B
    X0 = np.column_stack([oni, B])
    X1 = np.column_stack([oni, pna, B])
    X2 = np.column_stack([oni, pna, oni * pna, B])
    out = {}
    for o, y in ysets.items():
        T = fast_pois(y, X0)[0]
        b1 = fast_pois(y, X1)
        b2 = fast_pois(y, X2)
        out[o] = dict(T=T, D=b1[0], b=b1[1], gamma=b2[2], b_pna2=b2[1], b_oni2=b2[0])
    a = np.linalg.lstsq(X0, pna, rcond=None)[0][0]
    for o in ysets:
        out[o]["a"] = a
        out[o]["ind"] = a * out[o]["b"]
    return out


def wald_gamma(d, y, oni, pna):
    X2 = np.column_stack([oni, pna, oni * pna, d.B])
    b, cov = pois_fit(y, X2, d.G, d.n_per)
    return b[2], np.sqrt(cov[2, 2])


def season_blocks(x, G, n_per):
    return x.reshape(G, n_per)


def run_cell(d, name, nperm, rng, primary_outcome, comps=("hf", "all")):
    """All estimates for one analysis cell. primary_outcome = name of the HF-like column."""
    all_c, hf_c = "all", primary_outcome
    ys = {"hf": d.C[hf_c].values.astype(float), "all": d.C[all_c].values.astype(float)}
    oni, pna = d.oni, d.pna
    G, n_per = d.G, d.n_per
    res = dict(analysis=name, window=d.window, seasons=f"{d.s0}-{d.s1}", timing="now" if d.now else "lag",
               n_seasons=G, n_days=len(d.C), n_all=int(ys["all"].sum()), n_hf=int(ys["hf"].sum()),
               n_missing_index_days=d.n_missing, outcome_col=hf_c)
    est = point(d, ys, oni, pna)
    # clustered fits for SEs
    X0 = np.column_stack([oni, d.B])
    X1 = np.column_stack([oni, pna, d.B])
    X2 = np.column_stack([oni, pna, oni * pna, d.B])
    for o in ("hf", "all"):
        y = ys[o]
        b0, c0 = pois_fit(y, X0, G, n_per)
        b1, c1 = pois_fit(y, X1, G, n_per)
        b2, c2 = pois_fit(y, X2, G, n_per)
        res.update({f"T_{o}": b0[0], f"seT_{o}": np.sqrt(c0[0, 0]),
                    f"D_{o}": b1[0], f"seD_{o}": np.sqrt(c1[0, 0]),
                    f"b_{o}": b1[1], f"seb_{o}": np.sqrt(c1[1, 1]),
                    f"gamma_{o}": b2[2], f"segamma_{o}": np.sqrt(c2[2, 2])})
    ba, ca = ols_fit(pna, X0, G, n_per)
    res.update(a=ba[0], sea=np.sqrt(ca[0, 0]))
    # shares: stacked fits
    bs, cs = stacked_fit(ys["all"], ys["hf"], X2, G, n_per)
    k = X2.shape[1]
    res.update(gamma_share=bs[k + 2], segamma_share=np.sqrt(cs[k + 2, k + 2]))
    bs1, cs1 = stacked_fit(ys["all"], ys["hf"], X1, G, n_per)
    k1 = X1.shape[1]
    res.update(b_share=bs1[k1 + 1], seb_share=np.sqrt(cs1[k1 + 1, k1 + 1]))
    # point estimates and permutation / bootstrap targets
    obs = {o: est[o] for o in ("hf", "all")}
    res["ind_hf"] = obs["hf"]["ind"]
    res["ind_all"] = obs["all"]["ind"]
    res["ind_share"] = res["ind_hf"] - res["ind_all"]
    res["T_share"] = res["T_hf"] - res["T_all"]
    res["diffmethod_hf"] = res["T_hf"] - res["D_hf"]
    # ---- permutations
    pn = pna.reshape(G, n_per)
    on = oni.reshape(G, n_per)
    ge = dict(gamma_hf=0, gamma_all=0, gamma_share=0, ind_hf=0, ind_all=0, ind_share=0)
    zobs = dict(gamma_hf=abs(res["gamma_hf"] / res["segamma_hf"]), gamma_all=abs(res["gamma_all"] / res["segamma_all"]))
    for _ in range(nperm):
        p1 = rng.permutation(G)  # P1: both indices move together
        o1, q1 = on[p1].reshape(-1), pn[p1].reshape(-1)
        wz = {}
        gm = {}
        for o in ("hf", "all"):
            g_, s_ = wald_gamma(d, ys[o], o1, q1)
            wz[o] = abs(g_ / s_)
            gm[o] = g_
        ge["gamma_hf"] += wz["hf"] >= zobs["gamma_hf"] - 1e-12
        ge["gamma_all"] += wz["all"] >= zobs["gamma_all"] - 1e-12
        ge["gamma_share"] += abs(gm["hf"] - gm["all"]) >= abs(res["gamma_hf"] - res["gamma_all"]) - 1e-12
        p2 = rng.permutation(G)  # P2: ONI blocks move against PNA and counts
        o2 = on[p2].reshape(-1)
        e2 = point(d, ys, o2, pna)
        ge["ind_hf"] += abs(e2["hf"]["ind"]) >= abs(res["ind_hf"]) - 1e-12
        ge["ind_all"] += abs(e2["all"]["ind"]) >= abs(res["ind_all"]) - 1e-12
        ge["ind_share"] += abs(e2["hf"]["ind"] - e2["all"]["ind"]) >= abs(res["ind_share"]) - 1e-12
    for kx, v in ge.items():
        res[f"pperm_{kx}"] = (v + 1) / (nperm + 1)
    # ---- season-block bootstrap
    rows_by_s = [np.arange(i * n_per, (i + 1) * n_per) for i in range(G)]
    keys = ["T_hf", "D_hf", "b_hf", "gamma_hf", "ind_hf", "T_all", "gamma_all", "ind_all", "a"]
    bs_list = []
    for _ in range(nperm):
        pick = rng.integers(0, G, G)
        rows = np.concatenate([rows_by_s[s] for s in pick])
        dd = _Sub(d, rows)
        ysb = {o: ys[o][rows] for o in ys}
        e = point(dd, ysb, oni[rows], pna[rows])
        row = [e["hf"]["T"], e["hf"]["D"], e["hf"]["b"], e["hf"]["gamma"], e["hf"]["ind"],
               e["all"]["T"], e["all"]["gamma"], e["all"]["ind"], e["hf"]["a"]]
        bs_list.append(row)
    bsa = np.array(bs_list)
    for j, kx in enumerate(keys):
        res[f"boot_lo_{kx}"], res[f"boot_hi_{kx}"] = np.percentile(bsa[:, j], [2.5, 97.5])
    sh_g = bsa[:, 3] - bsa[:, 6]
    sh_i = bsa[:, 4] - bsa[:, 7]
    res["boot_lo_gamma_share"], res["boot_hi_gamma_share"] = np.percentile(sh_g, [2.5, 97.5])
    res["boot_lo_ind_share"], res["boot_hi_ind_share"] = np.percentile(sh_i, [2.5, 97.5])
    # share mediated: only if the total effect's bootstrap CI excludes 0
    t_excl = res["boot_lo_T_hf"] > 0 or res["boot_hi_T_hf"] < 0
    res["T_ci_excludes_0"] = bool(t_excl)
    res["share_mediated"] = res["ind_hf"] / res["T_hf"] if t_excl else np.nan
    # ---- leave one season out
    lo = {kx: [] for kx in ("gamma_hf", "ind_hf", "T_hf", "D_hf")}
    for s in range(G):
        rows = np.concatenate([rows_by_s[i] for i in range(G) if i != s])
        dd = _Sub(d, rows)
        e = point(dd, {"hf": ys["hf"][rows]}, oni[rows], pna[rows])["hf"]
        lo["gamma_hf"].append(e["gamma"]); lo["ind_hf"].append(e["ind"])
        lo["T_hf"].append(e["T"]); lo["D_hf"].append(e["D"])
    for kx, v in lo.items():
        res[f"loso_min_{kx}"], res[f"loso_max_{kx}"] = min(v), max(v)
    # ---- detectable effect (80% power, two-sided .05): 2.8 x clustered SE
    res["mde_gamma_hf"] = 2.8 * res["segamma_hf"]
    # SE of the indirect effect: bootstrap SD
    res["se_ind_hf_boot"] = float(np.std(bsa[:, 4], ddof=1))
    res["mde_ind_hf"] = 2.8 * res["se_ind_hf_boot"]
    res["se_T_hf"] = res["seT_hf"]
    res["corr_oni_pna"] = float(np.corrcoef(oni, pna)[0, 1])
    return res


class _Sub:
    """Row-subset view of a Data object (design matrix only)."""
    def __init__(self, d, rows):
        self.B = d.B[rows]


def decision(res, sesoi, key):
    b, lo, hi = res[f"{key}_hf"], res[f"boot_lo_{key}_hf"], res[f"boot_hi_{key}_hf"]
    excl = lo > 0 or hi < 0
    p2 = min(1.0, 2 * res[f"pperm_{key}_hf"])
    rr_lo, rr_hi = np.exp(lo), np.exp(hi)
    if p2 < 0.05 and excl:
        return "detected"
    if rr_lo > 1 / sesoi and rr_hi < sesoi:
        return "well-powered null"
    return "inconclusive"


# ------------------------------------------------------------- special cells
def asymmetry(d, hf_col, nperm, rng):
    """S3: PNA x N+ and PNA x N-. Wald (2 df) of the two product terms, clustered."""
    y = d.C[hf_col].values.astype(float)
    pna, npl, nmi = d.pna, d.nplus, d.nminus
    X = np.column_stack([npl, nmi, pna, npl * pna, nmi * pna, d.B])
    b, cov = pois_fit(y, X, d.G, d.n_per)
    idx = [3, 4]
    W = float(b[idx] @ np.linalg.solve(cov[np.ix_(idx, idx)], b[idx]))
    from scipy import stats
    out = dict(analysis="S3 asymmetry", seasons=f"{d.s0}-{d.s1}",
               g_plus=b[3], se_plus=np.sqrt(cov[3, 3]), g_minus=b[4], se_minus=np.sqrt(cov[4, 4]),
               wald=W, p_asym=float(stats.chi2.sf(W, 2)))
    # block permutation, both indices (N+, N-, PNA) moved together
    G, n = d.G, d.n_per
    A = np.stack([npl, nmi, pna], 1).reshape(G, n, 3)
    ge = 0
    for _ in range(nperm):
        Ap = A[rng.permutation(G)].reshape(-1, 3)
        Xp = np.column_stack([Ap[:, 0], Ap[:, 1], Ap[:, 2], Ap[:, 0] * Ap[:, 2], Ap[:, 1] * Ap[:, 2], d.B])
        bp, cp = pois_fit(y, Xp, G, n)
        ge += float(bp[idx] @ np.linalg.solve(cp[np.ix_(idx, idx)], bp[idx])) >= W - 1e-12
    out["pperm_asym"] = (ge + 1) / (nperm + 1)
    return out


def phase_table(d, hf_col):
    """S4: descriptive 3x3 ENSO phase x PNA phase table with additive predictions."""
    import statsmodels.api as sm
    y = d.C[hf_col].values.astype(float)
    en = np.where(np.isnan(d.oni_raw), 0, d.oni_raw)
    ens = np.where(en >= 0.5, 2, np.where(en <= -0.5, 0, 1))  # 0 La Nina, 1 neutral, 2 El Nino
    pn = np.where(d.pna > 0.5, 2, np.where(d.pna < -0.5, 0, 1))  # 0 -PNA, 1 neutral, 2 +PNA
    cell = ens * 3 + pn
    names = {0: "La Nina", 1: "neutral", 2: "El Nino"}
    pnames = {0: "-PNA", 1: "neutral", 2: "+PNA"}
    cols = [c for c in range(9) if c != 4]
    X = np.column_stack([(cell == c).astype(float) for c in cols] + [d.B])
    r = sm.GLM(y, X, family=sm.families.Poisson()).fit(cov_type="cluster", cov_kwds={"groups": d.season})
    beta, cov = r.params, r.cov_params()
    pos = {c: i for i, c in enumerate(cols)}
    rows = []
    for c in range(9):
        e, p = divmod(c, 3)
        ndays = int((cell == c).sum())
        nseas = int(len(np.unique(d.season[cell == c])))
        nhf = int(y[cell == c].sum())
        row = dict(enso=names[e], pna=pnames[p], days=ndays, seasons=nseas, hf=nhf,
                   rate_per_100d=100 * nhf / max(ndays, 1))
        if c != 4:
            row["rr"] = float(np.exp(beta[pos[c]]))
            row["log_se"] = float(np.sqrt(cov[pos[c], pos[c]]))
            if e != 1 and p != 1:
                ce, cp_ = 3 * e + 1, 3 + p  # (ENSO, neutral PNA), (neutral ENSO, PNA)
                v = np.zeros(len(beta)); v[pos[c]] = 1; v[pos[ce]] = -1; v[pos[cp_]] = -1
                lr = float(v @ beta)
                row["additive_pred_rr"] = float(np.exp(beta[pos[ce]] + beta[pos[cp_]]))
                row["excess_ratio"] = float(np.exp(lr))
                row["excess_log_se"] = float(np.sqrt(v @ cov @ v))
        rows.append(row)
    return rows


# ------------------------------------------------------------- main
CPC = REPO = None


def main():
    global CPC, REPO
    tracks, CPC, REPO, out = sys.argv[1:5]
    nperm = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(SEED)
    T = S.load_tracks(tracks)
    I = S.indices(CPC, REPO)
    cut, n_hf, n_cut = S.depth_cut(T, BASIN)
    T["dq"] = (T.minp <= cut).astype(int)
    R, extra = [], {}

    def cell(name, window, s0, s1, col, now=False, era=False):
        d = Data(T, I, window, s0, s1, ["all", col], now=now, era=era)
        r = run_cell(d, name, nperm, rng, col)
        R.append(r)
        print(name, r["seasons"], r["timing"], "done", flush=True)
        return d

    d_main = cell("P primary", "octapr", 2004, 2025, "hf")
    extra["S3"] = asymmetry(d_main, "hf", nperm, rng)
    extra["S4"] = phase_table(d_main, "hf")
    cell("S1 same-time", "octapr", 2004, 2025, "hf", now=True)
    cell("S2 Jun-May", "junmay", 2004, 2025, "hf")
    for (s0, s1) in ((1979, 2025), (1979, 2000), (2001, 2025)):
        cell("S5 depth", "octapr", s0, s1, "dq", era=(s0 == 1979 and s1 == 2025))
    cell("S6 gust 2001-25", "octapr", 2001, 2025, "hf")
    cell("S6 gust 1979-2025 era term", "octapr", 1979, 2025, "hf", era=True)
    cell("S6 gust 1979-2000 (within-era)", "octapr", 1979, 2000, "hf")
    # S3 also on the depth version and on 2001-2025? Not pre-registered: not run.
    meta = dict(depth_cut_pac=cut, n_hf_2004_2025_pac=n_hf, n_depth_2004_2025_pac=n_cut,
                nperm=nperm, nboot=nperm, seed=SEED, n_tracks_file=int(len(T)))
    df = pd.DataFrame(R)
    # Benjamini-Hochberg over every gamma and indirect permutation p (HF, all, share) of every cell
    ps, tags = [], []
    for i, r in df.iterrows():
        for kx in ("gamma_hf", "gamma_all", "gamma_share", "ind_hf", "ind_all", "ind_share"):
            ps.append(r[f"pperm_{kx}"]); tags.append((i, kx))
    ps.append(extra["S3"]["pperm_asym"]); tags.append(("S3", "asym"))
    q = bh(np.array(ps))
    meta["fdr_family_size"] = len(ps)
    for (i, kx), qq in zip(tags, q):
        if i == "S3":
            extra["S3"]["q_asym"] = float(qq)
        else:
            df.loc[i, f"q_{kx}"] = qq
    prim = df.iloc[0]
    meta["decision_P1_gamma"] = decision(prim, 1.05, "gamma")
    meta["decision_P2_ind"] = decision(prim, 1.03, "ind")
    df.to_csv(os.path.join(out, "results.csv"), index=False, float_format="%.5g")
    json.dump(dict(meta=meta, **extra), open(os.path.join(out, "extra.json"), "w"), indent=1, default=float)
    write_summary(df, meta, extra, os.path.join(out, "summary.txt"))


def bh(p):
    o = np.argsort(p)
    q = p[o] * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    r = np.empty_like(q)
    r[o] = np.minimum(q, 1)
    return r


def rr(b):
    return f"{np.exp(b):.3f}"


def write_summary(df, meta, extra, path):
    L = ["ERA5 PROXY, pipeline A (800 km ocean gust index; HF-equivalent >= 71.7 kt). Pacific basin.",
         "Daily genesis counts, Poisson, month FE + season trend (+ era term where stated); SEs clustered by season;",
         "p from season-block permutation; [lo, hi] season-block bootstrap 95%; LOSO = leave-one-season-out range.",
         "RR scale = exp(coefficient); gamma per SD x SD of the two standardised, lagged indices.",
         f"permutations = bootstrap draws = {meta['nperm']}; BH family size {meta['fdr_family_size']}.",
         f"Pre-registered decisions (primary window): P1 gamma {meta['decision_P1_gamma']}; "
         f"P2 indirect effect {meta['decision_P2_ind']}.", ""]
    for _, r in df.iterrows():
        L.append(f"[{r.analysis}] {r.window} {r.seasons} index {r.timing}; outcome column {r.outcome_col}")
        L.append(f"  seasons {r.n_seasons}, days {r.n_days}, cyclones {r.n_all}, HF-like {r.n_hf}, "
                 f"index days filled {r.n_missing_index_days}, corr(ONI,PNA) over days {r.corr_oni_pna:+.2f}")
        L.append(f"  P1 interaction gamma (ONI x PNA)")
        for o, lab in (("hf", "HF     "), ("all", "all    "), ("share", "share  ")):
            lo_, hi_ = r[f"boot_lo_gamma_{o}"], r[f"boot_hi_gamma_{o}"]
            L.append(f"    {lab} RR {rr(r[f'gamma_{o}'])} (log {r[f'gamma_{o}']:+.4f}, se {r[f'segamma_{o}']:.4f}, "
                     f"z {r[f'gamma_{o}'] / r[f'segamma_{o}']:+.2f}, p_perm {r[f'pperm_gamma_{o}']:.4f}, q {r[f'q_gamma_{o}']:.3f}) "
                     f"boot [{np.exp(lo_):.3f}, {np.exp(hi_):.3f}]")
        L.append(f"    LOSO (HF) [{rr(r.loso_min_gamma_hf)}, {rr(r.loso_max_gamma_hf)}]; detectable at 80% power: RR {rr(r.mde_gamma_hf)}")
        L.append(f"  P2 pathway (HF)")
        L.append(f"    total T (ONI alone) RR {rr(r.T_hf)} (se {r.seT_hf:.4f}) boot [{rr(r.boot_lo_T_hf)}, {rr(r.boot_hi_T_hf)}]")
        L.append(f"    direct D (PNA held) RR {rr(r.D_hf)} (se {r.seD_hf:.4f}) boot [{rr(r.boot_lo_D_hf)}, {rr(r.boot_hi_D_hf)}]")
        L.append(f"    path a (SD PNA per SD ONI) {r.a:+.3f} (se {r.sea:.3f}) boot [{r.boot_lo_a:+.3f}, {r.boot_hi_a:+.3f}]")
        L.append(f"    PNA effect b (ONI held) RR {rr(r.b_hf)} (se {r.seb_hf:.4f})")
        L.append(f"    indirect a x b RR {rr(r.ind_hf)} (log {r.ind_hf:+.4f}) boot [{rr(r.boot_lo_ind_hf)}, {rr(r.boot_hi_ind_hf)}] "
                 f"p_perm {r.pperm_ind_hf:.4f} q {r.q_ind_hf:.3f}; T - D = {r.diffmethod_hf:+.4f}")
        sm = "undefined (total effect CI includes 0)" if not r.T_ci_excludes_0 else f"{r.share_mediated:.2f}"
        L.append(f"    share mediated {sm}; LOSO ind [{rr(r.loso_min_ind_hf)}, {rr(r.loso_max_ind_hf)}]; detectable (80%) RR {rr(r.mde_ind_hf)}")
        L.append(f"    all cyclones: T {rr(r.T_all)}, indirect {rr(r.ind_all)} (p_perm {r.pperm_ind_all:.3f}); "
                 f"share: indirect {rr(r.ind_share)} boot [{rr(r.boot_lo_ind_share)}, {rr(r.boot_hi_ind_share)}] p_perm {r.pperm_ind_share:.3f}")
        L.append("")
    s3 = extra["S3"]
    L.append(f"[S3 asymmetry] 2004-2025: PNA x El Nino-side (N+) RR {rr(s3['g_plus'])} (se {s3['se_plus']:.4f}), "
             f"PNA x La Nina-side (N-) RR {rr(s3['g_minus'])} (se {s3['se_minus']:.4f}); joint Wald {s3['wald']:.2f}, "
             f"p {s3['p_asym']:.3f}, p_perm {s3['pperm_asym']:.4f}, q {s3['q_asym']:.3f}")
    L.append("")
    L.append("[S4 phase table, descriptive] HF genesis rate per 100 days; RR vs neutral/neutral; additive prediction = RR(ENSO, neutral PNA) x RR(neutral ENSO, PNA)")
    L.append(f"  {'ENSO':8s} {'PNA':8s} {'days':>5s} {'seas':>4s} {'HF':>4s} {'rate':>6s} {'RR':>6s} {'additive':>8s} {'excess':>7s} {'z':>6s}")
    for r in extra["S4"]:
        ex = f"{r['excess_ratio']:7.3f} {np.log(r['excess_ratio']) / r['excess_log_se']:+6.2f}" if "excess_ratio" in r else ""
        ad = f"{r['additive_pred_rr']:8.3f}" if "additive_pred_rr" in r else f"{'':8s}"
        L.append(f"  {r['enso']:8s} {r['pna']:8s} {r['days']:5d} {r['seasons']:4d} {r['hf']:4d} {r['rate_per_100d']:6.2f} "
                 f"{r.get('rr', 1.0):6.3f} {ad} {ex}")
    L.append("")
    L.append("Depth cut (Pacific): " + json.dumps({k: v for k, v in meta.items() if "depth" in k or k == "n_hf_2004_2025_pac"}))
    open(path, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
