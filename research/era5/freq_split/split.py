"""Split a teleconnection's effect on HF-low frequency into more cyclones and a
larger share of cyclones reaching hurricane force. ERA5 PROXY, pipeline A.

    RR(HF lows) = RR(all cyclones) x RR(share reaching HF)

Plan and hypotheses: PREREGISTRATION.md (committed before any fit).

Cyclones: every pipeline A track in research/era5/hf_history/results/all_tracks.csv.gz
(lows below 1010 hPa, >= 24 h, >= 2 fixes in a basin domain). HF-equivalent:
800 km ocean gust index >= 71.7 kt. Each track sits on its genesis day (first fix).

Index: CPC daily NAO / PNA. Lagged = mean of days -10..-4 before genesis (as Q1);
same-time = mean of days -3..+3. Standardised over the analysis days.

Counts: daily genesis counts per basin, Poisson GLM with log link, index + calendar
month fixed effects + linear season trend. All and HF counts are fitted as one
stacked model with every term interacted with an HF indicator, so the
index x HF coefficient is log RR(share) = b(HF) - b(all) exactly, and its
season-clustered SE accounts for HF being a subset of all. The HF-share odds ratio
comes from a track-level logistic with the same terms.

Inference: season-clustered SEs; season-block permutation (each season's index
series moved to another season, day-of-season kept); season-block bootstrap for
intervals on RR(all), RR(share) and f = log RR(share) / log RR(HF); leave-one-season-out.

usage: split.py ALL_TRACKS_CSV_GZ CPC_DIR REPO_ROOT OUT_DIR [NPERM]
CPC_DIR holds norm.daily.{nao,pna}.index.b500101.current.ascii and oni_*.txt;
REPO_ROOT supplies docs/data/teleconnections.json for ONI after 2002.
"""
import datetime as dt
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore")
THR = 71.7
DEEP_MINP, DEEP_NFIX = 1000.0, 8
SEED = 20261008
ONI_SEAS = ["DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ"]
WINDOWS = {"octapr": ((10, 1), 210), "junmay": ((6, 1), 364)}


# ---------------------------------------------------------------- inputs
def read_cpc(path):
    """CPC daily file, fixed width 'YYYY MM DD value'; -99.000 is missing and can
    run into the day field ('2006 10 26-99.000')."""
    vals = {}
    for ln in open(path):
        if len(ln) < 11 or not ln[:4].isdigit():
            continue
        v = float(ln[10:])
        vals[pd.Timestamp(int(ln[:4]), int(ln[5:7]), int(ln[8:10]))] = np.nan if v <= -98 else v
    return pd.Series(vals).sort_index().asfreq("D")


def read_oni(cpc, repo):
    oni = {}
    t = json.load(open(os.path.join(repo, "docs/data/teleconnections.json")))["oni"]
    y0, m0 = map(int, t["start"].split("-"))
    for i, v in enumerate(t["values"]):
        if v is not None:
            oni[(y0 + (m0 - 1 + i) // 12, (m0 - 1 + i) % 12 + 1)] = v
    for fn in sorted(os.listdir(cpc)):
        if fn.startswith("oni"):
            for ln in open(os.path.join(cpc, fn)):
                p = ln.split()
                if len(p) == 4 and p[0] in ONI_SEAS:
                    oni.setdefault((int(p[1]), ONI_SEAS.index(p[0]) + 1), float(p[3]))
    return oni


def indices(cpc, repo):
    """Daily series: <K>_lag (days -10..-4) and <K>_now (days -3..+3)."""
    out = {}
    for k in ("NAO", "PNA"):
        s = read_cpc(os.path.join(cpc, f"norm.daily.{k.lower()}.index.b500101.current.ascii"))
        out[k + "_lag"] = s.rolling(7, min_periods=7).mean().shift(4)
        out[k + "_now"] = s.rolling(7, min_periods=7, center=True).mean()
    oni = read_oni(cpc, repo)
    days = out["NAO_lag"].index
    lagday = days - pd.Timedelta(days=7)  # month of the lag window's middle day
    out["ONI_lag"] = pd.Series([oni.get((d.year, d.month), np.nan) for d in lagday], index=days)
    return pd.DataFrame(out)


def load_tracks(path):
    T = pd.read_csv(path, dtype={"start": str, "peak_time": str})
    T["gen"] = pd.to_datetime(T.start.str[:8], format="%Y%m%d")
    T["hf"] = (T.gust800_kt >= THR).astype(int)
    T["deep"] = ((T.minp <= DEEP_MINP) & (T.n_fix >= DEEP_NFIX)).astype(int)
    return T


# ---------------------------------------------------------------- tables
def window_days(window, s0, s1):
    (m, d), n = WINDOWS[window]
    rows = []
    for s in range(s0, s1 + 1):
        d0 = dt.date(s, m, d)
        for k in range(n):
            rows.append((s, k, pd.Timestamp(d0 + dt.timedelta(days=k))))
    D = pd.DataFrame(rows, columns=["season", "day", "date"])
    D["month"] = D.date.dt.month
    return D


def assign(T, window, s0, s1):
    (m, d), n = WINDOWS[window]
    g = T.gen
    season = np.where((g.dt.month > m) | ((g.dt.month == m) & (g.dt.day >= d)), g.dt.year, g.dt.year - 1)
    start = pd.to_datetime(dict(year=season, month=m, day=d))
    day = (g - start).dt.days.values
    T = T.assign(season=season, day=day)
    return T[(T.season >= s0) & (T.season <= s1) & (T.day >= 0) & (T.day < n)].copy()


def daily_counts(T, D, basin, outcomes):
    x = T[T.basin == basin]
    C = D[["season", "day", "month"]].copy()
    key = pd.MultiIndex.from_frame(C[["season", "day"]])
    for o in outcomes:
        c = x[x[o] == 1].groupby(["season", "day"]).size() if o != "all" else x.groupby(["season", "day"]).size()
        C[o] = c.reindex(key, fill_value=0).values
    return C


def design(C, z, idx_cols):
    """Index columns + month dummies + centred season trend + intercept."""
    X = pd.DataFrame({c: z[c] for c in idx_cols})
    for mo in sorted(C.month.unique())[1:]:
        X[f"m{mo}"] = (C.month.values == mo).astype(float)
    X["trend"] = (C.season.values - C.season.mean()) / 10.0
    X["const"] = 1.0
    return X


# ---------------------------------------------------------------- fits
def pois(y, X):
    return sm.GLM(y, X, family=sm.families.Poisson()).fit()


def stacked(C, X, a, b):
    """Stack outcome a (denominator) and b; every column interacted with is_b.
    Returns (coef of index_k for a, coef of index_k x is_b) and their clustered SEs."""
    n = len(C)
    Xa = X.values
    Z = np.zeros_like(Xa)
    XX = np.vstack([np.hstack([Xa, Z]), np.hstack([Xa, Xa])])
    y = np.concatenate([C[a].values, C[b].values])
    g = np.concatenate([C.season.values, C.season.values])
    r = sm.GLM(y, XX, family=sm.families.Poisson()).fit(cov_type="cluster", cov_kwds={"groups": g})
    return r, Xa.shape[1]


def fast_pois(y, X, iters=30):
    """Newton-Raphson Poisson MLE (same estimate as the GLM; used in resampling loops)."""
    b = np.zeros(X.shape[1])
    b[-1] = np.log(max(y.mean(), 1e-9))
    for _ in range(iters):
        mu = np.exp(X @ b)
        step = np.linalg.solve(X.T @ (X * mu[:, None]), X.T @ (y - mu))
        b += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b


def coefs(C, X, chain, k):
    """log RR per SD of index column k for every outcome in chain, plain Poisson fits."""
    Xv = np.asarray(X, float)
    return {o: fast_pois(np.asarray(C[o], float), Xv)[k] for o in chain}


def logit_share(T, Z, idx):
    """Track-level logistic HF ~ index (+ any joint index) + month FE + trend, season-clustered."""
    y = T.hf.values
    X = pd.DataFrame(Z, columns=list(Z.columns) if hasattr(Z, "columns") else [idx], index=T.index)
    for mo in sorted(T.gen.dt.month.unique())[1:]:
        X[f"m{mo}"] = (T.gen.dt.month.values == mo).astype(float)
    X["trend"] = (T.season.values - T.season.mean()) / 10.0
    X["const"] = 1.0
    r = sm.GLM(y, X, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": T.season.values})
    return r.params[idx], r.bse[idx]


# ---------------------------------------------------------------- one decomposition
def decompose(T, I, basin, idx, window, s0, s1, nperm, rng, chain=("all", "hf"), now=False, joint=None):
    """chain: outcomes from widest to narrowest; the last is HF (or a depth stand-in)."""
    D = window_days(window, s0, s1)
    Tw = assign(T, window, s0, s1)
    Tw = Tw[Tw.basin == basin]
    col = idx + ("_now" if now else "_lag")
    cols = [col] + ([j + ("_now" if now else "_lag") for j in joint] if joint else [])
    raw = I.reindex(D.date)[cols].reset_index(drop=True)
    ok = raw.notna().all(axis=1).values
    D, raw = D[ok].reset_index(drop=True), raw[ok].reset_index(drop=True)
    mu, sd = raw.mean(), raw.std(ddof=0)
    z = (raw - mu) / sd
    C = daily_counts(Tw, D, basin, chain)
    X = design(C, z, cols)
    k = 0
    top, last = chain[0], chain[-1]

    # point estimates and clustered SEs
    res = dict(basin=basin, idx=idx, window=window, seasons=f"{s0}-{s1}", timing="now" if now else "lag",
               joint="+".join([idx] + list(joint)) if joint else "", n_seasons=s1 - s0 + 1,
               n_days=len(C), n_tracks=int(C[top].sum()), n_hf=int(C[last].sum()),
               share=C[last].sum() / C[top].sum(), sd_index=float(sd.iloc[0]), outcome=last)
    b = coefs(C, X, chain, k)
    for o in chain:
        r = pois(C[o].values, X)
        rc = sm.GLM(C[o].values, X, family=sm.families.Poisson()).fit(cov_type="cluster", cov_kwds={"groups": C.season.values})
        res[f"b_{o}"], res[f"se_{o}"] = rc.params.iloc[k], rc.bse.iloc[k]
        res[f"disp_{o}"] = float(r.pearson_chi2 / r.df_resid)
    for a, bb in zip(chain[:-1], chain[1:]):
        r, p = stacked(C, X, a, bb)
        res[f"b_{bb}|{a}"], res[f"se_{bb}|{a}"] = r.params[p + k], r.bse[p + k]
    r, p = stacked(C, X, top, last)
    res["b_share"], res["se_share"] = r.params[p + k], r.bse[p + k]
    res["f_share"] = res["b_share"] / res[f"b_{last}"]
    # odds ratio, track level (only for the gust HF outcome on the full population)
    if chain == ("all", "hf") or chain == ("all", "deep", "hf"):
        zt = (I.reindex(Tw.gen)[cols].reset_index(drop=True) - mu) / sd
        okt = zt.notna().all(axis=1).values
        lb, lse = logit_share(Tw[okt], zt[okt].set_index(Tw.index[okt]), col)
        res["logOR_share"], res["se_logOR_share"] = lb, lse

    # season-block permutation: move each season's index block to another season
    seasons = np.arange(s0, s1 + 1)
    zz = {s: z[C.season.values == s].values for s in seasons}
    lens = {s: len(v) for s, v in zz.items()}
    obs = np.array([b[o] for o in chain] + [b[last] - b[top]])
    ge = np.zeros(len(obs))
    nvalid = 0
    for _ in range(nperm):
        perm = rng.permutation(seasons)
        blocks = []
        good = True
        for s, ps in zip(seasons, perm):
            v = zz[ps]
            if len(v) != lens[s]:  # different valid-day count (leap years, gaps): wrap
                v = np.resize(v, (lens[s], v.shape[1]))
            blocks.append(v)
        Zp = np.vstack(blocks)
        Xp = X.copy()
        Xp[cols] = Zp
        bp = coefs(C, Xp, chain, k)
        stat = np.array([bp[o] for o in chain] + [bp[last] - bp[top]])
        ge += np.abs(stat) >= np.abs(obs) - 1e-12
        nvalid += 1
    pnames = list(chain) + ["share"]
    for nm, g_ in zip(pnames, ge):
        res[f"pperm_{nm}"] = (g_ + 1) / (nvalid + 1)

    # season-block bootstrap
    nb = nperm
    bs = []
    idx_by_s = {s: np.where(C.season.values == s)[0] for s in seasons}
    for _ in range(nb):
        pick = rng.choice(seasons, len(seasons), replace=True)
        rows = np.concatenate([idx_by_s[s] for s in pick])
        Cb, Xb = C.iloc[rows].reset_index(drop=True), X.iloc[rows].reset_index(drop=True)
        try:
            bb = coefs(Cb, Xb, (top, last), k)
        except Exception:
            continue
        bs.append((bb[top], bb[last] - bb[top], bb[last]))
    bs = np.array(bs)
    fb = bs[:, 1] / bs[:, 2]
    for j, nm in enumerate(("all", "share", "hf")):
        lo, hi = np.percentile(bs[:, j], [2.5, 97.5])
        res[f"boot_lo_{nm}"], res[f"boot_hi_{nm}"] = lo, hi
    res["f_boot_lo"], res["f_boot_hi"] = np.percentile(fb, [2.5, 97.5])
    res["f_boot_median"] = float(np.median(fb))
    res["n_boot"] = len(bs)

    # leave one season out
    lo_ = {o: [] for o in ("top", "share", "last", "f")}
    for s in seasons:
        keep = C.season.values != s
        bl = coefs(C[keep], X[keep], (top, last), k)
        lo_["top"].append(bl[top]); lo_["last"].append(bl[last]); lo_["share"].append(bl[last] - bl[top])
        lo_["f"].append((bl[last] - bl[top]) / bl[last])
    for o, v in lo_.items():
        res[f"loso_min_{o}"], res[f"loso_max_{o}"] = min(v), max(v)
    return res


def depth_cut(T, basin):
    """minp cut matching the 71.7 kt count, primary window and seasons."""
    x = assign(T, "octapr", 2004, 2025)
    x = x[x.basin == basin]
    n = int(x.hf.sum())
    cut = np.sort(x.minp.values)[n - 1]
    return float(cut), n, int((x.minp <= cut).sum())


# ---------------------------------------------------------------- main
def fmt_rr(b):
    return f"{np.exp(b):.3f}"


def main():
    tracks, cpc, repo, out = sys.argv[1:5]
    nperm = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(SEED)
    T = load_tracks(tracks)
    I = indices(cpc, repo)
    T["deepor"] = ((T.deep == 1) | (T.hf == 1)).astype(int)
    R = []
    tag = lambda r, **kw: R.append({**r, **kw})

    prim = [("atl", "NAO"), ("pac", "PNA")]
    # primary
    for basin, idx in prim:
        tag(decompose(T, I, basin, idx, "octapr", 2004, 2025, nperm, rng), analysis="primary")
    # S1 same-time index
    for basin, idx in prim:
        tag(decompose(T, I, basin, idx, "octapr", 2004, 2025, nperm, rng, now=True), analysis="S1 same-time")
    # S2 June-May
    for basin, idx in prim:
        tag(decompose(T, I, basin, idx, "junmay", 2004, 2025, nperm, rng), analysis="S2 Jun-May")
    # S3 three-stage chain (deep stage = minp <= 1000 and >= 48 h, or HF)
    for basin, idx in prim:
        tag(decompose(T, I, basin, idx, "octapr", 2004, 2025, nperm, rng, chain=("all", "deepor", "hf")),
            analysis="S3 chain")
    # S4 depth-equivalent, 1979+
    cuts = {}
    for basin, idx in prim:
        cut, n_hf, n_cut = depth_cut(T, basin)
        cuts[basin] = dict(minp_cut=cut, n_hf_2004_2025=n_hf, n_depth_2004_2025=n_cut)
        T[f"dq_{basin}"] = (T.minp <= cut).astype(int)
        for s0, s1 in ((1979, 2025), (1979, 2000), (2004, 2025)):
            tag(decompose(T, I, basin, idx, "octapr", s0, s1, nperm, rng, chain=("all", f"dq_{basin}")),
                analysis="S4 depth", depth_cut=cut)
    # S5 cross-basin index and ONI
    for basin, idx in (("atl", "PNA"), ("pac", "NAO"), ("atl", "ONI"), ("pac", "ONI")):
        tag(decompose(T, I, basin, idx, "octapr", 2004, 2025, nperm, rng), analysis="S5 other index")
    # S6 NAO and PNA together
    for basin, idx, other in (("atl", "NAO", "PNA"), ("pac", "PNA", "NAO")):
        tag(decompose(T, I, basin, idx, "octapr", 2004, 2025, nperm, rng, joint=[other]), analysis="S6 joint")

    df = pd.DataFrame(R)
    df.to_csv(os.path.join(out, "results.csv"), index=False, float_format="%.5g")
    meta = dict(threshold_kt=THR, deep_stage=f"minp <= {DEEP_MINP} hPa and n_fix >= {DEEP_NFIX}, or HF",
                depth_cuts=cuts, nperm=nperm, nboot=nperm, seed=SEED, n_tracks_file=int(len(T)),
                lag="mean of days -10..-4 before genesis", now="mean of days -3..+3")
    json.dump(meta, open(os.path.join(out, "meta.json"), "w"), indent=1)
    write_summary(df, meta, os.path.join(out, "summary.txt"))


def write_summary(df, meta, path):
    L = ["ERA5 PROXY, pipeline A (800 km ocean gust index; HF-equivalent >= 71.7 kt).",
         "RR(HF) = RR(all cyclones) x RR(share). Per SD of the index. Daily genesis counts,",
         "Poisson with month FE + season trend; SEs clustered by season; p from season-block",
         "permutation; [lo, hi] season-block bootstrap 95%; LOSO = leave-one-season-out range.",
         "f = log RR(share) / log RR(HF): the fraction of the effect carried by the share.",
         f"permutations = bootstrap draws = {meta['nperm']}.", ""]
    for _, r in df.iterrows():
        last = r.outcome
        L.append(f"[{r.analysis}] {r.basin} {r.idx}{(' (with ' + r.joint + ')') if isinstance(r.joint, str) and r.joint else ''} "
                 f"{r.window} {r.seasons} index {r.timing}")
        L.append(f"  tracks {r.n_tracks}, outcome '{last}' {r.n_hf} (share {r.share:.4f}), days {r.n_days}, "
                 f"dispersion all {r.disp_all:.2f}")
        for nm, b, se, p, lo, hi, l0, l1 in (
                ("RR(all)  ", r.b_all, r.se_all, r.pperm_all, r.boot_lo_all, r.boot_hi_all, r.loso_min_top, r.loso_max_top),
                ("RR(share)", r.b_share, r.se_share, r.pperm_share, r.boot_lo_share, r.boot_hi_share, r.loso_min_share, r.loso_max_share),
                ("RR(HF)   ", r[f"b_{last}"], r[f"se_{last}"], r[f"pperm_{last}"], r.boot_lo_hf, r.boot_hi_hf, r.loso_min_last, r.loso_max_last)):
            L.append(f"  {nm} {np.exp(b):.3f}  (log {b:+.4f}, se {se:.4f}, z {b / se:+.2f}, p_perm {p:.4f})  "
                     f"boot [{np.exp(lo):.3f}, {np.exp(hi):.3f}]  LOSO [{np.exp(l0):.3f}, {np.exp(l1):.3f}]")
        L.append(f"  f(share) {r.f_share:.2f}  boot [{r.f_boot_lo:.2f}, {r.f_boot_hi:.2f}]  LOSO [{r.loso_min_f:.2f}, {r.loso_max_f:.2f}]")
        if not pd.isna(r.get("logOR_share", np.nan)):
            L.append(f"  OR(share, track logistic) {np.exp(r.logOR_share):.3f} (se of log {r.se_logOR_share:.4f}, z {r.logOR_share / r.se_logOR_share:+.2f})")
        if r.analysis == "S3 chain":
            L.append(f"  stage RR: deep|all {np.exp(r['b_deepor|all']):.3f} (z {r['b_deepor|all'] / r['se_deepor|all']:+.2f}),"
                     f" HF|deep {np.exp(r['b_hf|deepor']):.3f} (z {r['b_hf|deepor'] / r['se_hf|deepor']:+.2f})")
        L.append("")
    L.append("Depth-equivalent cuts (S4): " + json.dumps(meta["depth_cuts"]))
    open(path, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
