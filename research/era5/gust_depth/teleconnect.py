"""H3: teleconnection response of Pacific (and Atlantic) HF-low counts under G, D, D' and the archive.  ERA5 PROXY, pipeline A.
Plan: PREREGISTRATION.md (committed before this was run).  Machinery reused from ../freq_split/split.py.

usage: teleconnect.py ALL_TRACKS_CSV_GZ ARCHIVE_JSON CPC_DIR REPO_ROOT OUT_DIR [NDRAW]
"""
import os, sys, json, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "freq_split"))
import match, split
warnings.filterwarnings("ignore")
SEED = 20261010


def clim_cells(T):
    """cell key for D': basin, calendar month of peak time, 10-deg lat band, 30-deg lon band."""
    return list(zip(T.basin, T.peak_time.str[4:6], (T.peak_lat // 10).astype(int), (T.peak_lon // 30).astype(int)))


def add_outcomes(T, M):
    T = T.copy()
    T["hf"] = (T.gust800_kt >= split.THR).astype(int)
    T["psea"] = np.where(T.gen.dt.month >= 6, T.gen.dt.year, T.gen.dt.year - 1)
    ref = T[T.psea.between(2004, 2025)]
    key = pd.Series(clim_cells(T), index=T.index)
    kref = key[ref.index]
    cm = ref.groupby(kref).minp.agg(["mean", "size"])
    bm = ref.groupby([ref.basin, ref.peak_time.str[4:6]]).minp.mean()
    cm = cm[cm["size"] >= 25]["mean"]
    clim = key.map(cm)
    fb = pd.Series([bm.get((b, m), np.nan) for b, m in zip(T.basin, T.peak_time.str[4:6])], index=T.index)
    T["clim"] = clim.fillna(fb)
    T["dev"] = T.minp - T.clim
    cuts = {}
    for basin in ("pac", "atl"):
        x = split.assign(T, "octapr", 2004, 2025); x = x[x.basin == basin]
        n = int(x.hf.sum())
        c = split.depth_cut(T, basin)[0]
        cp = float(np.sort(x.dev.values)[n - 1])
        cuts[basin] = dict(depth_cut=c, dev_cut=cp, n_g=n)
        b = T.basin == basin
        T.loc[b, "D"] = (T.loc[b, "minp"] <= c).astype(int)
        T.loc[b, "Dp"] = (T.loc[b, "dev"] <= cp).astype(int)
    T["D"] = T.D.astype(int); T["Dp"] = T.Dp.astype(int)
    T["G"] = T.hf
    T["both"] = ((T.G == 1) & (T.D == 1)).astype(int)
    T["Gonly"] = ((T.G == 1) & (T.D == 0)).astype(int)
    T["Donly"] = ((T.G == 0) & (T.D == 1)).astype(int)
    # archive events on the track's genesis day; unmatched ones: archive start minus the median offset
    mm = M[M.cls == "low"].copy()
    gen = T.set_index("track").gen
    ev = {e["basin"] + ":" + e["id"]: e for e in match.load_archive(os.path.join(HERE, "..", "..", "..", "docs", "data", "hf-lows.json"))}
    rows, offs = [], []
    for r in mm.itertuples():
        e = ev[r.key]
        sd = pd.to_datetime(str(e["start"])[:8], format="%Y%m%d")
        if r.track is not None and not np.isnan(r.track):
            g = gen[int(r.track)]
            offs.append((sd - g).days); rows.append((r.basin, g, int(r.track), True))
        else:
            rows.append((r.basin, sd, None, False))
    off = float(np.median(offs))
    A = []
    for basin, d, tr, m in rows:
        A.append((basin, d if m else d - pd.Timedelta(days=off), m))
    A = pd.DataFrame(A, columns=["basin", "gen", "matched"])
    return T, A, cuts, off


def counts(T, A, D, basin, outs):
    x = T[T.basin == basin]
    C = D[["season", "day", "month"]].copy()
    key = pd.MultiIndex.from_frame(C[["season", "day"]])
    for o in outs:
        c = (x.groupby(["season", "day"]).size() if o == "all" else x[x[o] == 1].groupby(["season", "day"]).size())
        C[o] = c.reindex(key, fill_value=0).values
    return C


def arch_counts(A, D, basin, window, s0, s1):
    a = A[A.basin == basin].reset_index(drop=True)
    (m, d), n = split.WINDOWS[window]
    g = a.gen
    season = np.where((g.dt.month > m) | ((g.dt.month == m) & (g.dt.day >= d)), g.dt.year, g.dt.year - 1)
    start = pd.to_datetime(dict(year=season, month=m, day=d))
    day = (g - start).dt.days.values
    a = a.assign(season=season, day=day)
    a = a[(a.season >= s0) & (a.season <= s1) & (a.day >= 0) & (a.day < n)]
    c = a.groupby(["season", "day"]).size()
    key = pd.MultiIndex.from_frame(D[["season", "day"]])
    return c.reindex(key, fill_value=0).values


def run_basin(T, A, I, basin, idxs, window, s0, s1, ndraw, rng, era=False):
    D = split.window_days(window, s0, s1)
    Tw = split.assign(T, window, s0, s1); Tw = Tw[Tw.basin == basin]
    outs = ["all", "G", "D", "Dp", "both", "Gonly", "Donly"]
    C = counts(Tw, A, D, basin, outs)
    if not era:
        C["arch"] = arch_counts(A, D, basin, window, s0, s1); outs = outs + ["arch"]
    res, diffs = [], []
    seasons = np.arange(s0, s1 + 1)
    for idx in idxs:
        col = idx + "_lag"
        raw = I.reindex(D.date)[[col]].reset_index(drop=True)
        ok = raw.notna().all(axis=1).values
        Dk, rawk = D[ok].reset_index(drop=True), raw[ok].reset_index(drop=True)
        mu, sd = rawk.mean(), rawk.std(ddof=0)
        z = (rawk - mu) / sd
        Ck = C[ok].reset_index(drop=True)
        X = split.design(Ck, z, [col])
        if era:
            X.insert(1, "era", (Ck.season.values <= 2000).astype(float))
        k = 0
        Xv = np.asarray(X, float)
        b0 = {o: split.fast_pois(np.asarray(Ck[o], float), Xv)[k] for o in outs}
        se = {}
        for o in outs:
            r = sm.GLM(Ck[o].values, X, family=sm.families.Poisson()).fit(cov_type="cluster", cov_kwds={"groups": Ck.season.values})
            se[o] = r.bse.iloc[k]
        # permutation: move each season's index block
        zz = {s: z[Ck.season.values == s].values for s in seasons}
        lens = {s: len(v) for s, v in zz.items()}
        ge = {o: 0 for o in outs}
        for _ in range(ndraw):
            perm = rng.permutation(seasons)
            blocks = [np.resize(zz[p], (lens[s], 1)) if len(zz[p]) != lens[s] else zz[p] for s, p in zip(seasons, perm)]
            Xp = Xv.copy(); Xp[:, k] = np.vstack(blocks)[:, 0]
            for o in outs:
                bp = split.fast_pois(np.asarray(Ck[o], float), Xp)[k]
                ge[o] += abs(bp) >= abs(b0[o]) - 1e-12
        # season-block bootstrap, paired across outcomes
        by = {s: np.where(Ck.season.values == s)[0] for s in seasons}
        B = {o: [] for o in outs}
        for _ in range(ndraw):
            pick = rng.choice(seasons, len(seasons), replace=True)
            rows = np.concatenate([by[s] for s in pick])
            Xb = Xv[rows]
            for o in outs:
                try:
                    B[o].append(split.fast_pois(np.asarray(Ck[o], float)[rows], Xb)[k])
                except Exception:
                    B[o].append(np.nan)
        B = {o: np.array(v) for o, v in B.items()}
        for o in outs:
            lo, hi = np.nanpercentile(B[o], [2.5, 97.5])
            res.append(dict(basin=basin, idx=idx, window=window, seasons=f"{s0}-{s1}", era_term=era, outcome=o, n=int(Ck[o].sum()),
                            logRR=b0[o], se=se[o], RR=np.exp(b0[o]), lo=np.exp(lo), hi=np.exp(hi), p_perm=(ge[o] + 1) / (ndraw + 1), sd_index=float(sd.iloc[0])))
        def add(name, a, b, kind="diff"):
            est = (b0[a] - b0[b]) if kind == "diff" else None
            bs = B[a] - B[b]
            ok_ = ~np.isnan(bs)
            lo, hi = np.percentile(bs[ok_], [2.5, 97.5])
            p = 2 * min((bs[ok_] <= 0).mean(), (bs[ok_] >= 0).mean())
            diffs.append(dict(basin=basin, idx=idx, window=window, seasons=f"{s0}-{s1}", era_term=era, contrast=name, est_logRR_diff=est,
                              ratio_of_RR=float(np.exp(est)), lo=float(np.exp(lo)), hi=float(np.exp(hi)), p_boot=float(max(p, 1 / (ok_.sum() + 1)))))
        for a, b in (("D", "G"), ("Dp", "G"), ("D", "Dp")) + ((("arch", "G"), ("arch", "D"), ("arch", "Dp")) if not era else ()):
            add(f"{a} - {b}", a, b)
        # shares: outcome relative to all cyclones
        for o in ("G", "D", "Dp"):
            bs = B[o] - B["all"]; ok_ = ~np.isnan(bs)
            lo, hi = np.percentile(bs[ok_], [2.5, 97.5])
            est = b0[o] - b0["all"]
            p = 2 * min((bs[ok_] <= 0).mean(), (bs[ok_] >= 0).mean())
            res.append(dict(basin=basin, idx=idx, window=window, seasons=f"{s0}-{s1}", era_term=era, outcome=f"share_{o}", n=int(Ck[o].sum()),
                            logRR=est, se=np.nan, RR=np.exp(est), lo=np.exp(lo), hi=np.exp(hi), p_perm=float(max(p, 1 / (ok_.sum() + 1))), sd_index=float(sd.iloc[0])))
        # (a share difference between two outcomes equals their count difference: "all" cancels)
    return pd.DataFrame(res), pd.DataFrame(diffs), C


def main():
    tracks, archive, cpc, repo, out = sys.argv[1:6]
    ndraw = int(sys.argv[6]) if len(sys.argv) > 6 else 2000
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(SEED)
    match.MATCH_KM = 800
    T = split.load_tracks(tracks)
    ev = match.load_archive(archive)
    Tm = pd.read_csv(tracks, dtype={"start": str, "end": str, "peak_time": str})
    M = match.match_all(Tm, ev, 2004)
    T["peak_time"] = T.peak_time.astype(str)
    T, A, cuts, off = add_outcomes(T, M)
    I = split.indices(cpc, repo)
    print("cuts", cuts, "archive offset (days, archive start minus track genesis, median):", off)
    R, Df = [], []
    for basin, idx in (("pac", "PNA"), ("pac", "ONI"), ("atl", "NAO"), ("atl", "ONI")):
        r, d, C = run_basin(T, A, I, basin, [idx], "octapr", 2004, 2025, ndraw, rng)
        R.append(r); Df.append(d); print(basin, idx, "done", flush=True)
    # secondary: 47 seasons with an era term, D and D' only (G is not used across 2001)
    for basin, idx in (("pac", "PNA"), ("pac", "ONI")):
        r, d, C = run_basin(T, A, I, basin, [idx], "octapr", 1979, 2025, ndraw, rng, era=True)
        R.append(r); Df.append(d)
    pd.concat(R).to_csv(os.path.join(out, "tele_estimates.csv"), index=False, float_format="%.5g")
    pd.concat(Df).to_csv(os.path.join(out, "tele_contrasts.csv"), index=False, float_format="%.5g")
    json.dump(dict(cuts=cuts, archive_offset_days=off, ndraw=ndraw, seed=SEED), open(os.path.join(out, "tele_meta.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
