"""RA-22 contrasts: does a teleconnection (or the hemispheric pattern) raise brief HF lows (1 HF fix) as much as sustained ones
(>= 3 HF fixes)? ERA5 PROXY, pipeline A, and the archive. Plan: PREREGISTRATION.md (committed before this was run on k = 2, 3).

Three designs, each reusing the original code's setup unchanged and only changing the outcome:
  daily    PR 14 (freq_split): daily genesis counts, lagged CPC NAO (Atlantic) / PNA (Pacific), month FE + season trend
  weekly   PR 64 (hem_channels): weekly counts by first-fix date, leave-one-season-out pattern index of PR 41
  archive  PR 41 (hemispheric): weekly archive counts (class 'low') by first-fix date, same pattern index

Per design and basin it writes log-scale estimates for HF_k (k = 1, 2, 3), the brief class (exactly 1 HF fix) and the
paired differences (HF_k - HF_1; sustained - brief), each with a season-block bootstrap interval, a season-block permutation
p and a season-clustered SE (for the detectable effect 2.8 x SE).

usage: contrast.py DESIGN OUT_CSV NPERM NBOOT
"""
import datetime as dt, json, os, sys, warnings
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
REPO = os.path.abspath(os.path.join(ERA, "..", ".."))
sys.path.insert(0, os.path.join(ERA, "hem_channels"))
sys.path.insert(0, os.path.join(ERA, "freq_split"))
import chanlib as C       # noqa: E402
import channels as CH     # noqa: E402
import split as S         # noqa: E402

warnings.filterwarnings("ignore")
SEED = 20261008 + 22
CPC = os.environ.get("CPC_DIR", "/mnt/project-files/teleconnection-test/cpc_indices")
NHF = pd.read_csv(os.path.join(HERE, "results", "track_nhf.csv"))

# outcome name -> (description, predicate on (hf flag, n_hf))
CLASSES = {
    "all": lambda hf, n: np.ones(len(hf), bool),
    "h1": lambda hf, n: hf & (n >= 1),
    "h2": lambda hf, n: hf & (n >= 2),
    "h3": lambda hf, n: hf & (n >= 3),
    "b": lambda hf, n: hf & (n == 1),      # brief: a single 6 h HF fix
    "m": lambda hf, n: hf & (n == 2),      # middle class (S3): exactly 2 HF fixes
}
PAIRS = [("h2", "h1"), ("h3", "h1"), ("h3", "b")]     # (a, b) -> log RR(a) - log RR(b); h3 - b is the sustained-minus-brief contrast


def tracks_with_nhf(first, last):
    T = C.load_tracks(first=first, last=last)
    T = T.merge(NHF, on="track", how="left")
    T["n_hf"] = T.n_hf.fillna(0).values
    return T


# ------------------------------------------------------------------ models
class Weekly:
    """Weekly design of PR 64 / PR 41: Y[name] counts per week (660 rows), LP[name] log1p(previous-week count)."""
    kind = "weekly"

    def __init__(self, x, Y, LP):
        self.M = C.month_dummies()
        self.D = CH.Design(self.M, x, Y, LP)
        self.x = x
        self.ns = len(C.SEASONS)

    def rows(self, seas_idx):
        return np.concatenate([np.arange(i * C.NWEEK, (i + 1) * C.NWEEK) for i in seas_idx])

    def slopes(self, names, seas_idx, perm=None):
        rows = self.rows(seas_idx)
        x = self.x
        if perm is not None:
            x = self.x.copy()
            x[rows] = self.x[rows].reshape(len(seas_idx), C.NWEEK)[perm].reshape(-1)
        return self.D.slopes(rows, x, names)

    def clustered_se(self, names_signs, seas_idx):
        rows = self.rows(seas_idx)
        cl = np.repeat(np.arange(len(seas_idx)), C.NWEEK)
        outs = [(self.D.X(n, rows), self.D.Y[n][rows]) for n, _ in names_signs]
        return C.joint_se(outs, cl, [s for _, s in names_signs])


class Daily:
    """Daily design of PR 14 for one basin / index."""
    kind = "daily"

    def __init__(self, T, I, basin, idx, outcomes):
        D = S.window_days("octapr", 2004, 2025)
        Tw = S.assign(T, "octapr", 2004, 2025)
        Tw = Tw[Tw.basin == basin]
        col = idx + "_lag"
        raw = I.reindex(D.date)[[col]].reset_index(drop=True)
        ok = raw.notna().all(axis=1).values
        D, raw = D[ok].reset_index(drop=True), raw[ok].reset_index(drop=True)
        z = (raw - raw.mean()) / raw.std(ddof=0)
        for o in outcomes:
            Tw[o] = CLASSES[o](Tw.hf.values == 1, Tw.n_hf.values).astype(int)
        self.C = S.daily_counts(Tw, D, basin, outcomes)
        self.X = S.design(self.C, z, [col])
        self.season = self.C.season.values
        self.seasons = np.arange(2004, 2026)
        self.Xv = np.asarray(self.X, float)
        self.ns = len(self.seasons)
        self.n = {o: int(self.C[o].sum()) for o in outcomes}

    def rows(self, seas_idx):
        return np.concatenate([np.where(self.season == self.seasons[i])[0] for i in seas_idx])

    def slopes(self, names, seas_idx, perm=None):
        rows = self.rows(seas_idx)
        X = self.Xv[rows].copy()
        if perm is not None:
            blocks = []
            for j, i in enumerate(seas_idx):
                src = self.Xv[self.season == self.seasons[seas_idx[perm[j]]], 0]
                tgt = np.where(self.season == self.seasons[i])[0]
                blocks.append(np.resize(src, len(tgt)))
            X[:, 0] = np.concatenate(blocks)
        return {n: S.fast_pois(self.C[n].values[rows].astype(float), X)[0] for n in names}

    def clustered_se(self, names_signs, seas_idx):
        import statsmodels.api as sm
        rows = self.rows(seas_idx)
        # stacked Poisson with every column interacted with the outcome indicator, as PR 14; contrast of the index terms
        names = [n for n, _ in names_signs]
        k = self.Xv.shape[1]
        blocks, ys, gs = [], [], []
        for j, n in enumerate(names):
            Z = np.zeros((len(rows), k * len(names)))
            Z[:, j * k:(j + 1) * k] = self.Xv[rows]
            blocks.append(Z); ys.append(self.C[n].values[rows]); gs.append(self.season[rows])
        r = sm.GLM(np.concatenate(ys), np.vstack(blocks), family=sm.families.Poisson()).fit(
            cov_type="cluster", cov_kwds={"groups": np.concatenate(gs)})
        c = np.zeros(k * len(names))
        for j, (_, s) in enumerate(names_signs):
            c[j * k] = s
        return float(np.sqrt(c @ r.cov_params() @ c))


# ------------------------------------------------------------------ generic estimation
def analyse(M, names, nperm, nboot, rng, extra_pairs=()):
    """Estimates for every outcome in names and every pair in PAIRS (when both are present)."""
    allidx = list(range(M.ns))
    obs = M.slopes(names, allidx)
    pairs = [p for p in list(PAIRS) + list(extra_pairs) if p[0] in names and p[1] in names]
    tstat = lambda s: {**s, **{f"{a}-{b}": s[a] - s[b] for a, b in pairs}}
    t_obs = tstat(obs)
    keys = list(t_obs)
    null = {k: np.empty(nperm) for k in keys}
    for i in range(nperm):
        t = tstat(M.slopes(names, allidx, perm=rng.permutation(M.ns)))
        for k in keys:
            null[k][i] = t[k]
    bt = {k: np.empty(nboot) for k in keys}
    for i in range(nboot):
        pick = rng.integers(0, M.ns, M.ns)
        t = tstat(M.slopes(names, list(pick)))
        for k in keys:
            bt[k][i] = t[k]
    rows = []
    for k in keys:
        if "-" in k:
            a, b = k.split("-")
            se = M.clustered_se([(a, 1.0), (b, -1.0)], allidx)
        else:
            se = M.clustered_se([(k, 1.0)], allidx)
        lo, hi = np.percentile(bt[k], [2.5, 97.5])
        p = (1 + int((np.abs(null[k]) >= abs(t_obs[k]) - 1e-12).sum())) / (1 + nperm)
        rows.append(dict(stat=k, log_est=t_obs[k], est=np.exp(t_obs[k]), lo=np.exp(lo), hi=np.exp(hi), se_log=se, p_perm=p,
                         mde_rr=np.exp(2.8 * se)))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ builders
def proxy_weekly(basin, Tfull):
    x, _ = C.index(basin)
    t_all = Tfull[Tfull.basin == basin]
    win = C.assign_week(t_all)
    Y, LP = {}, {}
    for name, fn in CLASSES.items():
        m = pd.Series(fn(t_all.hf.values == 1, t_all.n_hf.values), index=t_all.index)
        Y[name] = C.weekly_counts(win, m.reindex(win.index))
        LP[name] = np.log1p(C.prev_week_counts(t_all.assign(_m=m), m, basin))
    return Weekly(x, Y, LP), {k: int(v.sum()) for k, v in Y.items()}


def archive_weekly(basin):
    x, _ = C.index(basin)
    d = json.load(open(os.path.join(REPO, "docs/data/hf-lows.json")))
    f = {k: i for i, k in enumerate(d["lowFields"])}
    # h1 keeps every archive event of class 'low' (4 of 1,829 in 2004-2025 have no HF-category fix), so k = 1 reproduces PR 41 exactly
    days = {k: {} for k in ("h1", "h2", "h3", "b", "m")}
    for r in d["lows"]:
        if r[f["cls"]] != "low" or r[f["basin"]] != basin:
            continue
        n = r[f["hfN"]]
        s = str(r[f["start"]])
        day = dt.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
        for k, ok in (("h1", True), ("h2", n >= 2), ("h3", n >= 3), ("b", n == 1), ("m", n == 2)):
            if ok:
                days[k][day] = days[k].get(day, 0) + 1
    Y, LP = {}, {}
    for k in days:
        y, lp = [], []
        for s in C.SEASONS:
            for w in range(C.NWEEK):
                S0 = C.week_start(s, w)
                y.append(sum(days[k].get(S0 + dt.timedelta(days=i), 0) for i in range(7)))
                lp.append(np.log1p(sum(days[k].get(S0 - dt.timedelta(days=i), 0) for i in range(1, 8))))
        Y[k], LP[k] = np.array(y, float), np.array(lp)
    return Weekly(x, Y, LP), {k: int(v.sum()) for k, v in Y.items()}


def main():
    design, out, nperm, nboot = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    rng = np.random.default_rng(SEED)
    rows = []
    if design == "daily":
        T = S.load_tracks(os.path.join(ERA, "hf_history/results/all_tracks.csv.gz"))
        T = T.merge(NHF, on="track", how="left")
        T["n_hf"] = T.n_hf.fillna(0)
        I = S.indices(CPC, REPO)
        for basin, idx in (("atl", "NAO"), ("pac", "PNA")):
            M = Daily(T, I, basin, idx, list(CLASSES))
            df = analyse(M, [c for c in CLASSES if c != "all"] + ["all"], nperm, nboot, rng)
            df.insert(0, "basin", basin); df.insert(1, "design", f"daily {idx}")
            df["n_tracks"] = M.n["all"]; rows.append(df)
            print(basin, M.n, flush=True)
    else:
        Tfull = tracks_with_nhf(2003, 2025)
        for basin in C.BASINS:
            M, n = proxy_weekly(basin, Tfull) if design == "weekly" else archive_weekly(basin)
            names = list(CLASSES) if design == "weekly" else ["h1", "h2", "h3", "b", "m"]
            df = analyse(M, names, nperm, nboot, rng)
            df.insert(0, "basin", basin); df.insert(1, "design", design)
            rows.append(df)
            print(basin, n, flush=True)
    pd.concat(rows).to_csv(out, index=False, float_format="%.6g")


if __name__ == "__main__":
    main()
