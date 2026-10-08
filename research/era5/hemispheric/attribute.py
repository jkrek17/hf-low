"""Attribution of the pattern to named indices, out-of-sample index, maps. See PREREGISTRATION.md.

usage: attribute.py WORKDIR [OUT]
Stages (all in one run, about 10 minutes):
  1. Leave-one-season-out pattern index over all 22 seasons (EOFs, lambda re-chosen inside each fold's 21 seasons).
  2. Attribution regressions A1 (NAO, PNA, ONI, MJO1, MJO2), A2 (+ three products), A3 (+ AO) with month effects;
     adjusted R-squared, 95% season-block bootstrap interval, standardised coefficients.
  3. Descriptive LOSO skill over all 22 seasons (pattern vs B1 vs named), not the pre-registered split.
  4. Maps: composite (top quintile minus bottom quintile of the out-of-sample index) of Z500, U250 and SST, with
     BH-FDR field significance from a season-block permutation; discovery-half vs held-out-half map correlation.
"""
import os, sys, json, time, datetime as dt
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hemlib as H
import run as R

ALL = list(range(2004, 2026))
NBOOT_A = 2000
NPERM_MAP = 1000


def fit_fold(ctx, days, b, train, win="lag1"):
    eof = H.EOFs(ctx.maps_for(win, train), ctx.mask)
    pcs = eof.scores(ctx.maps_for(win, train))
    ix = ctx.index_table(train, win)
    df, st = R.add_indices(H.weekly_table(days, train, b), ix)
    m = {"B1": H.fit_model(df, [], None, np.inf), "N": H.fit_model(df, R.NAMED, None, np.inf)}
    lam, tot = H.loso_lambda(df, [], pcs)
    m["P"] = H.fit_model(df, [], pcs, lam)
    return eof, m, st, lam


def oos_index(ctx, days, b, win="lag1"):
    """For each season: the pattern model's penalised-PC linear predictor from a fit on the other 21 seasons,
    plus the held-out deviance of P, B1 and N."""
    rows, dev = [], []
    for s in ALL:
        train = [x for x in ALL if x != s]
        eof, m, st, lam = fit_fold(ctx, days, b, train, win)
        pcs = eof.scores(ctx.maps_for(win, [s]))
        df = R.apply_stats(H.weekly_table(days, [s], b), ctx.index_table([s], win), st)
        base, extra, pc = R.parts(m["P"], df, pcs)
        y = df.y.values.astype(float)
        dev.append(dict(season=s, lam=lam, D_B1=H.deviance(y, H.predict(m["B1"], df, None)),
                        D_N=H.deviance(y, H.predict(m["N"], df, None)), D_P=H.deviance(y, H.predict(m["P"], df, pcs)),
                        events=int(y.sum())))
        for k in range(H.NWEEK):
            rows.append(dict(season=s, week=k, idx=pc[k], y=y[k]))
        print(f"  oos {b} {s} lam={lam}", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(dev)


def attribution(ctx, oos, win="lag1"):
    ix = ctx.index_table(ALL, win)
    df = pd.DataFrame(dict(season=oos.season.values, week=oos.week.values))
    df["month"] = [(H.week_start(s_, k_) + dt.timedelta(days=3)).month for s_, k_ in zip(df.season, df.week)]
    df, _ = R.add_indices(df, ix)
    out = {}
    M = pd.get_dummies(df.month).reindex(columns=[10, 11, 12, 1, 2, 3, 4], fill_value=0).values.astype(float)
    y = oos.idx.values
    specs = {"A1": R.NAMED, "A2": R.NAMED + R.PROD, "A3": R.NAMED + R.PROD + ["ao"]}
    def r2(cols, rows):
        X = np.hstack([M[rows]] + ([df[cols].values[rows]] if cols else []))
        yy = y[rows]
        beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
        res = yy - X @ beta
        n, p = X.shape
        tot = ((yy - yy.mean()) ** 2).sum()
        return 1 - (res @ res) / tot, 1 - (res @ res) / (n - p) / (tot / (n - 1)), beta
    month_r2 = r2([], np.arange(len(y)))[1]
    seas = np.array(ALL)
    rng = np.random.default_rng(H.SEED + 7)
    for k, cols in specs.items():
        full = np.arange(len(y))
        R2, adj, beta = r2(cols, full)
        bs = []
        for _ in range(NBOOT_A):
            pick = rng.integers(0, len(ALL), len(ALL))
            rows = np.concatenate([np.arange(i * H.NWEEK, (i + 1) * H.NWEEK) for i in pick])
            bs.append(r2(cols, rows)[1])
        sd = y.std()
        coef = {c: float(beta[7 + j] / sd) for j, c in enumerate(cols)}
        out[k] = dict(r2=float(R2), adj_r2=float(adj), adj_r2_ci=(float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))),
                      std_coef_per_sd_index=coef, adj_r2_month_only=float(month_r2))
    # simple correlations
    out["corr"] = {c: float(np.corrcoef(y, df[c].values)[0, 1]) for c in R.NAMED + ["ao"]}
    # within-month (month-demeaned) correlations
    ym = y - M @ np.linalg.lstsq(M, y, rcond=None)[0]
    out["corr_month_removed"] = {}
    for c in R.NAMED + ["ao"]:
        xc = df[c].values
        xm = xc - M @ np.linalg.lstsq(M, xc, rcond=None)[0]
        out["corr_month_removed"][c] = float(np.corrcoef(ym, xm)[0, 1])
    return out


def composite_maps(ctx, oos_by_basin, win="lag1"):
    maps = ctx.maps_for(win, ALL)
    res = {}
    rng = np.random.default_rng(H.SEED + 8)
    ns = len(ALL)
    for b, oos in oos_by_basin.items():
        idx = oos.idx.values
        out = {}
        for v in H.VARS:
            X = maps[v]
            def comp(ix):
                hi, lo = ix >= np.quantile(ix, 0.8), ix <= np.quantile(ix, 0.2)
                return np.nanmean(X[hi], 0) - np.nanmean(X[lo], 0)
            obs = comp(idx)
            null = np.empty((NPERM_MAP,) + obs.shape)
            for i in range(NPERM_MAP):
                pm = rng.permutation(ns)
                null[i] = comp(idx.reshape(ns, H.NWEEK)[pm].reshape(-1))
            p = (1 + (np.abs(null) >= np.abs(obs)[None]).sum(0)) / (1 + NPERM_MAP)
            valid = ~np.isnan(obs)
            q = np.full(obs.shape, np.nan)
            pv = p[valid]; order = np.argsort(pv); n = len(pv)
            adj = pv[order] * n / (np.arange(n) + 1)
            adj = np.minimum.accumulate(adj[::-1])[::-1]
            qq = np.empty(n); qq[order] = np.minimum(adj, 1)
            q[valid] = qq
            # halves
            half = {}
            for name, ss in (("disc", range(0, 11)), ("held", range(11, 22))):
                rows = np.concatenate([np.arange(i * H.NWEEK, (i + 1) * H.NWEEK) for i in ss])
                ix = idx[rows]; XX = X[rows]
                hi, lo = ix >= np.quantile(ix, 0.8), ix <= np.quantile(ix, 0.2)
                half[name] = np.nanmean(XX[hi], 0) - np.nanmean(XX[lo], 0)
            w = np.cos(np.radians(H.LAT))[:, None] * np.ones((1, 64))
            m = ~np.isnan(half["disc"]) & ~np.isnan(half["held"])
            a, c = half["disc"][m], half["held"][m]; ww = w[m]
            a = a - (a * ww).sum() / ww.sum(); c = c - (c * ww).sum() / ww.sum()
            pc = float((a * c * ww).sum() / np.sqrt((a * a * ww).sum() * (c * c * ww).sum()))
            out[v] = dict(diff=obs, q=q, disc=half["disc"], held=half["held"], pattern_corr_halves=pc,
                          frac_sig=float(np.nanmean(q < 0.05)))
        res[b] = out
    return res


def main():
    work = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(H.HERE, "results")
    ctx = R.Ctx(work)
    days = H.load_archive_counts(None)
    oos, devs, attr = {}, {}, {}
    for b in H.BASINS:
        t = time.time()
        oos[b], devs[b] = oos_index(ctx, days, b)
        oos[b].to_csv(os.path.join(out, f"oos_index_{b}.csv"), index=False)
        attr[b] = attribution(ctx, oos[b])
        d = devs[b]
        attr[b]["loso22"] = dict(SS_P_B1=float(1 - d.D_P.sum() / d.D_B1.sum()), SS_N_B1=float(1 - d.D_N.sum() / d.D_B1.sum()),
                                 seasons_P_beats_B1=int((d.D_P < d.D_B1).sum()), seasons_N_beats_B1=int((d.D_N < d.D_B1).sum()),
                                 events=int(d.events.sum()))
        print(b, json.dumps(attr[b], default=float)[:900], f"{time.time()-t:.0f}s", flush=True)
    json.dump(attr, open(os.path.join(out, "attribution.json"), "w"), indent=1, default=float)
    cm = composite_maps(ctx, oos)
    np.savez_compressed(os.path.join(out, "composite_maps.npz"),
                        **{f"{b}_{v}_{k}": cm[b][v][k] for b in cm for v in cm[b] for k in ("diff", "q", "disc", "held")})
    summ = {b: {v: dict(pattern_corr_halves=cm[b][v]["pattern_corr_halves"], frac_sig=cm[b][v]["frac_sig"]) for v in cm[b]} for b in cm}
    json.dump(summ, open(os.path.join(out, "composite_summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
