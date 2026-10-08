"""POST HOC checks, run after the pre-registered held-out look. None of these was in PREREGISTRATION.md.
Each re-scores the same held-out seasons (2015-16..2025-26), so they are further looks; the count is logged.

  PH1  storm-track baseline: B2 = B1 + log(1 + all pipeline A cyclones in the basin in the previous 7 days)
  PH2  seasonal-cycle baseline: B1 + first two annual harmonics of week of season
  PH3  one field at a time (Z500, U250, SST) and leave-one-field-out
  PH4  effect size: held-out HF rate by quintile of the pattern index, and RR per SD of the index
usage: posthoc.py WORKDIR [OUT]
"""
import os, sys, json, pickle, datetime as dt
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hemlib as H
import run as R


def all_cyclone_days():
    T = pd.read_csv(os.path.join(H.REPO, "research/era5/hf_history/results/all_tracks.csv.gz"), dtype={"start": str})
    days = {b: {} for b in H.BASINS}
    for b in H.BASINS:
        for st in T[T.basin == b].start:
            d = dt.date(int(st[:4]), int(st[4:6]), int(st[6:8]))
            days[b][d] = days[b].get(d, 0) + 1
    return days


def add_cols(df, seasons, b, allc):
    prev = []
    for s, k in zip(df.season, df.week):
        S = H.week_start(s, k)
        prev.append(sum(allc[b].get(S - dt.timedelta(days=i), 0) for i in range(1, 8)))
    df["lall"] = np.log1p(prev)
    w = (df.week.values + 0.5) / H.NWEEK * 2 * np.pi * (210 / 365.25)     # fraction of the Oct-Apr window as phase
    ph = 2 * np.pi * (df.week.values * 7 + 3) / 365.25
    df["h1c"], df["h1s"], df["h2c"], df["h2s"] = np.cos(ph), np.sin(ph), np.cos(2 * ph), np.sin(2 * ph)
    return df


def score(df_tr, df_te, pcs_tr, pcs_te, extra, rng, nseas, label):
    base = H.fit_model(df_tr, extra, None, np.inf)
    lam, tot = H.loso_lambda(df_tr, extra, pcs_tr)
    pat = H.fit_model(df_tr, extra, pcs_tr, lam)
    y = df_te.y.values.astype(float)
    mu0 = H.predict(base, df_te, None); muP = H.predict(pat, df_te, pcs_te if np.isfinite(lam) else None)
    D0, DP = H.deviance(y, mu0), H.deviance(y, muP)
    bP, eP, pP = R.parts(pat, df_te, pcs_te)
    obs, p, _ = R.perm_p(y, bP + eP, pP, D0, rng, nseas)
    return dict(label=label, lam=float(lam), ss=1 - DP / D0, p=p, cv_gain=1 - min(tot.values()) / tot[np.inf])


def main():
    work = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(H.HERE, "results")
    ctx = R.Ctx(work)
    fz = pickle.load(open(os.path.join(work, "frozen_primary.pkl"), "rb"))
    days = H.load_archive_counts(None)
    allc = all_cyclone_days()
    train, test = fz["train"], fz["test"]
    eof = fz["eof"]
    pcs_tr = eof.scores(ctx.maps_for("lag1", train)); pcs_te = eof.scores(ctx.maps_for("lag1", test))
    ix_tr, ix_te = ctx.index_table(train, "lag1"), ctx.index_table(test, "lag1")
    rng = np.random.default_rng(H.SEED + 11)
    res = {}
    for b in H.BASINS:
        dtr, st = R.add_indices(H.weekly_table(days, train, b), ix_tr)
        dte = R.apply_stats(H.weekly_table(days, test, b), ix_te, st)
        dtr, dte = add_cols(dtr, train, b, allc), add_cols(dte, test, b, allc)
        r = {}
        r["PH0_primary_check"] = score(dtr, dte, pcs_tr, pcs_te, [], rng, len(test), "primary (should reproduce P1)")
        r["PH1_stormtrack"] = score(dtr, dte, pcs_tr, pcs_te, ["lall"], rng, len(test), "baseline B1 + prior-week all-cyclone count")
        r["PH2_seasonal"] = score(dtr, dte, pcs_tr, pcs_te, ["h1c", "h1s", "h2c", "h2s"], rng, len(test), "baseline B1 + annual harmonics")
        for nm, sl in (("z500", slice(0, 10)), ("u250", slice(10, 20)), ("sst", slice(20, 30))):
            r[f"PH3_{nm}_only"] = score(dtr, dte, pcs_tr[:, sl], pcs_te[:, sl], [], rng, len(test), f"{nm} PCs only")
        for nm, drop in (("no_z500", slice(0, 10)), ("no_u250", slice(10, 20)), ("no_sst", slice(20, 30))):
            keep = np.ones(30, bool); keep[drop] = False
            r[f"PH3_{nm}"] = score(dtr, dte, pcs_tr[:, keep], pcs_te[:, keep], [], rng, len(test), f"leave out {nm}")
        # PH4 effect size from the frozen primary model
        m = fz["basins"][b]["models"]
        df = R.apply_stats(H.weekly_table(days, test, b), ix_te, fz["basins"][b]["stats"])
        y = df.y.values.astype(float)
        base, extra, pc = R.parts(m["P"], df, pcs_te)
        # index scaled by the discovery-period SD of the pattern index
        idx_tr = R.parts(m["P"], H.weekly_table(days, train, b).assign(**{c: 0.0 for c in R.NAMED}), pcs_tr)[2]
        sd = float(idx_tr.std()); z = (pc - idx_tr.mean()) / sd
        mu1 = H.predict(m["B1"], df, None)
        # RR per SD: Poisson GLM of y on z with offset log(mu1); season-block bootstrap
        def gam(rows):
            g = 0.0
            for _ in range(30):
                mu = mu1[rows] * np.exp(g * z[rows])
                grad = (z[rows] * (y[rows] - mu)).sum(); hes = (z[rows] ** 2 * mu).sum()
                g += grad / hes
            return g
        allrows = np.arange(len(y)); g0 = gam(allrows)
        bs = []
        for _ in range(2000):
            pick = rng.integers(0, len(test), len(test))
            rows = np.concatenate([np.arange(i * H.NWEEK, (i + 1) * H.NWEEK) for i in pick])
            bs.append(gam(rows))
        q = np.quantile(z, [0.2, 0.8])
        lo, hi = z <= q[0], z >= q[1]
        r["PH4_effect"] = dict(rr_per_sd=float(np.exp(g0)), ci=(float(np.exp(np.percentile(bs, 2.5))), float(np.exp(np.percentile(bs, 97.5)))),
                               top_quintile_rate=float(y[hi].mean()), bottom_quintile_rate=float(y[lo].mean()),
                               ratio=float(y[hi].mean() / y[lo].mean()), n_top=int(hi.sum()), n_bottom=int(lo.sum()),
                               index_sd_training=sd,
                               ratio_vs_B1_expected=float((y[hi].sum() / mu1[hi].sum()) / (y[lo].sum() / mu1[lo].sum())))
        res[b] = r
        print(b, json.dumps({k: (v if k == 'PH4_effect' else {a: (round(c, 4) if isinstance(c, float) else c) for a, c in v.items()}) for k, v in r.items()}, default=float, indent=0)[:3000], flush=True)
    json.dump(res, open(os.path.join(out, "posthoc.json"), "w"), indent=1, default=float)
    with open(os.path.join(out, "heldout_looks.log"), "a") as f:
        f.write(f"{dt.datetime.now(dt.timezone.utc).isoformat()} POSTHOC PH0-PH4 scored on 2015..2025 (primary window)\n")


if __name__ == "__main__":
    main()
