"""Stage B scoring (B.1, B.2) of PREREGISTRATION.md, exactly as written plus Deviations 1-6. Pipeline A is not involved;
the ERA5 index is the z+u-only discovery-fit index (Deviation 6), the ERA5 record is a proxy. Truth: OPC archive weekly counts.
hemispheric/hemlib.py is reused unchanged (load_archive_counts, weekly_table, fit_model/predict -> fit_pois, deviance).
Skill = 1 - D(model with index)/D(baseline), D the Poisson deviance summed over the held-out seasons (PR 41's SS definition,
run.py perm_p/boot_delta); leave-one-season-out over 15 seasons. The index enters as ONE unpenalised regressor, standardised
with the training-fold mean and sd, its coefficient fitted on the training seasons only; index weights are never refitted.
usage: b_score.py [--nsim N] [--selftest]"""
import os, sys, json, time
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hemispheric"))
import hemlib as H

SEASONS = list(range(2004, 2019))
SEED = 20261008
NBOOT = 2000
RR_PLANT = 1.18
RR_GRID = [1.0, 1.05, 1.10, 1.15, 1.18, 1.20, 1.25, 1.30, 1.40, 1.50, 1.70, 2.0, 2.5]


def loso_dev(df, y, z, seasons_arr):
    """Per-season held-out deviance of baseline (month+lprev) and with-index models, leave-one-season-out.
    df: month, lprev columns; y counts; z raw index (None for baseline only)."""
    seas = sorted(set(seasons_arr))
    D0 = np.zeros(len(seas)); D1 = np.zeros(len(seas))
    d = df.assign(y=y, season=seasons_arr)
    for i, s in enumerate(seas):
        tr, te = seasons_arr != s, seasons_arr == s
        m0 = H.fit_model(d[tr], [], None, np.inf)
        D0[i] = H.deviance(y[te], H.predict(m0, d[te], None))
        if z is not None:
            mu_, sd_ = z[tr].mean(), z[tr].std(ddof=1)
            dd = d.assign(gidx=(z - mu_) / sd_)
            m1 = H.fit_model(dd[tr], ["gidx"], None, np.inf)
            D1[i] = H.deviance(y[te], H.predict(m1, dd[te], None))
    return D0, D1


def skill(D0, D1):
    return 1 - D1.sum() / D0.sum()


def boot_skill(D0, D1, draws):
    bs = 1 - D1[draws].sum(1) / D0[draws].sum(1)
    return bs


def summarize(D0, D1, draws):
    obs = skill(D0, D1)
    bs = boot_skill(D0, D1, draws)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    p = (1 + int((bs <= 0).sum())) / (1 + len(bs))
    return dict(skill=float(obs), ci=[float(lo), float(hi)], p_one_sided=float(p), passes=bool(obs > 0 and lo > 0))


def _sim_one(args):
    """One planted-effect simulation: returns bool pass. Same pipeline as the real B.2."""
    df, mu, z, phi, rr, seed, seasons_arr, draws = args
    rng = np.random.default_rng(seed)
    zs = (z - z.mean()) / z.std(ddof=1)
    mean = mu * rr ** zs
    if phi > 1.0001:
        k = mean / (phi - 1.0)
        y = rng.poisson(rng.gamma(k, mean / k))
    else:
        y = rng.poisson(mean)
    y = y.astype(float)
    D0, D1 = loso_dev(df, y, z, seasons_arr)
    return summarize(D0, D1, draws)["passes"]


def power_curve(df, y, z, seasons_arr, draws, nsim, label, pool):
    m0 = H.fit_model(df.assign(y=y), [], None, np.inf)
    mu = H.predict(m0, df, None)
    n, p = len(y), len(m0["beta"])
    phi = max(1.0, float(np.sum((y - mu) ** 2 / mu) / (n - p)))
    res = {}
    for rr in RR_GRID:
        t = time.time()
        args = [(df, mu, z, phi, rr, SEED + 1000 * int(round(rr * 100)) + i, seasons_arr, draws) for i in range(nsim)]
        out = pool.map(_sim_one, args, chunksize=10)
        res[rr] = float(np.mean(out))
        print(f"  power {label} RR={rr}: {res[rr]:.3f} ({time.time()-t:.0f}s)", flush=True)
    pw = [res[r] for r in RR_GRID]
    mde = None
    for j in range(1, len(RR_GRID)):
        if pw[j] >= 0.8 > pw[j - 1] or (pw[j] >= 0.8 and j == 1):
            r0, r1, p0, p1 = RR_GRID[j - 1], RR_GRID[j], pw[j - 1], pw[j]
            mde = r0 + (0.8 - p0) / (p1 - p0) * (r1 - r0) if p1 > p0 else r1
            break
    return dict(dispersion=phi, power_by_rr_per_sd={str(k): v for k, v in res.items()}, power_at_1p18=res[RR_PLANT],
                false_positive_rate_rr1=res[1.0], mde_rr_per_sd_80pct=mde, mde_note=("not reached on grid up to 2.5" if mde is None else "linear interpolation of the power grid"))


def bh(ps):
    ps = np.asarray(ps); m = len(ps); o = np.argsort(ps)
    q = np.empty(m); run = 1.0
    for rank in range(m, 0, -1):
        i = o[rank - 1]; run = min(run, ps[i] * m / rank); q[i] = run
    return q


def main(nsim=500):
    idx = pd.read_csv(os.path.join(HERE, "results", "index_zu.csv"))
    assert len(idx) == 450 and sorted(idx.season.unique()) == SEASONS
    days = H.load_archive_counts(None)
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, len(SEASONS), size=(NBOOT, len(SEASONS)))
    out = dict(seeds=dict(bootstrap=SEED, n_boot=NBOOT, n_sim=nsim), basins={})
    pool = Pool(os.cpu_count())
    for b in H.BASINS:
        df = H.weekly_table(days, SEASONS, b)
        assert (df.season.values == idx.season.values).all() and (df.week.values == idx.week.values).all()
        y = df.y.values.astype(float); sea = df.season.values
        zg, ze = idx[f"{b}_gefs"].values, idx[f"{b}_era5"].values
        r = dict(n_weeks=len(y), total_events=int(y.sum()))
        # B.1
        rho = spearmanr(zg, ze).correlation
        bsr = np.empty(NBOOT)
        for i, dr in enumerate(draws):
            rows = np.concatenate([np.where(sea == SEASONS[j])[0] for j in dr])
            bsr[i] = spearmanr(zg[rows], ze[rows]).correlation
        lo, hi = np.percentile(bsr, [2.5, 97.5])
        p1 = (1 + int((bsr <= 0.5).sum())) / (1 + NBOOT)
        r["B1"] = dict(spearman=float(rho), ci=[float(lo), float(hi)], p_one_sided_le_0p5=float(p1), passes=bool(rho > 0.5),
                       passes_with_ci_above_0p5=bool(lo > 0.5))
        # B.2
        D0, D1 = loso_dev(df, y, zg, sea)
        r["B2"] = summarize(D0, D1, draws)
        r["B2"]["deviance_baseline"] = float(D0.sum()); r["B2"]["deviance_gefs"] = float(D1.sum())
        r["B2"]["per_season_skill"] = {str(s): float(1 - D1[i] / D0[i]) for i, s in enumerate(SEASONS)}
        # reference (not a test): ERA5 z+u index
        _, E1 = loso_dev(df, y, ze, sea)
        ref = dict(all_2004_2018=summarize(D0, E1, draws))
        for nm, sel in (("seasons_2004_2014", np.arange(11)), ("seasons_2015_2018", np.arange(11, 15))):
            d2 = rng.integers(0, len(sel), size=(NBOOT, len(sel)))
            ref[nm] = summarize(D0[sel], E1[sel], d2)
            ref[nm]["note"] = "skill from the same 15-fold leave-one-season-out fits, pooled over these seasons only"
        r["reference_era5_zu_B2"] = ref
        g_sel = np.arange(11, 15)
        r["reference_gefs_2015_2018_subset_skill"] = float(1 - D1[g_sel].sum() / D0[g_sel].sum())
        print(b, r["B1"], r["B2"]["skill"], r["B2"]["ci"], flush=True)
        # power
        r["power"] = power_curve(df, y, zg, sea, draws, nsim, b, pool)
        out["basins"][b] = r
    pool.close()
    ps = [out["basins"][b][t]["p_one_sided_le_0p5" if t == "B1" else "p_one_sided"] for t in ("B1", "B2") for b in H.BASINS]
    names = [f"{t}_{b}" for t in ("B1", "B2") for b in H.BASINS]
    q = bh(ps)
    out["raw_p"] = dict(zip(names, map(float, ps)))
    out["bh_q_within_these_4_only"] = dict(zip(names, map(float, q)))
    out["note_bh"] = "BH over these 4 only; the parent combines the 4 raw p into the 14-test family."
    json.dump(out, open(os.path.join(HERE, "results", "stageB.json"), "w"), indent=1)
    write_txt(out)
    return out


def write_txt(o):
    L = ["Stage B scoring (research/era5/gfs_test/b_score.py). Truth: OPC archive weekly counts, seasons 2004..2018 (450 weeks per basin).",
         "ERA5 index = z+u-only discovery-fit PR 41 index on ERA5 (a proxy; Deviation 6); GEFS index = same frozen weights on GEFS v12 control, model-climatology anomalies.",
         f"Bootstrap: {o['seeds']['n_boot']} whole-season resamples, seed {o['seeds']['bootstrap']}. Power sims: {o['seeds']['n_sim']} per effect level.",
         "Skill = 1 - (held-out Poisson deviance with index)/(held-out deviance of month dummies + log1p(prev week)), pooled over held-out seasons, leave-one-season-out.", ""]
    for b, r in o["basins"].items():
        B1, B2, P = r["B1"], r["B2"], r["power"]
        L += [f"== {b.upper()} (total archive events {r['total_events']}) ==",
              f"B.1 Spearman(GEFS index, ERA5 z+u index) = {B1['spearman']:.3f}  95% CI [{B1['ci'][0]:.3f}, {B1['ci'][1]:.3f}]  one-sided p(rho<=0.5) = {B1['p_one_sided_le_0p5']:.4f}  criterion rho>0.5: {'PASS' if B1['passes'] else 'FAIL'} (CI lower bound above 0.5: {B1['passes_with_ci_above_0p5']})",
              f"B.2 GEFS held-out deviance skill = {100*B2['skill']:.2f}%  95% CI [{100*B2['ci'][0]:.2f}, {100*B2['ci'][1]:.2f}]%  one-sided p(skill<=0) = {B2['p_one_sided']:.4f}  criterion (positive, interval above 0): {'PASS' if B2['passes'] else 'FAIL'}",
              f"    deviances: baseline {B2['deviance_baseline']:.2f}, with GEFS index {B2['deviance_gefs']:.2f}",
              f"    power of B.2 at the PR 41 effect (RR {RR_PLANT}/SD, post hoc): {100*P['power_at_1p18']:.1f}%; false-positive rate at RR 1.0: {100*P['false_positive_rate_rr1']:.1f}%; "
              f"minimum detectable RR/SD at 80% power: {P['mde_rr_per_sd_80pct'] if P['mde_rr_per_sd_80pct'] is None else round(P['mde_rr_per_sd_80pct'],3)} ({P['mde_note']}); dispersion phi {P['dispersion']:.3f}",
              "    power by planted RR per SD: " + ", ".join(f"{k}:{100*v:.0f}%" for k, v in P["power_by_rr_per_sd"].items()),
              "    REFERENCE, NOT A TEST: same skill with the ERA5 z+u index on the same weeks:"]
        for nm, s in r["reference_era5_zu_B2"].items():
            L.append(f"      {nm}: {100*s['skill']:.2f}%  95% CI [{100*s['ci'][0]:.2f}, {100*s['ci'][1]:.2f}]%  p(skill<=0)={s['p_one_sided']:.4f}")
        L += [f"    (GEFS-index skill on seasons 2015-2018 alone, from the same folds: {100*r['reference_gefs_2015_2018_subset_skill']:.2f}%)", ""]
    L += ["Raw one-sided p values (to be combined by the parent into the 14-test BH family): " + ", ".join(f"{k}={v:.4f}" for k, v in o["raw_p"].items()),
          "BH q over these 4 only (NOT the 14-test family): " + ", ".join(f"{k}={v:.4f}" for k, v in o["bh_q_within_these_4_only"].items()),
          "A B.2 null with power under 80% is 'cannot tell', not 'no effect' (preregistration)."]
    open(os.path.join(HERE, "results", "stageB.txt"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    n = int(sys.argv[sys.argv.index("--nsim") + 1]) if "--nsim" in sys.argv else 500
    main(n)
