"""Skill, shares, Shapley decomposition, bootstrap, permutation tests and FDR for the region ablation.
Reads the npz files written by ablate.py (no fields needed). See PREREGISTRATION.md for every rule.

usage: stats.py OUTDIR RESULTDIR
"""
import os, sys, json, itertools
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hemispheric"))
sys.path.insert(0, HERE)
import hemlib as H
import run as R
import regions as RG

NS = 22
NBOOT = 10000
NPERM = 10000
BASINS = H.BASINS


def season_dev(y, eta):
    mu = np.exp(np.clip(eta, -30, 30))
    return np.array([H.deviance(y[i * 30:(i + 1) * 30], mu[i * 30:(i + 1) * 30]) for i in range(NS)])


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    run = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        run = min(run, p[i] * n / rank); q[i] = run
    return q


def main(out, res):
    os.makedirs(res, exist_ok=True)
    rng = np.random.default_rng(H.SEED + 11)
    draws = rng.integers(0, NS, size=(NBOOT, NS))
    summary, tests, curve = {}, [], []
    for b in BASINS:
        ref = np.load(os.path.join(out, f"refs_{b}.npz"))
        y = ref["y"]
        D1 = season_dev(y, ref["eta_b1"]); D2 = season_dev(y, ref["eta_b2"])
        models = {}
        for sp in RG.specs(b):
            fn = os.path.join(out, f"{b}_{sp['name']}.npz")
            if not os.path.exists(fn):
                continue
            z = np.load(fn)
            assert np.array_equal(z["y"], y)
            models[sp["name"]] = dict(spec=sp, eb=z["eta_base"], ep=z["eta_pc"], lam=z["lam"],
                                      D=season_dev(y, z["eta_base"] + z["eta_pc"]))
        base = lambda m: D2 if m["spec"]["base"] == "b2" else D1
        gain = lambda m: base(m).sum() - m["D"].sum()                      # deviance units
        ss = lambda m: gain(m) / base(m).sum()
        bgain = lambda m: base(m)[draws].sum(1) - m["D"][draws].sum(1)

        info = {}
        # permutation p for SS > 0 (all models)
        for nm, m in models.items():
            D_ref = base(m).sum()
            obs, p, _ = R.perm_p(y, m["eb"], m["ep"], D_ref, rng, NS, n=NPERM)
            bs = bgain(m) / base(m)[draws].sum(1)
            info[nm] = dict(fam=m["spec"]["fam"], ss=ss(m), p_perm=p, ss_ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                            seasons_beat_ref=int((m["D"] < base(m)).sum()), lam_none=int(np.isinf(m["lam"]).sum()))
            tests.append(dict(basin=b, family=("F1" if m["spec"]["fam"] == "primary" else "F3"), test=f"SS>0 {nm}", stat=ss(m), p=p))
        if "FULL" not in models:
            continue
        F = models["FULL"]; gF = gain(F); bgF = bgain(F)
        # shares and Shapley (primary)
        prim = ["only_UP", "only_LOC", "only_DN", "no_UP", "no_LOC", "no_DN"]
        if all(k in models for k in prim):
            def G_of(idx_gain):
                return {frozenset(): 0.0, frozenset({"UP"}): idx_gain["only_UP"], frozenset({"LOC"}): idx_gain["only_LOC"],
                        frozenset({"DN"}): idx_gain["only_DN"], frozenset({"LOC", "DN"}): idx_gain["no_UP"],
                        frozenset({"UP", "DN"}): idx_gain["no_LOC"], frozenset({"UP", "LOC"}): idx_gain["no_DN"],
                        frozenset({"UP", "LOC", "DN"}): None}
            def shap(g):
                Gd = G_of(g); Gd[frozenset({"UP", "LOC", "DN"})] = g["FULL"]
                phi = {}
                import math
                for i in ("UP", "LOC", "DN"):
                    others = [p for p in ("UP", "LOC", "DN") if p != i]
                    v = 0.0
                    for r in range(3):
                        for S in itertools.combinations(others, r):
                            w = math.factorial(r) * math.factorial(2 - r) / 6.0
                            v += w * (Gd[frozenset(S) | {i}] - Gd[frozenset(S)])
                    phi[i] = v
                return phi
            g_obs = {k: gain(models[k]) for k in prim + ["FULL"]}
            g_bs = {k: bgain(models[k]) for k in prim + ["FULL"]}
            phi = shap(g_obs)
            share = {k: v / gF for k, v in phi.items()}
            ok = bgF > 0
            bsphi = {k: [] for k in phi}
            for d in np.where(ok)[0]:
                ph = shap({k: g_bs[k][d] for k in g_bs})
                for k in ph:
                    bsphi[k].append(ph[k] / g_bs["FULL"][d])
            share_ci = {k: ([float(np.percentile(v, 5)), float(np.percentile(v, 95))] if len(v) else None) for k, v in bsphi.items()}
            f = {k: gain(models["only_" + k]) / gF for k in ("UP", "LOC", "DN")}
            f_bs = {k: (g_bs["only_" + k] / bgF)[ok] for k in ("UP", "LOC", "DN")}
            f_ci = {k: ([float(np.percentile(v, 5)), float(np.percentile(v, 95))] if len(v) else [float('nan')] * 2) for k, v in f_bs.items()}
            drop = {k: gF - gain(models["no_" + k]) for k in ("UP", "LOC", "DN")}
            drop_bs = {k: bgF - g_bs["no_" + k] for k in drop}
            drop_share = {k: drop[k] / gF for k in drop}
            # pairwise keep-only differences (deviance units / D_B1 total -> skill points)
            pair = {}
            for a, c in (("UP", "LOC"), ("UP", "DN"), ("LOC", "DN")):
                d_obs = (gain(models["only_" + a]) - gain(models["only_" + c])) / D1.sum()
                d_bs = (g_bs["only_" + a] - g_bs["only_" + c]) / D1[draws].sum(1)
                p2 = min(1.0, 2 * min((1 + (d_bs <= 0).sum()) / (NBOOT + 1), (1 + (d_bs >= 0).sum()) / (NBOOT + 1)))
                pair[f"{a}-{c}"] = dict(diff=float(d_obs), ci95=[float(np.percentile(d_bs, 2.5)), float(np.percentile(d_bs, 97.5))],
                                         se=float(d_bs.std()), p=float(p2), mdd80=float(2.8 * d_bs.std()))
                tests.append(dict(basin=b, family="F2", test=f"keep-only {a}-{c}", stat=float(d_obs), p=float(p2)))
            drop_p = {}
            for k in drop:
                d_ = drop_bs[k] / D1[draws].sum(1)
                pp = (1 + (d_ <= 0).sum()) / (NBOOT + 1)
                drop_p[k] = dict(loss_skill=float(drop[k] / D1.sum()), ci95=[float(np.percentile(d_, 2.5)), float(np.percentile(d_, 97.5))],
                                 share_of_full=float(drop_share[k]), p=float(pp))
                tests.append(dict(basin=b, family="F2", test=f"drop {k}: SS(FULL)-SS(no_{k})>0", stat=float(drop[k] / D1.sum()), p=float(pp)))
            # split halves (descriptive)
            halves = {}
            for nmh, sl in (("2004-14", slice(0, 11)), ("2015-25", slice(11, 22))):
                gFh = D1[sl].sum() - F["D"][sl].sum()
                halves[nmh] = dict(ss_full=float(gFh / D1[sl].sum()),
                                   f_UP=float((D1[sl].sum() - models["only_UP"]["D"][sl].sum()) / gFh) if gFh > 0 else None,
                                   f_LOC=float((D1[sl].sum() - models["only_LOC"]["D"][sl].sum()) / gFh) if gFh > 0 else None)
            cat = None
            fu, fl = f["UP"], f["LOC"]
            if gF <= 0: cat = "NO SKILL"
            elif fu >= .5 and fl < .5: cat = "REMOTE"
            elif fl >= .5 and fu < .5: cat = "LOCAL"
            elif fu >= .5 and fl >= .5: cat = "SHARED"
            else: cat = "DISTRIBUTED"
            side = lambda lo, hi: "undefined" if not (lo == lo) else ("below" if hi < .5 else ("above" if lo >= .5 else "straddles"))
            resolved = {k: side(*f_ci[k]) for k in ("UP", "LOC")}
            summary[b] = dict(ss_full=ss(F), ss_full_ci95=info["FULL"]["ss_ci95"], p_full=info["FULL"]["p_perm"],
                              bootstrap_valid_frac=float(ok.mean()),
                              f_keep_only=f, f_ci90=f_ci, category=cat, interval_side_of_half=resolved,
                              shapley_share=share, shapley_ci90=share_ci, drop_one=drop_p, pairwise=pair, halves=halves,
                              f_DN_ge_half=bool(f["DN"] >= .5))
        # V1/V2/V3/imprint extras
        for nm, m in models.items():
            if m["spec"]["fam"] in ("V1", "V2", "V3", "imprint"):
                info[nm]["f_of_full"] = float(gain(m) / gF) if m["spec"]["base"] == "b1" else float(gain(m) / gain(models["FULL_b2"])) if "FULL_b2" in models else None
        summary.setdefault(b, {})["models"] = info
        for nm, m in models.items():
            if m["spec"]["fam"] == "V3":
                curve.append(dict(basin=b, centre=float(nm.split("_")[1]), ss=ss(m), p=info[nm]["p_perm"]))
    # FDR: within families, and across everything, per basin and overall
    T = pd.DataFrame(tests)
    T["q_family"] = np.nan; T["q_basin"] = np.nan
    for (b, fam), g in T.groupby(["basin", "family"]):
        T.loc[g.index, "q_family"] = bh(g.p.values)
    for b, g in T.groupby("basin"):
        T.loc[g.index, "q_basin"] = bh(g.p.values)
    T["q_all"] = bh(T.p.values)
    T.to_csv(os.path.join(res, "tests.csv"), index=False)
    pd.DataFrame(curve).to_csv(os.path.join(res, "window_curve.csv"), index=False)
    cnt = {b: dict(total=int(len(g)), pass_q_basin=int((g.q_basin < .05).sum()), pass_q_all=int((g.q_all < .05).sum())) for b, g in T.groupby("basin")}
    summary["test_counts"] = cnt
    summary["test_counts"]["all"] = dict(total=int(len(T)), pass_q_all=int((T.q_all < .05).sum()))
    json.dump(summary, open(os.path.join(res, "summary.json"), "w"), indent=1, default=float)
    print(json.dumps({b: {k: v for k, v in summary[b].items() if k != "models"} for b in BASINS if b in summary}, indent=1, default=float))
    print(cnt)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
