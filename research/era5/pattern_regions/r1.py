"""R1: replication on outcomes no region model has met. See PREREGISTRATION.md.

Each primary model (FULL, keep-only and drop-one sectors) is refitted on all 22 archive seasons (EOFs, lambda by inner
leave-one-season-out), its penalised-PC linear predictor is standardised over 1979-80..2000-01 and tested as in
hemispheric S5: pipeline A *depth* counts (proxy outcome; cuts 966.2 hPa Atlantic, 965.0 hPa Pacific), month + previous
week + gamma x index, one-sided season-block permutation (5,000). Within-era variation only (decision 1).

usage: r1.py WORKDIR RESULTDIR
"""
import os, sys, json, datetime as dt
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hemispheric"))
sys.path.insert(0, HERE)
import hemlib as H
import run as R
import regions as RG

ALL = list(range(2004, 2026))
ERA = list(range(1979, 2001))
NPERM = 5000
NBOOT = 10000


def main(work, res):
    os.makedirs(res, exist_ok=True)
    ctx = R.Ctx(work)
    days = H.load_archive_counts(None)
    depth = H.load_proxy_days("depth")
    rng = np.random.default_rng(H.SEED + 21)
    draws = rng.integers(0, len(ERA), size=(NBOOT, len(ERA)))
    out = {}
    for b in H.BASINS:
        dfa = H.weekly_table(days, ALL, b)
        dfe = H.weekly_table(depth, ERA, b)
        y = dfe.y.values.astype(float)
        ns = len(ERA)
        rows = {}
        for sp in RG.specs(b):
            if sp["fam"] != "primary" or sp["base"] != "b1":
                continue
            eof = RG.REOFs(ctx.maps_for("lag1", ALL), ctx.mask, sp["cells"], sp["vars"])
            pcs = eof.scores(ctx.maps_for("lag1", ALL))
            lam, _ = H.loso_lambda(dfa, [], pcs)
            mdl = H.fit_model(dfa, [], pcs, lam)
            if not np.isfinite(lam):
                rows[sp["name"]] = dict(lam=None, note="no PCs retained")
                continue
            pe = eof.scores(ctx.maps_for("lag1", ERA))
            _, idx = H.predict(mdl, dfe, pe, parts=True)
            idx = (idx - idx.mean()) / idx.std()
            d1 = dfe.assign(ix=idx)
            m0 = H.fit_model(d1, [], None, np.inf)
            m1 = H.fit_model(d1, ["ix"], None, np.inf)
            mu0, mu1 = H.predict(m0, d1, None), H.predict(m1, d1, None)
            gam = float(m1["beta"][8])
            gain = H.deviance(y, mu0) - H.deviance(y, mu1)
            sd0 = np.array([H.deviance(y[i * 30:(i + 1) * 30], mu0[i * 30:(i + 1) * 30]) for i in range(ns)])
            sd1 = np.array([H.deviance(y[i * 30:(i + 1) * 30], mu1[i * 30:(i + 1) * 30]) for i in range(ns)])
            null = []
            for _ in range(NPERM):
                pm = rng.permutation(ns)
                dp = d1.assign(ix=idx.reshape(ns, 30)[pm].reshape(-1))
                null.append(float(H.fit_model(dp, ["ix"], None, np.inf)["beta"][8]))
            p = (1 + int((np.array(null) >= gam).sum())) / (1 + NPERM)
            rows[sp["name"]] = dict(lam=float(lam), gamma=gam, rr_per_sd=float(np.exp(gam)), dev_gain=float(gain), p=float(p),
                                    season_gain=(sd0 - sd1).tolist(), events=int(y.sum()))
            print(b, sp["name"], f"lam={lam} RR/SD={np.exp(gam):.3f} gain={gain:.2f} p={p:.4f}", flush=True)
        if "FULL" in rows and rows["FULL"].get("dev_gain", 0) > 0:
            gF = np.array(rows["FULL"]["season_gain"])
            gFb = gF[draws].sum(1)
            for nm, r in rows.items():
                if "season_gain" not in r:
                    continue
                g = np.array(r["season_gain"])[draws].sum(1)
                ok = gFb > 0
                r["share_of_full"] = float(r["dev_gain"] / rows["FULL"]["dev_gain"])
                r["share_ci90"] = [float(np.percentile((g / gFb)[ok], 5)), float(np.percentile((g / gFb)[ok], 95))]
        out[b] = rows
    json.dump(out, open(os.path.join(res, "r1.json"), "w"), indent=1)
    with open(os.path.join(res, "looks.log"), "a") as f:
        f.write(f"{dt.datetime.utcnow().isoformat()}Z R1 scored on pipeline A depth counts, seasons {ERA[0]}..{ERA[-1]} (proxy outcome, within-era)\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
