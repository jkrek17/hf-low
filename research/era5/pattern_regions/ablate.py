"""Nested leave-one-season-out scoring of region-restricted pattern models. See PREREGISTRATION.md.

usage: ablate.py WORKDIR OUTDIR [NPROC] [name-filter]
  WORKDIR  fields from research/era5/hemispheric/fields.py (season_*.npz)
  OUTDIR   one npz per (basin, model) plus refs_<basin>.npz; resumable (existing files are skipped)

For each held-out season s of the 22 archive seasons 2004-05..2025-26: EOFs (region-restricted) and lambda (inner
leave-one-season-out) from the other 21 seasons, refit on those 21, score season s. Stored per week: archive count,
the baseline linear predictor, the penalised-PC part, and lambda per season. The reference B1 (and B2) are refitted per
fold the same way. Scheme identical to hemispheric/attribute.py; FULL must reproduce its +5.74% (Atlantic) and +3.40%
(Pacific) before anything else is reported.
"""
import os, sys, time, datetime as dt
import numpy as np, pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hemispheric"))
sys.path.insert(0, HERE)
import hemlib as H
import run as R
import posthoc as PH
import regions as RG

ALL = list(range(2004, 2026))
G = {}


def init(work):
    ctx = R.Ctx(work)
    ctx.maps("lag1")
    G["ctx"] = ctx
    G["days"] = H.load_archive_counts(None)
    G["allc"] = PH.all_cyclone_days()


def fast_loso(df, extra, pcs, lambdas=H.LAMBDAS):
    """Same result as H.loso_lambda, with the design matrix built once."""
    y = df["y"].values.astype(float)
    seas = df["season"].values
    Xp, nu = H.design(df, extra, pcs)
    X0, _ = H.design(df, extra, None)
    tot = {}
    for lam in lambdas:
        X = Xp if np.isfinite(lam) else X0
        pen = np.zeros(X.shape[1]); pen[nu:] = lam
        d = 0.0
        for s in np.unique(seas):
            tr, te = seas != s, seas == s
            b = H.fit_pois(X[tr], y[tr], pen)
            d += H.deviance(y[te], np.exp(np.clip(X[te] @ b, -30, 30)))
        tot[lam] = d if np.isfinite(d) else float("inf")
    return min(tot, key=tot.get), tot


def table(basin, base):
    df = H.weekly_table(G["days"], ALL, basin)
    return PH.add_cols(df, ALL, basin, G["allc"])


def refs(basin):
    """LOSO B1 and B2 linear predictors (no pattern)."""
    df = table(basin, "b2")
    out = {k: np.zeros(len(df)) for k in ("eta_b1", "eta_b2")}
    for s in ALL:
        tr, te = (df.season != s).values, (df.season == s).values
        for key, extra in (("eta_b1", []), ("eta_b2", ["lall"])):
            m = H.fit_model(df[tr], extra, None, np.inf)
            out[key][te] = np.log(H.predict(m, df[te], None))
    return df.y.values.astype(float), out


def run_model(spec):
    basin, name = spec["basin"], spec["name"]
    fn = os.path.join(G["out"], f"{basin}_{name}.npz")
    if os.path.exists(fn):
        return fn
    t0 = time.time()
    ctx = G["ctx"]
    extra = ["lall"] if spec["base"] == "b2" else []
    df = table(basin, spec["base"])
    n = len(df)
    eta_base, eta_pc, lams = np.zeros(n), np.zeros(n), []
    for s in ALL:
        train = [x for x in ALL if x != s]
        tr, te = (df.season != s).values, (df.season == s).values
        eof = RG.REOFs(ctx.maps_for("lag1", train), ctx.mask, spec["cells"], spec["vars"])
        pcs_tr = eof.scores(ctx.maps_for("lag1", train))
        pcs_te = eof.scores(ctx.maps_for("lag1", [s]))
        lam, _ = fast_loso(df[tr], extra, pcs_tr)
        mdl = H.fit_model(df[tr], extra, pcs_tr, lam)
        b, p = H.predict(mdl, df[te], pcs_te, parts=True)
        eta_base[te], eta_pc[te] = b, p
        lams.append(lam)
    np.savez(fn, y=df.y.values.astype(float), eta_base=eta_base, eta_pc=eta_pc,
             lam=np.array([np.inf if not np.isfinite(l) else l for l in lams]))
    print(f"{basin} {name} {time.time()-t0:.0f}s lam={lams}", flush=True)
    return fn


def run_refs(basin):
    fn = os.path.join(G["out"], f"refs_{basin}.npz")
    if not os.path.exists(fn):
        y, o = refs(basin)
        np.savez(fn, y=y, **o)
    return fn


def work(job):
    return run_refs(job[1]) if job[0] == "refs" else run_model(job[1])


if __name__ == "__main__":
    wd, out = sys.argv[1], sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    flt = sys.argv[4] if len(sys.argv) > 4 else ""
    os.makedirs(out, exist_ok=True)
    G["out"] = out
    init(wd)
    jobs = [("refs", b) for b in H.BASINS]
    for b in H.BASINS:
        jobs += [("model", s) for s in RG.specs(b) if flt in s["name"]]
    # FULL first so the reproduction check is available early
    jobs.sort(key=lambda j: 0 if j[0] == "refs" or j[1]["name"] == "FULL" else 1)
    with Pool(nproc) as p:
        for fn in p.imap_unordered(work, jobs):
            pass
    with open(os.path.join(out, "looks.log"), "a") as f:
        f.write(f"{dt.datetime.utcnow().isoformat()}Z nested LOSO scoring of {len(jobs)} jobs over archive seasons 2004..2025 "
                f"(filter={flt!r})\n")
