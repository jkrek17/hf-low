"""Hemispheric pattern search: stages. See PREREGISTRATION.md for every choice.

usage: run.py STAGE WORKDIR [OUT]
  freeze   fit every model on the training seasons only (outcomes of the training seasons), choose
           lambda by leave-one-season-out inside training, write frozen models. Does not read the
           outcomes of the test seasons. Splits: primary, swap, lag2, conc (conc = concurrent contrast).
  evaluate score the frozen models once on the test seasons (this is the look at the held-out outcomes).
"""
import os, sys, json, time, pickle, datetime as dt
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hemlib as H
from indices import Indices

NAMED = ["nao", "pna", "oni", "mjo1", "mjo2"]
PROD = ["nao_pna", "nao_oni", "pna_oni"]
WINDOWS = {"lag1": (-7, -1), "lag2": (-14, -8), "conc": (0, 6)}
SPLITS = {
    "primary": dict(train=list(range(2004, 2015)), test=list(range(2015, 2026)), win="lag1"),
    "swap":    dict(train=list(range(2015, 2026)), test=list(range(2004, 2015)), win="lag1"),
    "lag2":    dict(train=list(range(2004, 2015)), test=list(range(2015, 2026)), win="lag2"),
    "conc":    dict(train=list(range(2004, 2015)), test=list(range(2015, 2026)), win="conc"),
}
NPERM = 10000
NBOOT = 10000
NPOW = 2000
RR_LEVELS = [1.05, 1.10, 1.15, 1.20, 1.30, 1.50]


class Ctx:
    def __init__(self, work):
        self.work = work
        t = time.time()
        self.F = H.load_fields(work)
        self.A = H.anomalies(self.F)
        self.mask = H.sst_mask(self.A)
        self.idx = Indices(os.path.join(H.REPO, "docs/data/teleconnections.json"))
        self._maps = {}
        print(f"fields ready in {time.time()-t:.0f}s; sst cells {self.mask.sum()}", flush=True)

    def maps(self, win):
        if win not in self._maps:
            a, b = WINDOWS[win]
            self._maps[win] = H.weekly_maps(self.F, self.A, list(range(1979, 2026)), a, b)
        return self._maps[win]

    def rows(self, seasons):
        return np.array([(s - 1979) * H.NWEEK + k for s in seasons for k in range(H.NWEEK)])

    def maps_for(self, win, seasons):
        r = self.rows(seasons)
        return {v: self.maps(win)[v][r] for v in H.VARS}

    def index_table(self, seasons, win):
        a, b = WINDOWS[win]
        rows = []
        for s in seasons:
            for k in range(H.NWEEK):
                S = H.week_start(s, k)
                rows.append(self.idx.window(S + dt.timedelta(days=a), S + dt.timedelta(days=b)))
        return pd.DataFrame(rows)


def add_indices(df, ix, ref=None):
    """Standardise index columns with the reference (training) mean/sd; returns df and the stats."""
    ref = ix if ref is None else ref
    st = {c: (float(ref[c].mean()), float(ref[c].std(ddof=1))) for c in ix.columns}
    for c in ix.columns:
        df[c] = ((ix[c].values - st[c][0]) / st[c][1])
        df[c] = df[c].fillna(0.0)          # the source has no MJO value in two weeks (2021 wk 13, 2022 wk 13): neutral
    df["nao_pna"] = df.nao * df.pna; df["nao_oni"] = df.nao * df.oni; df["pna_oni"] = df.pna * df.oni
    return df, st


def apply_stats(df, ix, st):
    for c in ix.columns:
        df[c] = ((ix[c].values - st[c][0]) / st[c][1])
        df[c] = df[c].fillna(0.0)          # the source has no MJO value in two weeks (2021 wk 13, 2022 wk 13): neutral
    df["nao_pna"] = df.nao * df.pna; df["nao_oni"] = df.nao * df.oni; df["pna_oni"] = df.pna * df.oni
    return df


# ------------------------------------------------------------------ freeze
def freeze(ctx, name, days):
    sp = SPLITS[name]
    train, win = sp["train"], sp["win"]
    eof = H.EOFs(ctx.maps_for(win, train), ctx.mask)
    pcs_tr = eof.scores(ctx.maps_for(win, train))
    ix_tr = ctx.index_table(train, win)
    out = dict(name=name, train=train, test=sp["test"], win=win, eof=eof, basins={})
    for b in H.BASINS:
        df = H.weekly_table(days, train, b)
        df, st = add_indices(df, ix_tr)
        df0 = df.assign(lprev=0.0)
        m = {}
        m["B0"] = H.fit_model(df0, [], None, np.inf)
        m["B1"] = H.fit_model(df, [], None, np.inf)
        m["N"] = H.fit_model(df, NAMED, None, np.inf)
        m["NI"] = H.fit_model(df, NAMED + PROD, None, np.inf)
        cv = {}
        for key, d, extra in (("P", df, []), ("PN", df, NAMED), ("P0", df0, [])):
            lam, tot = H.loso_lambda(d, extra, pcs_tr)
            cv[key] = {str(k): v for k, v in tot.items()}
            m[key] = H.fit_model(d, extra, pcs_tr, lam)
        # SOM on the Z500 scores
        z = pcs_tr[:, :H.K]
        W = H.som_fit(z)
        node = H.som_assign(W, z)
        nodes = np.eye(9)[node][:, 1:]
        m["SOM"] = H.fit_model(df, [], nodes, 0.5)
        out["basins"][b] = dict(models=m, cv=cv, stats=st, W=W, node_tr=node,
                                dispersion=float(_phi(df, m["B1"])))
    return out


def _phi(df, model):
    mu = H.predict(model, df, None)
    n, p = len(df), len(model["beta"])
    return max(1.0, float(np.sum((df.y.values - mu) ** 2 / mu) / (n - p)))


def write_frozen_json(fz, path):
    """Coefficients, penalties and CV tables of the frozen models, committed before the held-out look."""
    o = dict(split=fz["name"], train=fz["train"], test=fz["test"], window=fz["win"],
             eof_variance_explained={k: round(v, 4) for k, v in fz["eof"].var_explained.items()}, basins={})
    for b, d in fz["basins"].items():
        o["basins"][b] = dict(
            models={k: dict(lam=m["lam"], extra=m["extra"], beta=[round(float(x), 6) for x in m["beta"]]) for k, m in d["models"].items()},
            loso_cv_deviance={k: v for k, v in d["cv"].items()}, index_stats=d["stats"], dispersion=d["dispersion"],
            som_train_node_counts=np.bincount(d["node_tr"], minlength=9).tolist())
    json.dump(o, open(path, "w"), indent=1, default=lambda x: None if x == np.inf else float(x))


def summarize_freeze(fz):
    rows = []
    for b, d in fz["basins"].items():
        r = dict(split=fz["name"], basin=b, n_train_weeks=len(d["node_tr"]), dispersion=round(d["dispersion"], 3))
        for key in ("P", "PN", "P0"):
            cv = d["cv"][key]
            r[f"lambda_{key}"] = d["models"][key]["lam"]
            base = cv["inf"]
            r[f"cvdev_{key}_best_vs_nopc"] = round(1 - min(cv.values()) / base, 4)
        rows.append(r)
    return rows


# ------------------------------------------------------------------ scoring pieces
def parts(model, df, pcs):
    X, nu = H.design(df, model["extra"], pcs if np.isfinite(model["lam"]) else None)
    b = model["beta"]
    base = X[:, :8] @ b[:8]
    extra = X[:, 8:nu] @ b[8:nu] if nu > 8 else np.zeros(len(df))
    pc = X[:, nu:] @ b[nu:] if X.shape[1] > nu else np.zeros(len(df))
    return base, extra, pc


def dev_mat(y, mu):
    t = np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0)
    return 2 * (t - (y - mu)).sum(-1)


def perm_p(y, eta_fixed, eta_blk, D_ref, rng, nseas, kind="ss", n=NPERM):
    """One-sided p for SS = 1 - D(model with permuted block)/D_ref, season-block permutation of eta_blk."""
    ns = nseas
    obs = 1 - H.deviance(y, np.exp(np.clip(eta_fixed + eta_blk, -30, 30))) / D_ref
    P = np.stack([rng.permutation(ns) for _ in range(n)])
    blk = eta_blk.reshape(ns, H.NWEEK)[P].reshape(n, -1)
    mu = np.exp(np.clip(eta_fixed[None] + blk, -30, 30))
    ss = 1 - dev_mat(y[None], mu) / D_ref
    return obs, (1 + int((ss >= obs - 1e-12).sum())) / (1 + n), ss


def boot_delta(y, mus, seasons_idx, rng, n=NBOOT):
    """Bootstrap over seasons of SS(P|B1) - SS(N|B1) = (D_N - D_P)/D_B1."""
    ns = len(seasons_idx)
    dP = np.array([H.deviance(y[i], mus["P"][i]) for i in seasons_idx])
    dN = np.array([H.deviance(y[i], mus["N"][i]) for i in seasons_idx])
    d1 = np.array([H.deviance(y[i], mus["B1"][i]) for i in seasons_idx])
    obs = (dN.sum() - dP.sum()) / d1.sum()
    draws = rng.integers(0, ns, size=(n, ns))
    bs = (dN[draws].sum(1) - dP[draws].sum(1)) / d1[draws].sum(1)
    return float(obs), (float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))), (1 + int((bs <= 0).sum())) / (1 + n)


# ------------------------------------------------------------------ evaluate
def evaluate(ctx, fz, days, log):
    sp = SPLITS[fz["name"]]
    test, win = fz["test"], fz["win"]
    eof = fz["eof"]
    pcs_te = eof.scores(ctx.maps_for(win, test))
    ix_te = ctx.index_table(test, win)
    res = {}
    rng = np.random.default_rng(H.SEED)
    for b in H.BASINS:
        d = fz["basins"][b]
        m = d["models"]
        df = H.weekly_table(days, test, b)
        df = apply_stats(df, ix_te, d["stats"])
        df0 = df.assign(lprev=0.0)
        y = df.y.values.astype(float)
        ns = len(test)
        sidx = [np.arange(i * H.NWEEK, (i + 1) * H.NWEEK) for i in range(ns)]
        z = pcs_te[:, :H.K]
        node = H.som_assign(d["W"], z)
        nodes = np.eye(9)[node][:, 1:]
        mu = {"B0": H.predict(m["B0"], df0, None), "B1": H.predict(m["B1"], df, None),
              "N": H.predict(m["N"], df, None), "NI": H.predict(m["NI"], df, None),
              "P": H.predict(m["P"], df, pcs_te), "PN": H.predict(m["PN"], df, pcs_te),
              "P0": H.predict(m["P0"], df0, pcs_te), "SOM": H.predict(m["SOM"], df, nodes)}
        D = {k: H.deviance(y, v) for k, v in mu.items()}
        r = dict(basin=b, n_weeks=len(y), n_seasons=ns, total_events=int(y.sum()), D=D,
                 lam=dict(P=m["P"]["lam"], PN=m["PN"]["lam"], P0=m["P0"]["lam"]))
        ss = lambda x, ref: 1 - D[x] / D[ref]
        r["SS"] = {"P|B1": ss("P", "B1"), "N|B1": ss("N", "B1"), "NI|B1": ss("NI", "B1"), "PN|N": ss("PN", "N"),
                   "P|B0": ss("P0", "B0"), "B1|B0": ss("B1", "B0"), "SOM|B1": ss("SOM", "B1"), "PN|B1": ss("PN", "B1")}
        # P1
        bP, eP, pP = parts(m["P"], df, pcs_te)
        obs, p1, _ = perm_p(y, bP + eP, pP, D["B1"], rng, ns)
        r["P1"] = dict(ss=obs, p=p1)
        # P4
        bPN, ePN, pPN = parts(m["PN"], df, pcs_te)
        obs, p4, _ = perm_p(y, bPN + ePN, pPN, D["N"], rng, ns)
        r["P4"] = dict(ss=obs, p=p4)
        # P3
        o, ci, p3 = boot_delta(y, mu, sidx, rng)
        r["P3"] = dict(delta=o, ci=ci, p=p3)
        # informational: named model vs B1 (permute the index rows)
        bN, eN, _ = parts(m["N"], df, None)
        obs, pn, _ = perm_p(y, bN, eN, D["B1"], rng, ns)
        r["named_perm"] = dict(ss=obs, p=pn)
        # SOM
        bS, eS, pS = parts(m["SOM"], df, nodes)
        obs, ps, _ = perm_p(y, bS + eS, pS, D["B1"], rng, ns)
        r["SOM"] = dict(ss=obs, p=ps)
        # SOM rate ratios (held-out and discovery), node counts
        r["som_nodes"] = dict(train_n=np.bincount(d["node_tr"], minlength=9).tolist(), test_n=np.bincount(node, minlength=9).tolist(),
                              beta=np.concatenate([[0.0], m["SOM"]["beta"][-8:]]).tolist(),
                              test_rate=[float(y[node == j].mean()) if (node == j).any() else None for j in range(9)])
        # power: planted effect on the frozen pattern index
        lam_f = m["P"]["lam"]
        r["power"] = power(y, df, pcs_te, m, d, bP + eP, pP, mu, rng, test)
        # per-season detail for reproduction
        r["per_season_dev"] = dict(B1=[H.deviance(y[i], mu["B1"][i]) for i in sidx], P=[H.deviance(y[i], mu["P"][i]) for i in sidx],
                                   N=[H.deviance(y[i], mu["N"][i]) for i in sidx])
        res[b] = r
        log(f"{fz['name']} {b}: SS(P|B1)={r['SS']['P|B1']:.4f} p={p1:.4f}  SS(N|B1)={r['SS']['N|B1']:.4f}  "
            f"SS(PN|N)={r['SS']['PN|N']:.4f} p={p4:.4f}  P3 delta={r['P3']['delta']:.4f} ci={ci}")
    return res


def power(y, df, pcs, m, d, eta_P, eta_pc, mu, rng, test):
    """Planted-effect power for P1. z = standardised pattern index of the lambda chosen (or best finite lambda)."""
    mP = m["P"]
    if not np.isfinite(mP["lam"]):
        fin = {float(k): v for k, v in d["cv"]["P"].items() if k != "inf"}
        lam = min(fin, key=fin.get)
        return dict(note=f"lambda=none chosen; planted on best finite lambda {lam}", rr=[], power=[], detectable=None, null_crit=None)
    z = (eta_pc - eta_pc.mean()) / eta_pc.std()
    phi = d["dispersion"]
    n = len(y)
    base1 = mu["B1"]
    muP0 = mu["P"]; mu1 = mu["B1"]
    def draw(rr, nsim):
        mean = np.repeat((base1 * (rr ** z))[None], nsim, axis=0)
        if phi > 1.0001:
            k = mean / (phi - 1.0)
            return rng.poisson(rng.gamma(k, mean / k))
        return rng.poisson(mean)
    def stat(ysim):
        ysim = ysim.astype(float)
        return 1 - dev_mat(ysim, muP0[None]) / dev_mat(ysim, mu1[None])
    null = stat(draw(1.0, NPOW))
    crit = float(np.percentile(null, 95))
    pw = []
    for rr in RR_LEVELS:
        s = stat(draw(rr, NPOW))
        pw.append(float((s > crit).mean()))
    det = next((rr for rr, p in zip(RR_LEVELS, pw) if p >= 0.8), None)
    return dict(rr=RR_LEVELS, power=pw, detectable=det, null_crit=crit, dispersion=phi)



# ------------------------------------------------------------------ S4 proxy transfer, S5 within-era
def proxy_transfer(ctx, fz, days_hf, log):
    """S4: the archive-fitted primary models scored on pipeline A HF-equivalent counts in the test seasons."""
    test, win = fz["test"], fz["win"]
    pcs_te = fz["eof"].scores(ctx.maps_for(win, test))
    ix_te = ctx.index_table(test, win)
    rng = np.random.default_rng(H.SEED + 4)
    out = {}
    for b in H.BASINS:
        d = fz["basins"][b]; m = d["models"]
        df = apply_stats(H.weekly_table(days_hf, test, b), ix_te, d["stats"])
        y = df.y.values.astype(float)
        mu1 = H.predict(m["B1"], df, None); muP = H.predict(m["P"], df, pcs_te)
        D1, DP = H.deviance(y, mu1), H.deviance(y, muP)
        base, extra, pc = parts(m["P"], df, pcs_te)
        obs, p, _ = perm_p(y, base + extra, pc, D1, rng, len(test))
        out[b] = dict(ss=1 - DP / D1, p=p, total_events=int(y.sum()))
        log(f"S4 proxy {b}: SS(P|B1)={1-DP/D1:.4f} p={p:.4f} events={int(y.sum())}")
    return out


def within_era(ctx, fz, days_depth, log, nperm=5000):
    """S5: 1979-80..2000-01, pipeline A depth counts. month + previous week + gamma * frozen pattern index."""
    seasons = list(range(1979, 2001))
    pcs = fz["eof"].scores(ctx.maps_for(fz["win"], seasons))
    rng = np.random.default_rng(H.SEED + 5)
    out = {}
    for b in H.BASINS:
        m = fz["basins"][b]["models"]["P"]
        df = H.weekly_table(days_depth, seasons, b)
        if not np.isfinite(m["lam"]):
            out[b] = dict(note="no pattern retained"); continue
        X, nu = H.design(df.assign(), [], pcs)
        idx = X[:, nu:] @ m["beta"][nu:]
        idx = (idx - idx.mean()) / idx.std()
        y = df.y.values.astype(float)
        def gain(ix):
            d1 = df.assign(ix=ix)
            m0 = H.fit_model(d1, [], None, np.inf)
            m1 = H.fit_model(d1, ["ix"], None, np.inf)
            return H.deviance(y, H.predict(m0, d1, None)) - H.deviance(y, H.predict(m1, d1, None)), float(m1["beta"][8])
        g0, gamma = gain(idx)
        ns = len(seasons)
        null = []
        for _ in range(nperm):
            pm = rng.permutation(ns)
            null.append(gain(idx.reshape(ns, H.NWEEK)[pm].reshape(-1))[1])
        null = np.array(null)
        p = (1 + int((null >= gamma).sum())) / (1 + nperm)
        out[b] = dict(gamma=gamma, rr_per_sd=float(np.exp(gamma)), deviance_gain=g0, p=p, total_events=int(y.sum()))
        log(f"S5 within-era {b}: gamma={gamma:.4f} RR/SD={np.exp(gamma):.3f} p={p:.4f} events={int(y.sum())}")
    return out


def main():
    stage, work = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(H.HERE, "results")
    os.makedirs(out, exist_ok=True)
    ctx = Ctx(work)
    if stage == "freeze":
        names = sys.argv[4].split(",") if len(sys.argv) > 4 else ["primary", "swap", "lag2", "conc"]
        days = H.load_archive_counts(None)
        allrows = []
        for nm in names:
            t = time.time()
            # training outcomes only: test-season rows are never built here
            fz = freeze(ctx, nm, days)
            pickle.dump(fz, open(os.path.join(ctx.work, f"frozen_{nm}.pkl"), "wb"))
            write_frozen_json(fz, os.path.join(out, f"frozen_{nm}.json"))
            rows = summarize_freeze(fz)
            allrows += rows
            print(nm, f"{time.time()-t:.0f}s", rows, flush=True)
        pd.DataFrame(allrows).to_csv(os.path.join(out, "freeze_" + "_".join(names) + ".csv"), index=False)
    elif stage == "evaluate":
        names = sys.argv[4].split(",")
        days = H.load_archive_counts(None)
        lines = []
        def log(x):
            print(x, flush=True); lines.append(x)
        for nm in names:
            fz = pickle.load(open(os.path.join(ctx.work, f"frozen_{nm}.pkl"), "rb"))
            res = evaluate(ctx, fz, days, log)
            if nm == "primary":
                res["S4_proxy"] = proxy_transfer(ctx, fz, H.load_proxy_days("hf"), log)
                res["S5_within_era"] = within_era(ctx, fz, H.load_proxy_days("depth"), log)
            json.dump(res, open(os.path.join(out, f"heldout_{nm}.json"), "w"), indent=1, default=float)
            with open(os.path.join(out, "heldout_looks.log"), "a") as f:
                f.write(f"{dt.datetime.utcnow().isoformat()}Z scored split={nm} test={fz['test'][0]}..{fz['test'][-1]} window={fz['win']}\n")
        open(os.path.join(out, "heldout_" + "_".join(names) + ".txt"), "w").write("\n".join(lines) + "\n")
    else:
        raise SystemExit("unknown stage")


if __name__ == "__main__":
    main()
