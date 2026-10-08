"""RA-11: do archive HF fixes per event drift, or step at a scatterometer change? Plan: PREREGISTRATION.md.

Archive HF events (basin, ID), seasons 2004-05..2025-26. ERA5 pipeline A appears only as a proxy contrast (S7).
Writes results/fixes_per_event-result.txt and results/tests.csv.   usage: python3 analyse.py
"""
import csv, io, math, os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
ARCH = os.path.join(ROOT, "data", "hf_lows")
LIFE = os.path.join(HERE, "..", "hf_history", "results", "lifecycle_events.csv")
LINKS = os.path.join(HERE, "..", "count_reconcile", "results", "linked_ids_for_review.csv")
FIRST, LAST = 2004, 2025
SEASONS = np.arange(FIRST, LAST + 1)
STEPS = [("P2 ASCAT-A launch 2006-10", "2006-10-01"), ("P3 QuikSCAT end 2009-11-23", "2009-11-23"),
         ("P4 ASCAT-B launch 2012-09", "2012-09-01"), ("P5 Metop-C launch 2018-11", "2018-11-01")]
B = 2000
SEED = 11


def t_cdf(t, df):
    if t == 0:
        return 0.5
    c = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
    f = lambda x: c * (1 + x * x / df) ** (-(df + 1) / 2)
    n = 4000
    h = abs(t) / n
    s = f(0) + f(abs(t)) + sum((4 if i % 2 else 2) * f(i * h) for i in range(1, n))
    a = s * h / 3
    return 0.5 + a if t > 0 else 0.5 - a


def t_p(t, df):
    return 2 * (1 - t_cdf(abs(t), df))


def t_crit(df):
    lo, hi = 0.0, 20.0
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if t_cdf(m, df) < 0.975 else (lo, m)
    return (lo + hi) / 2


def ols(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x); xm = x - x.mean(); sxx = (xm ** 2).sum()
    b = (xm * (y - y.mean())).sum() / sxx
    r = y - y.mean() - b * xm
    se = math.sqrt((r ** 2).sum() / (n - 2) / sxx)
    t = b / se
    dw = (np.diff(r) ** 2).sum() / (r ** 2).sum()
    return b, se, t, t_p(t, n - 2), dw, r


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    prev = 1.0
    for rank in range(n, 0, -1):
        i = o[rank - 1]
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q


# ------------------------------------------------------------------ data
def load_archive(merge=False, dedup=False):
    rows = []
    for b in ("Atl", "Pac"):
        d = pd.read_csv(os.path.join(ARCH, "HF_Data_-_%s.csv" % b), dtype=str)
        d["basin"] = b.lower()
        rows.append(d)
    d = pd.concat(rows, ignore_index=True)
    d["ID"] = d.ID.str.strip()
    d["key"] = d.basin + ":" + d.ID
    d = d[d.Category.str.strip() == "HF"].copy()
    # one HF row (pac:2024202506) has an 11-digit stamp 20241101018 between 12Z and 00Z fixes; read as 2024110118
    d["date"] = d.date.where(d.date.str.len() == 10, d.date.str[:8] + d.date.str[9:])
    d["t"] = pd.to_datetime(d.date, format="%Y%m%d%H", errors="coerce")
    assert d.t.notna().all(), d[d.t.isna()]
    d["lat"] = pd.to_numeric(d.Latitude, errors="coerce")
    d["p"] = pd.to_numeric(d.Pressure, errors="coerce")
    if merge:
        lk = pd.read_csv(LINKS)
        m = dict(zip(lk.id_b, lk.id_a))
        # chain-safe: resolve a->b->c
        def res(k):
            seen = set()
            while k in m and k not in seen:
                seen.add(k); k = m[k]
            return k
        d["key"] = d.key.map(res)
        dedup = True
    if dedup:
        d = d.drop_duplicates(["key", "t"])
    g = d.sort_values("t").groupby("key")
    ev = g.agg(n_hf=("t", "size"), t0=("t", "min"), minp=("p", "min"), lat0=("lat", "first"),
               basin=("basin", "first")).reset_index()
    ev["season"] = ev.key.str.split(":").str[1].str[:4].astype(int)
    ev = ev[ev.season.between(FIRST, LAST)].reset_index(drop=True)
    return ev


def season_means(ev, col="n_hf"):
    s = ev.groupby("season")[col].mean().reindex(SEASONS)
    return s.values


def step_stats(ev, date):
    """per-season sums/counts before and after the date (event first HF fix)."""
    after = (ev.t0 >= pd.Timestamp(date)).values
    S = {}
    for nm, msk in (("b", ~after), ("a", after)):
        sub = ev[msk]
        S["S" + nm] = sub.groupby("season").n_hf.sum().reindex(SEASONS, fill_value=0).values.astype(float)
        S["N" + nm] = sub.groupby("season").n_hf.size().reindex(SEASONS, fill_value=0).values.astype(float)
    return S


def step_diff(S, w=None):
    if w is None:
        w = np.ones(len(SEASONS))
    w = np.atleast_2d(w)
    Na, Nb = w @ S["Na"], w @ S["Nb"]
    with np.errstate(divide="ignore", invalid="ignore"):
        d = (w @ S["Sa"]) / Na - (w @ S["Sb"]) / Nb
    return d


def boot_weights(rng, nb=B):
    idx = rng.integers(0, len(SEASONS), size=(nb, len(SEASONS)))
    W = np.zeros((nb, len(SEASONS)))
    for i in range(nb):
        W[i] = np.bincount(idx[i], minlength=len(SEASONS))
    return W


def boot_p(est, boots):
    boots = boots[np.isfinite(boots)]
    pl = (np.sum(boots <= 0) + 1) / (len(boots) + 1)
    pg = (np.sum(boots >= 0) + 1) / (len(boots) + 1)
    return min(1.0, 2 * min(pl, pg)), np.percentile(boots, 2.5), np.percentile(boots, 97.5), len(boots)


def frac_after(ev, date):
    a = (ev.t0 >= pd.Timestamp(date)).astype(float)
    return a.groupby(ev.season).mean().reindex(SEASONS).values


def aic(y, X):
    n = len(y)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    rss = ((y - X @ beta) ** 2).sum()
    return n * math.log(rss / n) + 2 * X.shape[1]


def form_table(y, x_year, fr):
    """AIC for null, linear, step with season-mean y and step fraction fr."""
    one = np.ones(len(y))
    a0 = aic(y, one[:, None])
    a1 = aic(y, np.c_[one, x_year])
    a2 = aic(y, np.c_[one, fr])
    return a0, a1, a2


def prefer(a1, a2):
    if a2 <= a1 - 2: return "step"
    if a1 <= a2 - 2: return "linear"
    return "tie"


# ------------------------------------------------------------------ main
def run(ev, label, P, full=True):
    """primary tests on ev; returns dict. P = print function."""
    out = {}
    y = season_means(ev)
    P("%s: %d events, %d seasons, mean n_hf %.3f" % (label, len(ev), len(SEASONS), ev.n_hf.mean()))
    b, se, t, p, dw, _ = ols(SEASONS, y)
    P("  P1 linear: %+.3f /decade (se %.3f) t=%+.2f p=%.4f DW=%.2f" % (b * 10, se * 10, t, p, dw))
    out["P1"] = (b * 10, p)
    rng = np.random.default_rng(SEED)
    W = boot_weights(rng)
    # event-level bootstrap interval for the P1 slope: resample seasons (pairs)
    idx = rng.integers(0, len(SEASONS), size=(B, len(SEASONS)))
    sl = []
    for r in idx:
        xs = SEASONS[r]
        if len(set(xs)) < 3: continue
        sl.append(ols(xs, y[r])[0] * 10)
    sl = np.array(sl)
    P("     season-bootstrap 95%% interval for the P1 slope: [%+.3f, %+.3f]" % tuple(np.percentile(sl, [2.5, 97.5])))
    out["P1_ci"] = tuple(np.percentile(sl, [2.5, 97.5]))
    for nm, dt in STEPS:
        S = step_stats(ev, dt)
        est = float(step_diff(S)[0])
        bt = step_diff(S, W)
        pb, lo, hi, nbok = boot_p(est, bt)
        nb_, na_ = int(S["Nb"].sum()), int(S["Na"].sum())
        sb = int((S["Nb"] > 0).sum()); sa = int((S["Na"] > 0).sum())
        P("  %-28s step %+.3f [%+.3f, %+.3f] p=%.4f  (events %d before / %d after; seasons with events %d / %d; %d usable draws)" % (
            nm, est, lo, hi, pb, nb_, na_, sb, sa, nbok))
        out[nm.split()[0]] = (est, pb, lo, hi)
    return out


def main():
    buf = io.StringIO()

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True); buf.write(s + "\n")

    P("RA-11: archive HF fixes per event, drift vs scatterometer steps (analyse.py)")
    P("=" * 80)
    P("Archive HF events (basin, ID), seasons 2004-05..2025-26 (June-May, 22 seasons). ERA5 is a PROXY and")
    P("appears only in S7 (pipeline A = research/era5/hf_history). Plan: PREREGISTRATION.md.")
    P("")
    ev = load_archive()
    allp = []  # (name, family, p)

    P("== PRIMARY (family P, primary dataset: raw IDs) ==")
    prim = run(ev, "raw IDs", P)
    allp.append(("P1 linear drift", "P", prim["P1"][1]))
    for nm, _ in STEPS:
        allp.append((nm, "P", prim[nm.split()[0]][1]))
    qP = bh([a[2] for a in allp[:5]])
    P("  BH q within P: " + ", ".join("%s q=%.4f" % (a[0].split()[0], q) for a, q in zip(allp[:5], qP)))
    P("")

    # ---------------- S1 scan
    y = season_means(ev)
    P("== S1 unknown-date scan (season-mean step, split at season boundaries, >=3 seasons each side) ==")
    cands = list(range(2007, 2023))  # first season of the 'after' block
    def scan(yv):
        best = (0, None)
        for c in cands:
            k = c - FIRST
            if k < 3 or len(SEASONS) - k < 3: continue
            d = yv[k:].mean() - yv[:k].mean()
            s2 = (((yv[:k] - yv[:k].mean()) ** 2).sum() + ((yv[k:] - yv[k:].mean()) ** 2).sum()) / (len(yv) - 2)
            tt = d / math.sqrt(s2 * (1 / k + 1 / (len(yv) - k)))
            if abs(tt) > abs(best[0]): best = (tt, c)
        return best
    tobs, cbest = scan(y)
    rng = np.random.default_rng(SEED)
    perm = np.array([abs(scan(rng.permutation(y))[0]) for _ in range(5000)])
    ps1 = (np.sum(perm >= abs(tobs)) + 1) / (len(perm) + 1)
    k = cbest - FIRST
    P("  best split: seasons %d.. vs <%d; step %+.3f, t=%+.2f; permutation p (max|t| over %d splits)=%.4f" % (
        cbest, cbest, y[k:].mean() - y[:k].mean(), tobs, len(cands), ps1))
    allp.append(("S1 scan best step", "S", ps1))
    P("")

    # ---------------- S2 form comparison
    P("== S2 form comparison (AIC on season means; step regressor = fraction of the season's events after the date) ==")
    P("   rule: lower AIC by >=2 preferred, else tie")
    a0, a1, _ = form_table(y, SEASONS, np.zeros(len(y)))
    P("  null AIC %.2f, linear AIC %.2f" % (a0, a1))
    forms = {}
    for nm, dt in STEPS:
        fr = frac_after(ev, dt)
        _, _, a2 = form_table(y, SEASONS, fr)
        forms[nm] = a2
        P("  %-28s step AIC %.2f  (step - linear = %+.2f) -> %s" % (nm, a2, a2 - a1, prefer(a1, a2)))
    frb = (SEASONS >= cbest).astype(float)
    _, _, a2 = form_table(y, SEASONS, frb)
    P("  S1 best step (%d.., chosen on this data, so optimistic) step AIC %.2f (step - linear = %+.2f) -> %s" % (
        cbest, a2, a2 - a1, prefer(a1, a2)))
    P("")

    # ---------------- S3 basins
    P("== S3 basins ==")
    rng = np.random.default_rng(SEED); W = boot_weights(rng)
    for bs in ("atl", "pac"):
        eb = ev[ev.basin == bs]
        yb = season_means(eb)
        m = np.isfinite(yb)
        b_, se_, t_, p_, dw_, _ = ols(SEASONS[m], yb[m])
        P("  %s P1 linear %+.3f /decade (se %.3f) t=%+.2f p=%.4f  (n events %d, mean %.3f, seasons %d)" % (
            bs, b_ * 10, se_ * 10, t_, p_, len(eb), eb.n_hf.mean(), m.sum()))
        allp.append(("S3 %s linear" % bs, "S", p_))
        S = step_stats(eb, "2009-11-23"); est = float(step_diff(S)[0])
        pb, lo, hi, _ = boot_p(est, step_diff(S, W))
        P("  %s P3 QuikSCAT step %+.3f [%+.3f, %+.3f] p=%.4f" % (bs, est, lo, hi, pb))
        allp.append(("S3 %s QuikSCAT step" % bs, "S", pb))
    P("")

    # ---------------- S4 shape
    P("== S4 shape: season share of one-fix and >=3-fix events ==")
    for nm, f in (("one-fix share", lambda d: (d.n_hf == 1).mean()), (">=3-fix share", lambda d: (d.n_hf >= 3).mean())):
        ys = ev.groupby("season").apply(f, include_groups=False).reindex(SEASONS).values
        b_, se_, t_, p_, _, _ = ols(SEASONS, ys)
        P("  %-14s mean %.3f, slope %+.4f /decade (se %.4f) t=%+.2f p=%.4f" % (nm, ys.mean(), b_ * 10, se_ * 10, t_, p_))
        allp.append(("S4 " + nm, "S", p_))
    P("")

    # ---------------- S5 Atlantic >60N
    P("== S5 Atlantic first HF fix north vs south of 60N: QuikSCAT step ==")
    ea = ev[ev.basin == "atl"]
    rng = np.random.default_rng(SEED); W = boot_weights(rng)
    est = {}
    bts = {}
    for nm, sub in (("north", ea[ea.lat0 >= 60]), ("south", ea[ea.lat0 < 60])):
        S = step_stats(sub, "2009-11-23")
        est[nm] = float(step_diff(S)[0]); bts[nm] = step_diff(S, W)
        pb, lo, hi, _ = boot_p(est[nm], bts[nm])
        P("  %s of 60N: events %d (before %d / after %d), step %+.3f [%+.3f, %+.3f] p=%.4f" % (
            nm, len(sub), S["Nb"].sum(), S["Na"].sum(), est[nm], lo, hi, pb))
        allp.append(("S5 Atl %s step" % nm, "S", pb))
    dif = est["north"] - est["south"]; bd = bts["north"] - bts["south"]
    pb, lo, hi, _ = boot_p(dif, bd)
    P("  north minus south step %+.3f [%+.3f, %+.3f] p=%.4f" % (dif, lo, hi, pb))
    allp.append(("S5 north-south difference", "S", pb))
    P("")

    # ---------------- S6 covariates
    P("== S6 covariates: slope of n_hf per decade with event minimum pressure, basin and first-fix month ==")
    e6 = ev.dropna(subset=["minp"]).copy()
    e6["mon"] = e6.t0.dt.month
    def design(d, adj):
        cols = [np.ones(len(d)), ((d.season.values - 2014.5) / 10.0)]
        if adj:
            cols += [(d.minp.values - 960) / 10.0, (d.basin.values == "pac").astype(float)]
            for m in sorted(range(1, 13)):
                if m != 1:
                    cols.append((d.mon.values == m).astype(float))
        return np.column_stack(cols)
    def slope(d, adj):
        X = design(d, adj); beta, *_ = np.linalg.lstsq(X, d.n_hf.values.astype(float), rcond=None)
        return beta[1]
    rng = np.random.default_rng(SEED)
    res = {}
    for adj in (False, True):
        s0 = slope(e6, adj)
        bs_ = []
        groups = {s: e6[e6.season == s] for s in SEASONS}
        for _ in range(1000):
            pick = rng.integers(0, len(SEASONS), len(SEASONS))
            d = pd.concat([groups[SEASONS[i]] for i in pick])
            bs_.append(slope(d, adj))
        bs_ = np.array(bs_)
        # centre test on 0 using percentile-p
        pb, lo, hi, _ = boot_p(s0, bs_)
        res[adj] = (s0, lo, hi, pb)
        P("  event-level slope %s: %+.3f /decade, season-bootstrap 95%% [%+.3f, %+.3f], p=%.4f (n=%d events with pressure)" % (
            "adjusted  " if adj else "unadjusted", s0, lo, hi, pb, len(e6)))
    allp.append(("S6 adjusted slope", "S", res[True][3]))
    ym = e6.groupby("season").minp.mean().reindex(SEASONS).values
    b_, se_, t_, p_, _, _ = ols(SEASONS, ym)
    P("  event minimum pressure trend: %+.2f hPa/decade (se %.2f) t=%+.2f p=%.4f (mean %.1f hPa)" % (b_ * 10, se_ * 10, t_, p_, e6.minp.mean()))
    allp.append(("S6 minp trend", "S", p_))
    P("")

    # ---------------- S7 proxy
    P("== S7 ERA5 pipeline A proxy (HF-equivalent tracks, n_hf = in-domain 6-hourly fixes >=71.7 kt gust) ==")
    L = pd.read_csv(LIFE)
    L = L[L.season.between(FIRST, LAST)]
    yp = L.groupby("season").n_hf.mean().reindex(SEASONS).values
    b_, se_, t_, p_, dw_, _ = ols(SEASONS, yp)
    P("  proxy fixes per event: mean %.3f (%d tracks), slope %+.3f /decade (se %.3f) t=%+.2f p=%.4f" % (
        L.n_hf.mean(), len(L), b_ * 10, se_ * 10, t_, p_))
    allp.append(("S7 proxy slope", "S", p_))
    rng = np.random.default_rng(SEED)
    diffs = []
    for _ in range(B):
        r = rng.integers(0, len(SEASONS), len(SEASONS))
        if len(set(r)) < 3: continue
        diffs.append((ols(SEASONS[r], y[r])[0] - ols(SEASONS[r], yp[r])[0]) * 10)
    d0 = (ols(SEASONS, y)[0] - ols(SEASONS, yp)[0]) * 10
    pb, lo, hi, _ = boot_p(d0, np.array(diffs))
    P("  archive slope minus proxy slope %+.3f /decade [%+.3f, %+.3f], p=%.4f" % (d0, lo, hi, pb))
    allp.append(("S7 archive-proxy difference", "S", pb))
    P("  proxy QuikSCAT step (event unit = track season, no dates here): not computed (the table has no event date).")
    P("")

    # ---------------- post hoc (labelled; not in any family)
    P("== POST HOC (not pre-registered; not in the FDR family) ==")
    P("  season means of n_hf (events): " + "  ".join("%d:%.2f(%d)" % (s_, v, int(c)) for s_, v, c in
        zip(SEASONS, y, ev.groupby("season").size().reindex(SEASONS).values)))
    k = 2009 - FIRST + 1  # seasons 2004-2009 vs 2010-2025 (QuikSCAT ended inside season 2009)
    for nm, sl_ in (("2004-2009", slice(0, k)), ("2010-2025", slice(k, None))):
        b_, se_, t_, p_, _, _ = ols(SEASONS[sl_], y[sl_])
        P("  within %s: slope %+.3f /decade (se %.3f) t=%+.2f p=%.3f, mean %.3f" % (nm, b_ * 10, se_ * 10, t_, p_, y[sl_].mean()))
    fr = frac_after(ev, "2009-11-23")
    X = np.c_[np.ones(len(y)), SEASONS - 2014.5, fr]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res_ = y - X @ beta; s2_ = (res_ ** 2).sum() / (len(y) - 3)
    cov = s2_ * np.linalg.inv(X.T @ X)
    tl, ts = beta[1] / math.sqrt(cov[1, 1]), beta[2] / math.sqrt(cov[2, 2])
    P("  linear + QuikSCAT step together: linear %+.3f /decade (t=%+.2f, p=%.3f), step %+.3f (t=%+.2f, p=%.3f), 19 df" % (
        beta[1] * 10, tl, t_p(tl, 19), beta[2], ts, t_p(ts, 19)))
    P("  step at QuikSCAT vs step at ASCAT-B launch together: see AIC lines in S2 (P3 %.2f vs P4 %.2f)" % (forms[STEPS[1][0]], forms[STEPS[2][0]]))
    P("")

    # ---------------- sensitivities
    P("== SENSITIVITY: 21 possible duplicate archive IDs merged (not in FDR family) ==")
    evm = load_archive(merge=True)
    P("  events %d (raw %d); merged %d pairs" % (len(evm), len(ev), len(ev) - len(evm)))
    sm = run(evm, "merged IDs", P)
    P("")
    P("== SENSITIVITY: identical timestamps within an ID counted once (not in FDR family) ==")
    evd = load_archive(dedup=True)
    sd = run(evd, "dedup timestamps", P)
    P("")

    # ---------------- multiplicity
    P("== MULTIPLICITY ==")
    names = [a[0] for a in allp]; fam = [a[1] for a in allp]; pv = np.array([a[2] for a in allp])
    qa = bh(pv)
    qs = np.full(len(pv), np.nan)
    for f in ("P", "S"):
        ix = [i for i, x in enumerate(fam) if x == f]
        qs[ix] = bh(pv[ix])
    rows = []
    for n, f, p, qf, qall in zip(names, fam, pv, qs, qa):
        P("  %-34s fam %s  p=%.4f  q(within family)=%.4f  q(all %d)=%.4f" % (n, f, p, qf, len(pv), qall))
        rows.append(dict(test=n, family=f, p=p, q_family=qf, q_all=qall))
    P("  passing q<0.05 across all: %d of %d" % (int((qa < 0.05).sum()), len(pv)))
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "results", "tests.csv"), index=False)
    open(os.path.join(HERE, "results", "fixes_per_event-result.txt"), "w").write(buf.getvalue())


if __name__ == "__main__":
    main()
