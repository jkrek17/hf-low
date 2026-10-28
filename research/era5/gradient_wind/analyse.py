"""Tables (A all cells, C storm-core cells) and the three pre-registered tests. Usage: python3 analyse.py PAIRS.npz OUTDIR"""
import os, sys
import numpy as np

BINS = [0, 2, 4, 6, 8, 10, 12, 99]
LATB = [("30-37.5N", 30.0, 37.5), ("37.5-45N", 37.5, 45.01)]
STAB = [("stable (<0 C)", -99, 0), ("weakly unstable (0 to <6 C)", 0, 6), ("strongly unstable (6 C+)", 6, 99)]
CUT = [38.0, 54.0, 71.7]
B = 2000
rng = np.random.default_rng(20261028)


def winboot(win, num, den):
    """winter-bootstrap 95% interval of sum(num)/sum(den)"""
    u = np.unique(win); n_ = np.array([num[win == w].sum() for w in u], float); d_ = np.array([den[win == w].sum() for w in u], float)
    ix = rng.integers(0, len(u), (B, len(u)))
    r = n_[ix].sum(1) / np.maximum(d_[ix].sum(1), 1)
    return n_.sum() / d_.sum(), np.percentile(r, [2.5, 97.5])


def row(d, m):
    n = int(m.sum()); w = d["winter"][m]
    if n < 200 or len(np.unique(w)) < 8: return None
    g = d["gust"][m]; out = dict(n=n, winters=len(np.unique(w)), med=np.median(g), p10=np.percentile(g, 10), p90=np.percentile(g, 90),
                                 ratio=np.median(g / d["G"][m]))
    for c in CUT:
        p, ci = winboot(w, (g >= c).astype(float), np.ones(n)); out[f"p{c}"] = (p, ci)
    return out


def fmt(r):
    if r is None: return "thin"
    s = f"n={r['n']:>6} w={r['winters']:>2} med={r['med']:5.1f} p10-p90={r['p10']:4.1f}-{r['p90']:4.1f} gust/G={r['ratio']:.2f}"
    for c in CUT:
        p, ci = r[f"p{c}"]; s += f" | >={c:g}: {100*p:5.1f}% ({100*ci[0]:.1f}-{100*ci[1]:.1f})"
    return s


def table(d, name, extra):
    lines = [f"== {name} =="]
    for lb, lo, hi in LATB:
        for sn, slo, shi in ([("all stability", -999, 999)] + STAB):
            lines.append(f"-- {lb}, {sn}")
            for a, b in zip(BINS[:-1], BINS[1:]):
                m = extra & (d["lat"] >= lo) & (d["lat"] < hi) & (d["grad"] >= a) & (d["grad"] < b)
                if sn != "all stability": m &= np.isfinite(d["dT"]) & (d["dT"] >= slo) & (d["dT"] < shi)
                lines.append(f"  {a:>2}-{b if b < 99 else '+':<3} hPa/100nm  " + fmt(row(d, m)))
    return lines


def wmed(d, m):
    """winter-level medians of gust/G -> dict winter -> value"""
    r = d["gust"][m] / d["G"][m]; w = d["winter"][m]
    return {int(x): float(np.median(r[w == x])) for x in np.unique(w) if (w == x).sum() >= 20}


def diff_boot(a, b):
    ws = sorted(set(a) & set(b)); x = np.array([a[w] - b[w] for w in ws])
    ix = rng.integers(0, len(x), (B, len(x))); bs = x[ix].mean(1)
    return x.mean(), np.percentile(bs, [2.5, 97.5]), len(ws), float(2 * min((bs <= 0).mean(), (bs >= 0).mean()))


def ranks(a): return np.argsort(np.argsort(a)).astype(float)


def main(pairs, out):
    d = dict(np.load(pairs)); os.makedirs(out, exist_ok=True)
    ok = np.isfinite(d["gust"]) & np.isfinite(d["grad"])
    d = {k: v[ok] for k, v in d.items()}
    L = [f"cell-times {len(d['gust'])}, winters {len(np.unique(d['winter']))} ({d['winter'].min()}-{d['winter'].max()}), with temperatures {int(np.isfinite(d['dT']).sum())}",
         f"core cell-times {int(d['core'].sum())}; share with no real gradient-wind root: core, grad>=4: "
         f"{100*np.isnan(d['V'][d['core'] & (d['grad'] >= 4)]).mean():.1f}%; non-core, grad>=4: {100*np.isnan(d['V'][~d['core'] & (d['grad'] >= 4)]).mean():.1f}%", ""]
    L += table(d, "A: all ocean cells", np.ones(len(d["gust"]), bool)) + [""]
    L += table(d, "C: storm-core cells (within 216 nm of a surface low below 1000 hPa)", d["core"]) + [""]
    L += table(d, "Non-core cells (comparison for C)", ~d["core"]) + [""]
    # P1 stability, gradient 6-10, winters with temperatures
    g610 = (d["grad"] >= 6) & (d["grad"] < 10) & np.isfinite(d["dT"])
    p1 = diff_boot(wmed(d, g610 & (d["dT"] >= 6)), wmed(d, g610 & (d["dT"] < 0)))
    # P2 core vs non-core, gradient >= 6, all winters
    g6 = d["grad"] >= 6
    p2 = diff_boot(wmed(d, g6 & ~d["core"]), wmed(d, g6 & d["core"]))              # non-core minus core: expected > 0
    # P3 curvature: Spearman(gust, V) - Spearman(gust, G) in core cells, grad >= 4, cells with a real root
    m = d["core"] & (d["grad"] >= 4) & np.isfinite(d["V"])
    gu, G, V, w = d["gust"][m], d["G"][m], d["V"][m], d["winter"][m]
    rg, rG, rV = ranks(gu), ranks(G), ranks(V)
    def corr(a, b, i): return np.corrcoef(a[i], b[i])[0, 1]
    base = corr(rg, rV, slice(None)) - corr(rg, rG, slice(None))
    u = np.unique(w); idx = [np.nonzero(w == x)[0] for x in u]; bs = []
    for _ in range(B):
        i = np.concatenate([idx[j] for j in rng.integers(0, len(u), len(u))]); bs.append(corr(rg, rV, i) - corr(rg, rG, i))
    ci3 = np.percentile(bs, [2.5, 97.5]); p3 = float(2 * min((np.array(bs) <= 0).mean(), (np.array(bs) >= 0).mean()))
    ps = np.array([p1[3], p2[3], max(p3, 1 / B)]); o = np.argsort(ps); q = np.empty(3); prev = 1.0
    for rk, i in zip(range(3, 0, -1), o[::-1]): prev = min(prev, ps[i] * 3 / rk); q[i] = prev
    L += ["== Planned tests (winters are the unit; 2000 winter bootstraps; q = Benjamini-Hochberg over the three) ==",
          f"P1 strongly unstable minus stable, median gust/G, grad 6-10: {p1[0]:+.3f} (95% {p1[1][0]:+.3f} to {p1[1][1]:+.3f}), {p1[2]} winters, p {p1[3]:.4f}, q {q[0]:.4f}",
          f"P2 non-core minus core, median gust/G, grad>=6: {p2[0]:+.3f} (95% {p2[1][0]:+.3f} to {p2[1][1]:+.3f}), {p2[2]} winters, p {p2[3]:.4f}, q {q[1]:.4f}",
          f"P3 Spearman(gust, gradient wind) minus Spearman(gust, geostrophic), core, grad>=4, real root: {base:+.4f} (95% {ci3[0]:+.4f} to {ci3[1]:+.4f}), n {int(m.sum())}, {len(u)} winters, p {p3:.4f}, q {q[2]:.4f}"]
    # post hoc (added after P2 was read, see README): same contrasts within gradient bins, so the mix of bins cannot drive them
    def strat(m1, m2, edges, wins):
        per = {}
        for w in wins:
            ds = []
            for a, b in zip(edges[:-1], edges[1:]):
                x = m1 & (d["winter"] == w) & (d["grad"] >= a) & (d["grad"] < b); y = m2 & (d["winter"] == w) & (d["grad"] >= a) & (d["grad"] < b)
                if x.sum() >= 20 and y.sum() >= 20: ds.append(np.median(d["gust"][x] / d["G"][x]) - np.median(d["gust"][y] / d["G"][y]))
            if ds: per[w] = float(np.mean(ds))
        v = np.array(list(per.values())); ix = rng.integers(0, len(v), (B, len(v))); bsx = v[ix].mean(1)
        return v.mean(), np.percentile(bsx, [2.5, 97.5]), len(v)
    pb = np.isfinite(d["dT"])
    h2 = strat(~d["core"], d["core"], [6, 8, 10, 12, 99], np.unique(d["winter"]))
    h1 = strat(pb & (d["dT"] >= 6), pb & (d["dT"] < 0), [2, 4, 6], np.unique(d["winter"][pb]))
    h1b = strat(pb & (d["dT"] >= 6), pb & (d["dT"] < 0), [6, 8, 10], np.unique(d["winter"][pb]))
    L += ["== Post hoc (added after the pooled P2 was read; not in the plan) ==",
          f"P2b non-core minus core, median gust/G averaged over gradient bins 6-8, 8-10, 10-12, 12+: {h2[0]:+.3f} (95% {h2[1][0]:+.3f} to {h2[1][1]:+.3f}), {h2[2]} winters",
          f"P1b strongly unstable minus stable, averaged over bins 6-8, 8-10: {h1b[0]:+.3f} (95% {h1b[1][0]:+.3f} to {h1b[1][1]:+.3f}), {h1b[2]} winters",
          f"P1c strongly unstable minus stable, averaged over bins 2-4, 4-6: {h1[0]:+.3f} (95% {h1[1][0]:+.3f} to {h1[1][1]:+.3f}), {h1[2]} winters"]
    open(os.path.join(out, "gradient_wind_tables.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
