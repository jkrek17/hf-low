"""Held-out test of a GEFS + ECMWF ENS blend (see PREREGISTRATION.md). usage: blend.py JOINT.npz ../results/gefs_primary.json OUT_PREFIX"""
import sys, json, datetime as dt
import numpy as np
sys.path.insert(0, ".."); import stats as S

SEED = 20261029; NBOOT = 2000; WGRID = np.round(np.arange(0, 1.01, 0.1), 1)
BANDS = [0, 0.05, 0.20, 0.50, 0.80, 1.0001]


def logit_cal(n, y, N):
    """-> p[k] for k = 0..N: logistic on logit(clip(k/N)) fitted over k >= 1; k = 0 keeps its own training rate (half-case prior)."""
    a, b = S.slope(n, y, N); k = np.arange(N + 1); f = np.clip(k / N, 0.5 / N, 1 - 0.5 / N); x = np.log(f / (1 - f))
    if not np.isfinite(b): p = np.full(N + 1, y[1:].sum() / max(n[1:].sum(), 1))
    else: p = 1 / (1 + np.exp(-(a + b * x)))
    p[0] = (y[0] + 0.5 * y.sum() / max(n.sum(), 1)) / (n[0] + 0.5); return np.clip(p, 0, 1)


def brier2(n, y, p):
    return float((y * (1 - p) ** 2 + (n - y) * p ** 2).sum())


def run_cell(dates, n, y, pG, N_ens=51):
    """dates: run dates; n, y: (dates, 32, 52). Returns per-week Brier sums for methods, pooled held-out calibration bins, fitted weights."""
    wk = np.array([dt.datetime.strptime(d, "%Y%m%d").isocalendar()[1] for d in dates]); weeks = sorted(set(wk))
    res = {m: {} for m in ("gefs", "blend", "blend_half", "gefs_recal", "blend_vs_recal_base")}; clim = {}; ws = {}
    bins = np.zeros((2, len(BANDS) - 1))
    base = y.sum() / n.sum()
    for h in weeks:
        tr, te = wk != h, wk == h
        ntr, ytr = n[tr].sum(0), y[tr].sum(0); nte, yte = n[te].sum(0), y[te].sum(0)
        pE = logit_cal(ntr.sum(0), ytr.sum(0), 51)[None, :]             # ENS calibration over ke
        pGr = logit_cal(ntr.sum(1), ytr.sum(1), 31)[:, None]            # GEFS recalibrated on training weeks only
        pg = pG[:, None]
        bs = [brier2(ntr, ytr, w * pg + (1 - w) * pE) for w in WGRID]; w = float(WGRID[int(np.argmin(bs))]); ws[int(h)] = w
        pb = w * pg + (1 - w) * pE
        res["gefs"][h] = brier2(nte, yte, pg + 0 * pE); res["blend"][h] = brier2(nte, yte, pb); res["blend_half"][h] = brier2(nte, yte, 0.5 * pg + 0.5 * pE)
        res["gefs_recal"][h] = brier2(nte, yte, pGr + 0 * pE)
        bw = [brier2(ntr, ytr, w2 * pGr + (1 - w2) * pE) for w2 in WGRID]; w2 = float(WGRID[int(np.argmin(bw))])
        res["blend_vs_recal_base"][h] = brier2(nte, yte, w2 * pGr + (1 - w2) * pE)       # blend built on the recalibrated GEFS
        clim[h] = brier2(nte, yte, np.full_like(pb, base))
        for j in range(len(BANDS) - 1):
            m = (pb >= BANDS[j]) & (pb < BANDS[j + 1]) & np.ones_like(nte, bool)
            bins[0, j] += nte[m].sum(); bins[1, j] += yte[m].sum()
    return weeks, res, clim, ws, bins, int(y.sum()), int(n.sum())


if __name__ == "__main__":
    z = np.load(sys.argv[1]); K, N, Y = z["keys"], z["n"], z["y"]
    tab = json.load(open(sys.argv[2]))["primary"]; out = sys.argv[3]; rng = np.random.default_rng(SEED)
    cells = {}; skipped = []
    for basin in ("natl", "npac"):
        for T, truth in ((48, "ascat"), (64, "opc")):
            for lead in (24, 48, 72, 96):
                sel = [i for i, k in enumerate(K) if k[2] == basin and int(k[3]) == T and k[4] == truth and int(k[1]) == lead]
                dates = [K[i][0] for i in sel]; name = f"{basin}_{T}_{truth}_{lead}"
                if int(Y[sel].sum()) < S.MIN_EVENTS: skipped.append((name, int(Y[sel].sum()))); continue
                pG = np.array(tab[name]["cal_table"])
                cells[name] = run_cell(dates, N[sel], Y[sel], pG)
    names = list(cells); allw = sorted(set(w for c in cells.values() for w in c[0]))
    def bss_diff(wsel, A, B):
        """mean over cells of BSS(A) - BSS(B) on the resampled weeks; wsel has repeats."""
        d = []
        for nm in names:
            weeks, res, clim = cells[nm][:3]; cnt = {w: int((np.array(wsel) == w).sum()) for w in weeks}
            den = sum(cnt[w] * clim[w] for w in weeks)
            if den <= 0: continue
            sa = sum(cnt[w] * res[A][w] for w in weeks); sb = sum(cnt[w] * res[B][w] for w in weeks)
            d.append((sb - sa) / den)
        return np.mean(d)
    def cell_diff(wsel, nm, A, B):
        weeks, res, clim = cells[nm][:3]; cnt = {w: int((np.array(wsel) == w).sum()) for w in weeks}; den = sum(cnt[w] * clim[w] for w in weeks)
        return (sum(cnt[w] * res[B][w] for w in weeks) - sum(cnt[w] * res[A][w] for w in weeks)) / den if den > 0 else np.nan
    full = list(allw); comps = [("blend", "gefs"), ("blend_half", "gefs"), ("blend_vs_recal_base", "gefs_recal")]
    lines = [f"cells scored {len(names)} of 16; skipped (fewer than 30 events): {skipped}"]; summary = {}
    boot = {c: [] for c in comps}; cellboot = {nm: [] for nm in names}
    for _ in range(NBOOT):
        pick = list(rng.choice(allw, len(allw)))
        for c in comps: boot[c].append(bss_diff(pick, *c))
        for nm in names: cellboot[nm].append(cell_diff(pick, nm, "blend", "gefs"))
    for c in comps:
        d0 = bss_diff(full, *c); lo, hi = np.percentile(boot[c], [2.5, 97.5]); summary[f"{c[0]}_minus_{c[1]}"] = dict(D=d0, lo=lo, hi=hi)
        lines.append(f"PRIMARY" if c == comps[0] else "SECONDARY"); lines[-1] += f" mean BSS {c[0]} minus {c[1]}: {d0:+.4f} (95% {lo:+.4f} to {hi:+.4f})" + ("   RULE: blend beats GEFS alone" if (c == comps[0] and lo > 0) else "")
    pv = []; rows = []
    for nm in names:
        weeks, res, clim, ws, bins, ev, cs = cells[nm]; d0 = cell_diff(full, nm, "blend", "gefs"); b = np.array(cellboot[nm]); p = (np.sum(b <= 0) + 1) / (len(b) + 1); pv.append(p)
        den = sum(clim.values()); bss = {m: 1 - sum(res[m].values()) / den for m in res}
        rows.append((nm, ev, cs, bss, d0, np.percentile(b, [2.5, 97.5]), p, ws))
    rej = S.bh(pv)
    for (nm, ev, cs, bss, d0, ci, p, ws), r in zip(rows, rej):
        lines.append(f"{nm}: cases {cs:,} events {ev:,} | BSS gefs {bss['gefs']:+.3f} ens-blend {bss['blend']:+.3f} half {bss['blend_half']:+.3f} gefs-recal {bss['gefs_recal']:+.3f} blend-on-recal {bss['blend_vs_recal_base']:+.3f} | blend minus gefs {d0:+.3f} ({ci[0]:+.3f} to {ci[1]:+.3f}) p={p:.3f} BH={r} | fitted w by held-out week: median {np.median(list(ws.values())):.1f}, range {min(ws.values()):.1f}-{max(ws.values()):.1f}")
    for T in (48, 64):
        bn = sum(cells[nm][4] for nm in names if f"_{T}_" in nm)
        if isinstance(bn, np.ndarray):
            lines.append(f"blend reliability {T} kt (held-out chance band: cases, observed): " + "; ".join(f"{BANDS[j]*100:.0f}-{min(BANDS[j+1],1)*100:.0f}%: {int(bn[0,j]):,}, {bn[1,j]/bn[0,j]*100:.1f}%" if bn[0, j] else f"{BANDS[j]*100:.0f}%: none" for j in range(len(BANDS) - 1)))
    open(out + ".txt", "w").write("\n".join(lines) + "\n"); json.dump(summary, open(out + ".json", "w"), default=float, indent=1); print("\n".join(lines))
