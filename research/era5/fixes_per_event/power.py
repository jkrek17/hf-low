"""RA-11 power: planted-effect simulation (plan: PREREGISTRATION.md, 'Power plan').

Real archive events, linear trend removed (event-level OLS), whole seasons resampled with replacement into the 22 slots
(event dates = target-season June 1 + the event's day offset), then a planted linear decline (0.181 /decade) or a step
(0.40 at each sensor date) is added. Tests as in analyse.py: P1 t-test on season means; steps by season-block bootstrap
(500 draws, shared); S2 AIC rule. 1000 runs per truth scenario; MDE grid 400 runs per cell.
usage: python3 power.py
"""
import math, os
import numpy as np, pandas as pd
import analyse as A

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(7)
ev = A.load_archive()
S = A.SEASONS; NS = len(S)
x = (ev.season.values - ev.season.mean())
b = (x * (ev.n_hf - ev.n_hf.mean())).sum() / (x ** 2).sum()
ev["res"] = ev.n_hf - b * x - 0 * ev.n_hf.mean() + 0  # detrended, keeps the mean via -b*x only
ev["off"] = (ev.t0 - pd.to_datetime(ev.season.astype(str) + "-06-01")).dt.total_seconds() / 86400.0
src = {s: (ev.res.values[ev.season == s], ev.off.values[ev.season == s]) for s in S}
DATES = [d for _, d in A.STEPS]
NAMES = [n.split()[0] for n, _ in A.STEPS]
W = A.boot_weights(np.random.default_rng(5), 500)
start = {s: pd.Timestamp("%d-06-01" % s) for s in S}
cut = [(pd.Timestamp(d) - pd.Timestamp("1970-01-01")).total_seconds() / 86400 for d in DATES]
base_day = np.array([(start[s] - pd.Timestamp("1970-01-01")).days for s in S], float)


def one(lin=0.0, step=0.0, k=None):
    vals, sidx, days = [], [], []
    for j in range(NS):
        r, o = src[S[rng.integers(NS)]]
        vals.append(r); days.append(o + base_day[j]); sidx.append(np.full(len(r), j))
    v = np.concatenate(vals); si = np.concatenate(sidx); dy = np.concatenate(days)
    v = v + 0.0 * dy - lin / 10.0 * (S[si] - S.mean())
    if k is not None:
        v = v - step * (dy >= cut[k])
    out = {}
    sm = np.bincount(si, v, NS) / np.bincount(si, minlength=NS)
    b_, se_, t_, p_, _, _ = A.ols(S, sm)
    out["P1"] = p_ < 0.05
    pick = {}
    for q in range(4):
        aft = dy >= cut[q]
        Sa = np.bincount(si[aft], v[aft], NS); Na = np.bincount(si[aft], minlength=NS).astype(float)
        Sb = np.bincount(si[~aft], v[~aft], NS); Nb = np.bincount(si[~aft], minlength=NS).astype(float)
        st = dict(Sa=Sa, Na=Na, Sb=Sb, Nb=Nb)
        est = float(A.step_diff(st)[0])
        pb = A.boot_p(est, A.step_diff(st, W))[0]
        out[NAMES[q]] = pb < 0.05
        fr = np.bincount(si[aft], minlength=NS) / np.bincount(si, minlength=NS)
        a0, a1, a2 = A.form_table(sm, S, fr)
        pick[NAMES[q]] = A.prefer(a1, a2)
    return out, pick


def scen(label, n=1000, **kw):
    R, Pk = [], []
    for _ in range(n):
        o, p = one(**kw); R.append(o); Pk.append(p)
    df = pd.DataFrame(R)
    line = "%-34s P1 %.2f | steps " % (label, df.P1.mean()) + " ".join("%s %.2f" % (c, df[c].mean()) for c in NAMES)
    return line, df, Pk


lines = ["RA-11 power (planted effects on detrended real events; alpha 0.05; 22 seasons; steps by season-block bootstrap)"]
lines.append("rejection rates (fraction of 1000 runs with p<0.05):")
l, _, _ = scen("null (nothing planted)"); lines.append(l)
l, _, pk = scen("linear -0.181 /decade", lin=0.181); lines.append(l)
conf = {}
for k, nm in enumerate(NAMES):
    l, _, pk = scen("step -0.40 at %s" % nm, step=0.40, k=k); lines.append(l)
    c = pd.Series([p[nm] for p in pk]).value_counts(normalize=True)
    conf[nm] = "S2 picks (truth=step %s): step %.2f, linear %.2f, tie %.2f" % (nm, c.get("step", 0), c.get("linear", 0), c.get("tie", 0))
l, _, pk = scen("linear -0.181 (S2 truth=linear)", lin=0.181)
for nm in NAMES:
    c = pd.Series([p[nm] for p in pk]).value_counts(normalize=True)
    conf["lin_" + nm] = "S2 vs %s step when truth is linear: step %.2f, linear %.2f, tie %.2f" % (nm, c.get("step", 0), c.get("linear", 0), c.get("tie", 0))
lines.append("")
lines += list(conf.values())
lines.append("")
lines.append("MDE: smallest planted step (fixes per event) detected in >=80% of 400 runs, per date")
grid = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80]
for k, nm in enumerate(NAMES):
    mde = None; rates = []
    for g in grid:
        r = np.mean([one(step=g, k=k)[0][nm] for _ in range(400)])
        rates.append("%.2f:%.2f" % (g, r))
        if mde is None and r >= 0.8: mde = g
    lines.append("  %s MDE %s   (size:power %s)" % (nm, mde if mde else ">0.80", " ".join(rates)))
lm = None
for g in [0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30]:
    r = np.mean([one(lin=g)[0]["P1"] for _ in range(400)])
    if lm is None and r >= 0.8: lm = g
lines.append("  linear slope MDE (per decade, P1 t-test) %s" % lm)
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
open(os.path.join(HERE, "results", "power.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
