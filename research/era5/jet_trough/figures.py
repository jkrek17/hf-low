"""Held-out curves: observed frequency by fit-season decile (points, 90% season-bootstrap bars) and the
adjusted curves (lines: linear, hinge, spline; g-computation on held-out covariates). ERA5 proxy."""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OD = sys.argv[1] if len(sys.argv) > 1 else "results"
T = pd.read_csv(f"{OD}/thresholds.csv"); C = pd.read_csv(f"{OD}/curves.csv"); R = pd.read_csv(f"{OD}/raw_deciles.csv")
NAME = {"jet250": "250 hPa jet maximum within 1000 km (kt)", "trough_up": "Upstream 500 hPa trough depth (m below the monthly mean)"}
for var, lab in NAME.items():
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), sharex="col")
    for i, basin in enumerate(("atl", "pac")):
        for j, out in enumerate(("BOMB", "HFON")):
            a = ax[i, j]
            q = lambda d: d[(d.basin == basin) & (d.outcome == out) & (d["var"] == var) & (d.tag == "main")]
            r, c, t = q(R), q(C), q(T).iloc[0]
            xm = c.x.values
            mid = [np.nan] * 10
            # decile midpoints in x: use the x_lo of the next decile, mean of neighbours
            lo = r.x_lo.values.copy(); lo[0] = c.x.min() - (c.x.iloc[1] - c.x.iloc[0])
            ed = np.append(lo, c.x.max() + (c.x.iloc[-1] - c.x.iloc[-2]))
            mid = 0.5 * (ed[:-1] + ed[1:])
            a.errorbar(mid, r.freq_test, yerr=[r.freq_test - r.lo_test, r.hi_test - r.freq_test], fmt="o", ms=4, color="0.35",
                       capsize=2, label="observed, held-out deciles")
            a.plot(xm, c.lin, color="#1b9e77", lw=1.3, label="straight line")
            a.plot(xm, c.spl, color="#7570b3", lw=1.3, label="smooth spline")
            a.plot(xm, c.hinge, color="#d95f02", lw=1.6, label="one knot")
            a.axvline(t.knot, color="#d95f02", ls=":", lw=1)
            a.set_title(f"{'Atlantic' if basin == 'atl' else 'Pacific'}, {'bomb' if out == 'BOMB' else 'HF onset'} within 24 h", fontsize=10)
            a.set_ylim(bottom=0)
            if j == 0:
                a.set_ylabel("probability")
            if i == 1:
                a.set_xlabel(lab, fontsize=9)
            if i == 0 and j == 0:
                a.legend(fontsize=7, frameon=False)
    fig.suptitle(f"{var}: held-out seasons 2015-16..2025-26 (ERA5 proxy, pipeline A)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OD}/curves_{var}.png", dpi=110)
    plt.close(fig)
