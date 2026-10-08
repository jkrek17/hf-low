"""Readable tables and one figure from results/*.csv.   python3 make_report.py results"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1]
m = pd.read_csv(os.path.join(R, "monthly_tables.csv"))
t = pd.read_csv(os.path.join(R, "tests.csv"))
ORDER = [6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5]
NAME = "JFMAMJJASOND"
mn = {1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "Jun", 7: "Jul", 8: "Aug", 9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec"}
out = []


def col(q, b, field="mean_per_season"):
    d = m[(m.quantity == q) & (m.basin == b)].set_index("month")[field]
    return d.reindex(ORDER)


for b, nm in (("pac", "Pacific"), ("atl", "Atlantic")):
    df = pd.DataFrame({
        "archive HF/season": col("Q1", b),
        "A HF/season": col("Q2", b),
        "A cyclones/season": col("Q3", b),
        "A HF share %": col("Q4", b, "value") * 100,
        "A real cyc/season": col("Q5", b),
        "depth-strong/season (47 s)": col("Q7", b),
        "v250 eddy (m2/s2)": col("Q9", b),
        "MSLP eddy (Pa2)": col("Q10", b),
    })
    df.index = [mn[i] for i in df.index]
    out.append(f"{nm}: monthly means over seasons (archive and A: 2004-05..2025-26; depth 1979-80..2025-26; v250 2004-05..2020-21; MSLP 1979-80..2020-21)")
    out.append(df.round(2).to_string())
    out.append("")
lat = []
for b in ("pac", "atl"):
    lat.append(pd.DataFrame({f"{b} archive first-HF lat": col("Q11", b, "value"), f"{b} A onset lat": col("Q12", b, "value"), f"{b} A peak lat (all tracks)": col("Q13", b, "value")}).round(1))
L = pd.concat(lat, axis=1)
L.index = [mn[i] for i in L.index]
out.append("Mean latitude by month (deg N)")
out.append(L.to_string())
out.append("")
c = t.copy()
c["x"] = np.where(c.log, (np.exp(c.estimate) - 1) * 100, np.nan)
c["x_lo"] = np.where(c.log, (np.exp(c.lo) - 1) * 100, np.nan)
c["x_hi"] = np.where(c.log, (np.exp(c.hi) - 1) * 100, np.nan)
out.append("All 56 tests (log quantities: estimate is a log ratio, % is exp(estimate)-1; latitudes in degrees). q_bh is across all 56.")
out.append(c[["test", "estimate", "lo", "hi", "x", "x_lo", "x_hi", "p", "q_bh", "p_bonf7", "loo_min", "loo_max", "detectable_80", "n_seasons"]].round(3).to_string(index=False))
open(os.path.join(R, "summary.txt"), "w").write("\n".join(out) + "\n")

fig, axs = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True)
xs = np.arange(12)
labs = [mn[i] for i in ORDER]
for k, (b, nm) in enumerate((("pac", "Pacific"), ("atl", "Atlantic"))):
    a = axs[0, k]
    for q, lab, st, c_ in (("Q1", "archive HF count", "-o", "C0"), ("Q2", "pipeline A HF count (proxy)", "--s", "C1"), ("Q3", "all cyclones (A)", ":", "C2")):
        v = col(q, b, "mean_per_season").values
        v = v / np.nanmean(v[[4, 5, 6, 7, 8, 9, 10]])
        a.plot(xs, v, st, color=c_, label=lab, ms=4)
    v = col("Q9", b).values
    a.plot(xs, v / np.nanmean(v[[4, 5, 6, 7, 8, 9, 10]]), "-^", color="k", ms=4, label="v250 2-6 d eddy activity")
    a.set_title(f"{nm}: relative to Oct-Apr mean")
    a.set_ylabel("monthly value / Oct-Apr mean")
    if k == 0:
        a.legend(fontsize=7)
    a = axs[1, k]
    a.plot(xs, col("Q4", b, "value").values * 100, "-o", color="C3", ms=4, label="HF share of cyclones, % (A)")
    a.set_ylabel("HF share, %", color="C3")
    a2 = a.twinx()
    keep = np.array([i in (4, 5, 6, 7, 8, 9, 10) for i in range(12)])      # Oct-Apr: other months have too few events
    a2.plot(xs[keep], col("Q12", b, "value").values[keep], "-s", color="C4", ms=4, label="A mean HF onset latitude")
    a2.plot(xs[keep], col("Q11", b, "value").values[keep], "-o", color="C5", ms=4, label="archive mean first-HF latitude")
    a2.set_ylabel("latitude, deg N")
    a2.legend(fontsize=7, loc="lower right")
    a.set_xticks(xs)
    a.set_xticklabels(labs)
for a in axs.ravel():
    a.grid(alpha=0.3)
fig.suptitle("Seasonal cycle of HF lows: archive and ERA5 proxy (pipeline A), with eddy activity", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(R, "seasonal_cycle.png"), dpi=110)
