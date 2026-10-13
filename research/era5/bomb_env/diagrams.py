"""Observed rates in the two upper-level diagrams (6 x 6 quantile bins, all 2004-05..2025-26 fixes pooled).
Display only: skill is scored by score.py (models A and B). ERA5 proxy, pipeline A.
usage: diagrams.py [OUTDIR]
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score as S

PAIRS = {"jet_trough": ("jet250", "vadv500"), "jet_baroclinic": ("jet250", "eady")}
NB, MINN = 6, 30


def wilson(k, n, z=1.645):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results")
    D = S.load()
    rows = []
    for name, (xa, ya) in PAIRS.items():
        ex = np.nanquantile(D[xa], np.linspace(0, 1, NB + 1))
        ey = np.nanquantile(D[ya], np.linspace(0, 1, NB + 1))
        D["bx"] = np.clip(np.digitize(D[xa], ex[1:-1]), 0, NB - 1)
        D["by"] = np.clip(np.digitize(D[ya], ey[1:-1]), 0, NB - 1)
        for (bx, by), g in D.groupby(["bx", "by"]):
            r = g[g.cls >= 0]
            n, nr = len(g), len(r)
            k = int((r.cls == 4).sum())
            lo, hi = wilson(k, nr) if nr else (np.nan, np.nan)
            sh = np.bincount(g.gcls, minlength=4) / n
            rows.append(dict(diagram=name, xvar=xa, yvar=ya, bx=bx, by=by,
                             x_lo=ex[bx], x_hi=ex[bx + 1], y_lo=ey[by], y_hi=ey[by + 1],
                             n=n, n_scored=nr, rd_rate=k / nr if nr else np.nan, rd_lo90=lo, rd_hi90=hi,
                             p_below_gale=sh[0], p_gale=sh[1], p_storm=sh[2], p_hf=sh[3], shown=n >= MINN))
    C = pd.DataFrame(rows).round(5)
    C.to_csv(os.path.join(out, "diagram_cells.csv"), index=False)
    base = (D.cls[D.cls >= 0] == 4).mean()
    lines = [f"Rapid-deepening base rate {base:.4f}"]
    for name in PAIRS:
        c = C[C.diagram == name]
        lines.append(f"{name}: rd rate min {c.rd_rate.min():.3f}, max {c.rd_rate.max():.3f}; cells shown {int(c['shown'].sum())}/36; min n {c.n.min()}")
        top = c.sort_values("rd_rate").iloc[-1]
        lines.append(f"  highest cell: x {top.x_lo:.1f}..{top.x_hi:.1f} {top.xvar}, y {top.y_lo:.2f}..{top.y_hi:.2f} {top.yvar}, rate {top.rd_rate:.3f} (n {int(top.n_scored)})")
        bot = c.sort_values("rd_rate").iloc[0]
        lines.append(f"  lowest cell: rate {bot.rd_rate:.3f} (n {int(bot.n_scored)})")
    open(os.path.join(out, "diagrams.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
