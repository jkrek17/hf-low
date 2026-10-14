"""Step 3: calibration of gust band -> OPC-agreed share, leave-one-season-out (PREREGISTRATION.md). GFS gust vs OPC archive; 5 seasons, Oct-Apr.
usage: analyse.py results/point_days.csv"""
import sys, os, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
BL = ["<40", "40-48", "48-56", "56-64", "64-72", "72-80", ">=80"]
d = pd.read_csv(sys.argv[1], dtype={"cycle": str})
d["season"] = [int(c[:4]) if int(c[4:6]) >= 10 else int(c[:4]) - 1 for c in d.cycle]
seasons = sorted(d.season.unique()); print("seasons", seasons, "cycles", d.cycle.nunique())
d["basin2"] = d.basin
pool = d.copy(); pool["basin"] = "both"; D = pd.concat([d, pool], ignore_index=True)
G = D.groupby(["basin", "day", "band", "season"], as_index=False).agg(n=("n", "sum"), hf=("hf", "sum"), storm=("storm", "sum"))
ev = D.assign(hfc=(D.hf > 0).astype(int), stc=(D.storm > 0).astype(int)).groupby(["basin", "day", "band", "season"], as_index=False).agg(cyc_hf=("hfc", "sum"), cyc_st=("stc", "sum"))
G = G.merge(ev, on=["basin", "day", "band", "season"])
lines = ["GFS daily-max gust band vs OPC agreed (HF = Category HF fix, storm = S or HF) within 324 nm during the forecast day; 1-degree ocean points, Oct-Apr, "
         f"{len(seasons)} seasons {seasons[0]}-{str(seasons[0] + 1)[-2:]} to {seasons[-1]}-{str(seasons[-1] + 1)[-2:]}. GFS gust is a PROXY for sustained wind; OPC's call is an analyst judgement.", ""]
rows = []
for (bas, day), g in G.groupby(["basin", "day"]):
    for b in range(7):
        x = g[g.band == b]
        n, hf, st = x.n.sum(), x.hf.sum(), x.storm.sum()
        ph = hf / n if n else np.nan; ps = st / n if n else np.nan
        per_hf = [(r.hf / r.n if r.n else np.nan) for r in x.itertuples()]
        # leave-one-season-out training tables for the "usable" rule
        tr_n, tr_c = [], []
        for s in seasons:
            t = x[x.season != s]; tr_n.append(t.n.sum()); tr_c.append(t.cyc_hf.sum())
        rows.append(dict(basin=bas, day=day, band=b, band_kt=BL[b], n=n, hf_share=ph, storm_share=ps, cycles_hf=int(x.cyc_hf.sum()), cycles_storm=int(x.cyc_st.sum()),
                         season_min_hf=np.nanmin(per_hf) if per_hf else np.nan, season_max_hf=np.nanmax(per_hf) if per_hf else np.nan,
                         usable=bool(min(tr_n) >= 100 and min(tr_c) >= 5)))
T = pd.DataFrame(rows); T.to_csv(os.path.join(RES, "bands_pooled.csv"), index=False)
# LOSO Brier
out = []
for (bas, day), g in G.groupby(["basin", "day"]):
    for tgt in ("hf", "storm"):
        tot_b = tot_c = tot_n = 0; per = []
        for s in seasons:
            tr = g[g.season != s].groupby("band")[["n", tgt]].sum(); te = g[g.season == s].set_index("band")
            p = (tr[tgt] / tr["n"]).reindex(range(7)).fillna(0.0)
            n = te["n"]; y = te[tgt]; pp = p.reindex(te.index)
            br = float(((y * (1 - pp) ** 2 + (n - y) * pp ** 2).sum()) / n.sum())
            base = y.sum() / n.sum(); brc = float(base * (1 - base))
            btr = tr[tgt].sum() / tr["n"].sum(); brt = float((y * (1 - btr) ** 2 + (n - y) * btr ** 2).sum() / n.sum())
            per.append((s, br, brc, brt)); tot_b += br * n.sum(); tot_c += brc * n.sum(); tot_n += n.sum()
        bss = [1 - a / b for _, a, b, _ in per]; bss_tr = [1 - a / c for _, a, _, c in per]
        out.append(dict(basin=bas, day=day, target=tgt, brier=tot_b / tot_n, bss_vs_heldout_base=1 - tot_b / tot_c, bss_min=min(bss), bss_max=max(bss),
                        bss_vs_training_base_min=min(bss_tr), bss_vs_training_base_max=max(bss_tr)))
L = pd.DataFrame(out); L.to_csv(os.path.join(RES, "loso_brier.csv"), index=False)
for bas in ("atl", "pac", "both"):
    lines.append(f"== {bas}: OPC HF-agreed share (storm-agreed) by day and gust band; n point-days; cycles with an HF event; per-season range of HF share ==")
    for day in (1, 2, 3, 4):
        s = T[(T.basin == bas) & (T.day == day)].sort_values("band")
        lines.append(f"day {day}: " + " | ".join(f"{r.band_kt}: {100 * r.hf_share:.1f}% ({100 * r.storm_share:.1f}%) n={r.n} c={r.cycles_hf}{'' if r.usable else ' (too few)'}" for r in s.itertuples()))
    lines.append("")
lines.append("== Leave-one-season-out Brier skill of the banded table (vs the held-out season's own base rate); range over held-out seasons ==")
for r in L.sort_values(["target", "basin", "day"]).itertuples():
    lines.append(f"{r.target:5s} {r.basin:4s} day {r.day}: BSS {r.bss_vs_heldout_base:.3f} (seasons {r.bss_min:.3f} to {r.bss_max:.3f}; vs training base {r.bss_vs_training_base_min:.3f} to {r.bss_vs_training_base_max:.3f})")
open(os.path.join(RES, "calibration.txt"), "w").write("\n".join(lines)); print("\n".join(lines))
