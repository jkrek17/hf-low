"""Step 3: season-block bootstrap of coverage / shading cost / precision by radius and lead, decision rule (PREREGISTRATION.md).
ERA5 PROXY gust, position-error-only upper bound. Reads results/season_sums_*.csv."""
import os, glob, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
d = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(RES, "season_sums_*.csv"))])
d = d.groupby(["season", "basin", "cls", "lead", "R_nm"], as_index=False).sum(numeric_only=True)
seasons = np.sort(d.season.unique()); print(len(seasons), "seasons", seasons.min(), seasons.max())
rng = np.random.default_rng(20261014); B = 2000
idx = rng.integers(0, len(seasons), (B, len(seasons)))
out = []
for (basin, cls, lead, R), g in d.groupby(["basin", "cls", "lead", "R_nm"]):
    g = g.set_index("season").reindex(seasons).fillna(0)
    def stats(w):
        T, Tc, sh, F, Fn, dr = (np.tensordot(w, g[c].to_numpy(), axes=(-1, 0)) for c in ["T", "T_covered", "shaded", "F", "F_near_T", "draws"])
        return Tc / T, sh / np.maximum(F, 1), Fn / np.maximum(F, 1)   # coverage, shaded cells per forecast cell, precision
    cnt = np.zeros((B, len(seasons)))
    for b in range(B): cnt[b] = np.bincount(idx[b], minlength=len(seasons))
    pt = stats(np.ones(len(seasons))); bs = stats(cnt)
    row = dict(basin=basin, cls=cls, lead=lead, R_nm=R, seasons=len(seasons))
    for name, p, b in zip(["coverage", "shade_ratio", "precision"], pt, bs):
        row[name] = p; row[name + "_lo"], row[name + "_hi"] = np.percentile(b, [2.5, 97.5])
    out.append(row)
r = pd.DataFrame(out); r.to_csv(os.path.join(RES, "radius_results.csv"), index=False)
lines = ["Radius test (ERA5 PROXY gust; displacement-only upper bound; 22 Oct-Apr seasons; nautical miles)", ""]
for cls in ["hf", "storm", "gale"]:
    lines.append(f"== {cls.upper()} coverage (share of truth cells inside shading), 95% season-bootstrap interval ==")
    for basin in ["atl", "pac"]:
        for lead in sorted(r.lead.unique()):
            s = r[(r.basin == basin) & (r.cls == cls) & (r.lead == lead)].sort_values("R_nm")
            lines.append(f"{basin} {lead:>2} h: " + "  ".join(f"{x.R_nm}:{x.coverage:.2f} ({x.coverage_lo:.2f}-{x.coverage_hi:.2f})" for x in s.itertuples()))
    lines.append("")
lines.append("== Shading cost (shaded cells per forecast-wind cell), HF ==")
for basin in ["atl", "pac"]:
    for lead in sorted(r.lead.unique()):
        s = r[(r.basin == basin) & (r.cls == "hf") & (r.lead == lead)].sort_values("R_nm")
        lines.append(f"{basin} {lead:>2} h: " + "  ".join(f"{x.R_nm}:{x.shade_ratio:.1f}x" for x in s.itertuples()))
lines += ["", "== Decision rule: smallest R with HF coverage >= 0.80 in BOTH basins (point estimate; lower bound in brackets) =="]
for lead in sorted(r.lead.unique()):
    h = r[(r.cls == "hf") & (r.lead == lead)]
    ok = [R for R in sorted(h.R_nm.unique()) if (h[h.R_nm == R].coverage >= 0.8).all()]
    okl = [R for R in sorted(h.R_nm.unique()) if (h[h.R_nm == R].coverage_lo >= 0.8).all()]
    lines.append(f"{lead:>2} h: point {ok[0] if ok else 'none within 432'}; lower bound {okl[0] if okl else 'none within 432'}")
open(os.path.join(RES, "radius_results.txt"), "w").write("\n".join(lines)); print("\n".join(lines))
