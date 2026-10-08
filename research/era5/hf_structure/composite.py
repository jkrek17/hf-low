"""Storm-relative climatology of hurricane-force-strength gusts (ERA5 proxy, pipeline A).

Reads extract.py's per-fix statistics and boxes, attaches a life-cycle stage
from pipeline A's tracks and Hart's phase from the near-storm framework, and
writes:

  results/fixes.csv        one row per HF-strength fix (the committed table)
  results/structure.txt    the numbers, by basin, stage and Hart phase
  results/composite.npz    composite boxes (mean gust, HF frequency, MSLP)
  results/composite_*.png  figures

Definitions (the sample is every in-domain point of a pipeline A event track
with gust index g800 >= 71.7 kt, seasons 2004-05 to 2025-26):

  HF-equivalent area  owned ocean area within 1200 km where the ERA5 gust is
                      >= 71.7 kt, pipeline A's calibrated threshold applied
                      point by point. The threshold was calibrated on the
                      storm maximum, so this is an index of extent, not an
                      area of observed 64 kt wind.
  core area           area where the gust is >= 90% of the storm's own maximum;
                      independent of ERA5's absolute bias.
  RMG                 radius of the maximum gust (the gust analogue of RMW).
  quadrants           motion-relative: front-right (0-90 deg clockwise from the
                      direction of motion), rear-right, rear-left, front-left.
                      Earth-relative: NE, SE, SW, NW of the centre.
  stage               hours from the track's minimum MSLP (in-domain points):
                      deepening <= -12 h, mature -6..+6 h, filling >= +12 h.
                      Provisional until the life-cycle thread commits its
                      definitions.
  Hart phase          at 00/12 UTC only (env_2004.csv.gz): B >= 10 m
                      asymmetric, VTL > 0 warm core below 600 hPa.
  terrain flag        maximum gust within 100 km of land whose nearest land
                      point is Greenland or Iceland. Flagged, not studied.

Uncertainty: 90% intervals from resampling storms (the independent unit for
structure), 1000 draws.

usage: composite.py [ENV_2004_CSV_GZ]
"""
import os, sys, glob
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES_A = os.path.join(HERE, "..", "hf_history", "results")
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "hf_structure")
OUT = os.path.join(HERE, "results")
HF = 71.7
BOX_N, BOX_D = 121, 25.0
NM = 1.852
rng = np.random.default_rng(7)


def load():
    F = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(WORK, "stats", "*.csv")))],
                  ignore_index=True)
    P = pd.read_csv(os.path.join(RES_A, "era5_hf_catalog_tracks.csv"))
    P = P[P.track.isin(F.track.unique()) & P.basin.notna()].copy()
    P["dt"] = pd.to_datetime(P.time.astype(str), format="%Y%m%d%H")
    g = P.groupby("track")
    tmin = P.loc[g.msl.idxmin(), ["track", "dt"]].set_index("track").dt
    pk = P.loc[g.g800.idxmax(), ["track", "dt"]].set_index("track").dt
    hfp = P[P.g800 >= HF].groupby("track").dt
    first, last = hfp.min(), hfp.max()
    F["dt"] = pd.to_datetime(F.time.astype(str), format="%Y%m%d%H")
    h = lambda s: (F.dt - F.track.map(s)).dt.total_seconds() / 3600
    F["h_pmin"], F["h_peak"], F["h_first"] = h(tmin), h(pk), h(first)
    F["h_last"] = (F.track.map(last) - F.dt).dt.total_seconds() / 3600
    F["stage"] = np.select([F.h_pmin <= -12, F.h_pmin.abs() <= 6, F.h_pmin >= 12],
                           ["deepening", "mature", "filling"], "")
    F["episode"] = np.select([(F.h_first == 0) & (F.h_last == 0), F.h_first == 0, F.h_last == 0],
                             ["single", "first", "last"], "interior")
    # 12 h pressure change before the fix, from the track
    P = P.sort_values(["track", "dt"])
    P["dp12"] = P.groupby("track").msl.diff(2)
    P.loc[P.groupby("track").dt.diff(2) != pd.Timedelta(hours=12), "dp12"] = np.nan
    F = F.merge(P[["track", "time", "dp12"]], on=["track", "time"], how="left")
    grn = (F.gmax_coast_lat >= 59) & (F.gmax_coast_lon >= 285) & (F.gmax_coast_lon <= 350)
    F["terrain"] = (F.gmax_coast_km <= 100) & grn
    return F


def add_hart(F, path):
    if not path or not os.path.exists(path):
        F["phase"] = ""
        return F
    E = pd.read_csv(path, usecols=["track", "time", "B", "VTL"])
    F = F.merge(E, on=["track", "time"], how="left")
    F["phase"] = np.where(F.B.isna(), "",
                          np.where(F.B >= 10, np.where(F.VTL > 0, "asym-warm", "asym-cold"),
                                   np.where(F.VTL > 0, "sym-warm", "sym-cold")))
    return F


def boot(F, fns, n=1000):
    """90% interval of each fn(F) from resampling storms (one set of draws for all fns)."""
    tr = F.track.unique()
    idx = F.groupby("track").indices
    vals = []
    for _ in range(n):
        pick = rng.choice(tr, len(tr))
        d = F.iloc[np.concatenate([idx[t] for t in pick])]
        vals.append([fn(d) for fn in fns])
    return np.nanpercentile(np.array(vals, float), [5, 95], axis=0).T


QS = ["q_fr", "q_rr", "q_rl", "q_fl"]
ES = ["e_ne", "e_se", "e_sw", "e_nw"]


def summary(F):
    """The quoted numbers for one subset, each with a storm-resampled 90% interval."""
    out = {"fixes": len(F), "storms": F.track.nunique()}
    sh = lambda cols: (lambda d: 100 * d[cols].to_numpy().sum() / d.a_hf.sum())
    defs = {
        "RMG km, median": lambda d: d.gmax_r.median(),
        "max gust kt, median": lambda d: d.gmax.median(),
        "HF-eq area 1e3 km2, median": lambda d: d.a_hf.median() / 1e3,
        "HF-eq area 1e3 km2, mean": lambda d: d.a_hf.mean() / 1e3,
        "core (>=90% of max) area 1e3 km2, median": lambda d: d.a_p90.median() / 1e3,
        "gust>=48 kt area 1e3 km2, median": lambda d: d.a_g48.median() / 1e3,
        "HF-eq outer radius km, median": lambda d: d.hf_rmax.median(),
        "max gust in right half, %": lambda d: 100 * ((d.gmax_rel > 0) & (d.gmax_rel < 180)).mean(),
        "max gust in rear-right, %": lambda d: 100 * ((d.gmax_rel >= 90) & (d.gmax_rel < 180)).mean(),
        "max gust south of centre, %": lambda d: 100 * ((d.gmax_brg > 90) & (d.gmax_brg < 270)).mean(),
    }
    for q in QS + ES:
        defs[f"share of HF-eq area {q[2:].upper()}, %"] = sh([q])
    defs["share of HF-eq area south of centre, %"] = sh(["e_se", "e_sw"])
    defs["share of HF-eq area right of motion, %"] = sh(["q_fr", "q_rr"])
    ci = boot(F, list(defs.values()))
    for (k, fn), (lo, hi) in zip(defs.items(), ci):
        out[k] = (fn(F), lo, hi)
    return out


def fmt(out):
    lines = [f"  fixes {out['fixes']}, storms {out['storms']}"]
    for k, v in out.items():
        if isinstance(v, tuple):
            lines.append(f"  {k:<44s} {v[0]:8.1f}   [{v[1]:.1f}, {v[2]:.1f}]")
    return lines


def composites(F):
    keys = {}
    for b in ("atl", "pac"):
        for st in ("all", "deepening", "mature", "filling"):
            keys[(b, st)] = F[(F.basin == b) & ((F.stage == st) | (st == "all"))]
    acc = {}
    for t, d in F.groupby("time"):
        z = np.load(os.path.join(WORK, "boxes", f"{t}.npz"))
        for _, r in d.iterrows():
            tr = int(r.track)
            for (b, st), sub in keys.items():
                if r.basin != b or (st != "all" and r.stage != st):
                    continue
                a = acc.setdefault((b, st), {k: np.zeros((BOX_N, BOX_N)) for k in
                                             ("gm", "gm_n", "hm", "hm_n", "om", "om_n", "pm", "n")})
                gm, om = z[f"{tr}_gust_m"].astype(float), z[f"{tr}_own_m"].astype(float)
                gn, on = z[f"{tr}_gust_n"].astype(float), z[f"{tr}_own_n"].astype(float)
                a["gm"] += gm * om; a["om"] += om; a["hm"] += (gm >= HF) * om
                a["gm_n"] += gn * on; a["om_n"] += on; a["hm_n"] += (gn >= HF) * on
                a["pm"] += z[f"{tr}_msl_m"].astype(float); a["n"] += 1
    res = {}
    for k, a in acc.items():
        n = a["n"].max()
        res[k] = dict(gust_m=a["gm"] / np.maximum(a["om"], 1), hf_m=a["hm"] / n,
                      gust_n=a["gm_n"] / np.maximum(a["om_n"], 1), hf_n=a["hm_n"] / n,
                      msl_m=a["pm"] / n, n=n)
    return res


def figure(res, frame, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax_ = (np.arange(BOX_N) - BOX_N // 2) * BOX_D
    stages = ("all", "deepening", "mature", "filling")
    fig, axs = plt.subplots(2, 4, figsize=(15, 7.6), sharex=True, sharey=True)
    for i, b in enumerate(("atl", "pac")):
        for j, st in enumerate(stages):
            ax = axs[i, j]
            r = res[(b, st)]
            hf = r["hf_m" if frame == "m" else "hf_n"] * 100
            im = ax.pcolormesh(ax_, ax_, hf, cmap="magma_r", vmin=0, vmax=40, shading="auto")
            cs = ax.contour(ax_, ax_, r["gust_m" if frame == "m" else "gust_n"], [40, 50, 60],
                            colors="steelblue", linewidths=0.8)
            ax.clabel(cs, fontsize=7, fmt="%d")
            if frame == "m":
                ax.contour(ax_, ax_, r["msl_m"], np.arange(940, 1020, 4), colors="0.5", linewidths=0.5)
            ax.plot(0, 0, "k+", ms=8)
            ax.set_xlim(-1200, 1200); ax.set_ylim(-1200, 1200); ax.set_aspect("equal")
            ax.set_title(f"{'Atlantic' if b == 'atl' else 'Pacific'} {st}, n={int(r['n'])}", fontsize=9)
            if frame == "m":
                ax.annotate("", (300, -1100), (-300, -1100), arrowprops=dict(arrowstyle="->"))
    lab = ("km along motion (motion to the right)", "km left of motion") if frame == "m" \
        else ("km east", "km north")
    for ax in axs[1]: ax.set_xlabel(lab[0], fontsize=8)
    for ax in axs[:, 0]: ax.set_ylabel(lab[1], fontsize=8)
    fig.colorbar(im, ax=axs, shrink=0.6, label="% of fixes with ERA5 gust >= 71.7 kt (owned, ocean)")
    fig.suptitle("ERA5 proxy (pipeline A), HF-strength fixes 2004-05 to 2025-26: frequency of "
                 "HF-equivalent gust; blue = mean gust (kt)" + ("; grey = mean MSLP" if frame == "m" else ""),
                 fontsize=10)
    fig.savefig(path, dpi=80)


def main():
    os.makedirs(OUT, exist_ok=True)
    F = add_hart(load(), sys.argv[1] if len(sys.argv) > 1 else "")
    L = ["Storm-relative structure of HF-strength gusts. ERA5 PROXY, pipeline A (research/era5/hf_history).",
         "Sample: every in-domain point of a pipeline A event track with g800 >= 71.7 kt, seasons 2004-05..2025-26.",
         "ERA5 at 31 km under-resolves peak winds and narrow jets (e.g. sting jets): absolute areas at a",
         "fixed wind speed are biased low, and the HF-equivalent area uses the calibrated 71.7 kt GUST",
         "threshold, not an observed 64 kt sustained wind. Read shapes and relative numbers first.",
         "Intervals: 90%, resampling storms (1000 draws). Area shares are area-weighted over all fixes in the subset.", ""]
    L.append(f"Reproduction of pipeline A's index: g800 recomputed = catalog value at "
             f"{(F.g800 == F.g800_cat).mean() * 100:.2f}% of {len(F)} fixes; max |diff| "
             f"{(F.g800 - F.g800_cat).abs().max():.1f} kt; centre match > 0 km at {(F.match_km > 0).sum()} fixes.")
    L.append(f"Terrain flag (max gust within 100 km of Greenland/Iceland): {int(F.terrain.sum())} fixes "
             f"({F.terrain.mean() * 100:.1f}%), Atlantic {int(F[F.basin == 'atl'].terrain.sum())}. Flagged, not studied;"
             f" results below include them, and a line per basin repeats the headline without them.")
    L.append("")
    for b, name in (("atl", "ATLANTIC"), ("pac", "PACIFIC")):
        B = F[F.basin == b]
        L.append(f"== {name} ==")
        L += fmt(summary(B))
        if b == "atl":
            T = B[~B.terrain]
            s = summary(T)
            L.append("  -- without terrain-flagged fixes --")
            L += [l for l in fmt(s) if any(x in l for x in ("fixes", "RMG", "HF-eq area 1e3", "south of centre", "right of motion"))]
        for st in ("deepening", "mature", "filling"):
            L.append(f"  -- stage {st} --")
            L += fmt(summary(B[B.stage == st]))
        for ep in ("first", "interior", "last"):
            L.append(f"  -- episode position {ep} (multi-fix episodes) --")
            L += [l for l in fmt(summary(B[B.episode == ep])) if any(x in l for x in ("fixes", "RMG", "HF-eq area 1e3 km2, median", "max gust kt"))]
        L.append("")
    if (F.phase != "").any():
        L.append("== HART PHASE (00/12 UTC fixes with env_2004 match, both basins) ==")
        for ph in ("asym-cold", "sym-cold", "sym-warm", "asym-warm"):
            L.append(f"  -- {ph} --")
            L += [l for l in fmt(summary(F[F.phase == ph])) if any(x in l for x in
                  ("fixes", "RMG", "HF-eq area 1e3 km2, median", "core", "right of motion", "south of centre", "rear-right"))]
        L.append("")
    W = F[F.wmax.notna()]
    L.append("== SUSTAINED 10 m WIND (WeatherBench2 0.25 deg, one random HF fix per storm, before 2023-01-10) ==")
    for b in ("atl", "pac"):
        d = W[W.basin == b]
        L.append(f"  {b}: fixes {len(d)}; ERA5 max sustained wind median {d.wmax.median():.1f} kt "
                 f"(IQR {d.wmax.quantile(.25):.1f}-{d.wmax.quantile(.75):.1f}); "
                 f"gust factor gmax/wmax median {(d.gmax / d.wmax).median():.2f}")
        L.append(f"     fixes with any ERA5 sustained wind >= 64 kt: {(d.a_w64 > 0).mean() * 100:.1f}%;"
                 f" >= 51.2 kt (64/1.25): {(d.a_w51 > 0).mean() * 100:.1f}%;"
                 f" >= 48 kt: {(d.a_w48 > 0).mean() * 100:.1f}%")
        L.append(f"     area >= 51.2 kt median {d.a_w51.median() / 1e3:.1f}e3 km2; >= 34 kt median "
                 f"{d.a_w34.median() / 1e3:.0f}e3 km2; RMW median {d.wmax_r.median():.0f} km; "
                 f"RMW within 100 km of RMG at {((d.wmax_r - d.gmax_r).abs() <= 100).mean() * 100:.0f}%")
    L.append("")
    L.append("== QuikSCAT comparison (Von Ahn et al. 2006, Fig. 11: 11 Pacific + 6 Atlantic, earth-relative) ==")
    L.append("  Published: HF winds in a crescent south of the centre. No radius, area or quadrant numbers given.")
    for b in ("atl", "pac"):
        hfa = F[(F.basin == b) & (F.a_hf > 0)]
        s = 100 * (hfa.e_se.sum() + hfa.e_sw.sum()) / hfa.a_hf.sum()
        L.append(f"  {b}: share of HF-equivalent area south of the centre {s:.1f}%; "
                 f"median distance of HF-equivalent points {F[F.basin == b].hf_r50.median():.0f} km "
                 f"({F[F.basin == b].hf_r50.median() / NM:.0f} n mi)")
    open(os.path.join(OUT, "structure.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))

    cols = ["time", "track", "basin", "season", "lat", "lon", "msl", "g800", "g800_cat", "heading", "speed_kt",
            "gmax", "gmax_r", "gmax_rel", "gmax_brg", "gmax_lat", "gmax_lon", "gmax_coast_km", "terrain",
            "a_hf", "a_g64", "a_g48", "a_g34", "a_p90", "a_p80"] + QS + ES + \
           ["hf_r50", "hf_rmax", "g48_rmax", "hf_coast50", "own_ocean_frac", "wmax", "wmax_r", "wmax_rel",
            "ws_at_gmax", "a_w64", "a_w51", "a_w48", "a_w34", "h_pmin", "h_peak", "h_first", "h_last",
            "stage", "episode", "dp12", "phase"]
    F[[c for c in cols if c in F]].to_csv(os.path.join(OUT, "fixes.csv"), index=False)
    res = composites(F)
    np.savez_compressed(os.path.join(OUT, "composite.npz"),
                        **{f"{b}_{st}_{k}": v.astype(np.float32) for (b, st), r in res.items()
                           for k, v in r.items() if k != "n"})
    figure(res, "m", os.path.join(OUT, "composite_motion.png"))
    figure(res, "n", os.path.join(OUT, "composite_northup.png"))


if __name__ == "__main__":
    main()
