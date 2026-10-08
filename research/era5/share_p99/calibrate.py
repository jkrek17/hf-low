"""RA-28 calibration: the 00/12 UTC maximum, 98th and 99th percentile gust-index cuts against the archive, once, before any outcome.

Same recipe as hf_history/calibrate.py (forecast count = archive count on a 0.1 kt grid, 400 km / 800 km matching via hf_history/calib.match),
run on the 00/12 UTC fixes of fixes_2004.csv.gz. No pattern index and no cyclone-level outcome enters. ERA5 PROXY, pipeline A.
usage: python3 -I calibrate.py REPO EXTRACT_DIR
Writes results/cuts.json, results/calibration.txt, results/pull_check.txt, results/fix_p99.csv.gz (pulled percentiles, small).
"""
import glob, json, os, sys
import numpy as np, pandas as pd

root, xdir = sys.argv[1], sys.argv[2]
HERE = os.path.join(root, "research/era5/share_p99")
sys.path.insert(0, os.path.join(root, "research/era5/hf_history"))
import calib  # noqa: E402

FLOOR = 60.0
CAL = calib.CAL
GRID = np.round(np.arange(55, 85, 0.1), 1)
out = []
say = lambda s="": (print(s), out.append(s))

# ---- pulled percentiles and the g800 gate
X = pd.concat([pd.read_csv(f, dtype={"time": str}) for f in glob.glob(os.path.join(xdir, "*.csv"))], ignore_index=True)
X = X.drop_duplicates(["track", "time"])
X.to_csv(os.path.join(HERE, "results/fix_p99.csv.gz"), index=False, float_format="%.3f")
ok = X.g800_re.notna()
d = (X.g800_re - X.g800_cat).abs()
gate = f"fixes pulled {len(X)}; with owned cells {int(ok.sum())}; re-derived g800 within 0.5 kt of catalog: {(d[ok] <= 0.5).mean():.4%} ({int((d[ok] <= 0.5).sum())} of {int(ok.sum())}); max diff {d[ok].max():.2f} kt; matched centre > 0 km: {int((X.match_km > 0).sum())}"
open(os.path.join(HERE, "results/pull_check.txt"), "w").write(gate + "\n")
say("PULL GATE: " + gate)
assert (d[ok] <= 0.5).mean() >= 0.99, "g800 not reproduced"

# ---- fixes table with the three indices
F = pd.read_csv(os.path.join(root, "research/era5/intensity/results/fixes_2004.csv.gz"), dtype={"time": str})
F = F[F.basin.isin(["atl", "pac"])].merge(X[["track", "time", "p98", "p99"]], on=["track", "time"], how="left")
F["g_max"] = F.g800
F["g_p99"] = F.p99.fillna(0.0)
F["g_p98"] = F.p98.fillna(0.0)
miss = F[(F.g800 >= FLOOR) & F.p99.isna() & (F.season.between(2021, 2025))]
say(f"calibration-season fixes with g800 >= {FLOOR:.0f} kt and no pulled percentile: {len(miss)} (must be 0)")
assert len(miss) == 0

A = pd.read_csv(os.path.join(root, "research/era5/hf_history/results/all_tracks.csv.gz"))
A = A[A.basin.isin(["atl", "pac"])].set_index("track")


def track_table(metric):
    g = F.groupby("track")[metric].max().rename("index")
    T = A[["basin", "season"]].join(g, how="left")
    T["index"] = T["index"].fillna(0.0)
    return T


d_ = json.load(open(calib.ARCH))
ev = [dict(zip(d_["lowFields"], r)) for r in d_["lows"]]
for e in ev:
    e["fixes"] = [dict(zip(d_["fixFields"], f)) for f in e["fixes"]]
P = F[F.season.between(2021, 2025)]
M = calib.match(P, ev, set(CAL))
inS = {k: v for k, v in M.items() if v[1]["season"] in CAL}
say(f"Calibration seasons {CAL[0]}-{(CAL[0] + 1) % 100:02d}..{CAL[-1]}-{(CAL[-1] + 1) % 100:02d}, matching on 00/12 UTC fixes: "
    f"{len(inS)} archive events, {sum(t is not None for t, _ in inS.values())} matched to a track")
KEYS = ["n_obs", "n_fc", "bias", "pod", "far", "csi", "hss"]
fmt = lambda r: "  ".join(f"{k}={r[k]:.2f}" if isinstance(r[k], float) else f"{k}={r[k]}" for k in KEYS)
cuts, tabs = {}, {}
for name, metric in (("M", "g_max"), ("P99", "g_p99"), ("P98", "g_p98")):
    T = track_table(metric)
    tabs[name] = T
    f = calib.fit(T, M, CAL, GRID)
    cuts[name] = float(f["thr"])
    say(f"{name}: cut {f['thr']:.1f} kt  " + fmt(f))
    for b in ("atl", "pac"):
        r = calib.score(T[T.basin == b], {k: v for k, v in M.items() if v[1]["basin"] == b}, CAL, f["thr"])
        say(f"   {b}  " + fmt(r))
    if f["thr"] < FLOOR + 1.0:
        raise SystemExit(f"{name} cut {f['thr']} kt is within 1 kt of the pull floor {FLOOR}: stop and extend the pull")
    say("   leave-one-season-out cut: " + ", ".join(f"{s}: {calib.fit(T, M, [z for z in CAL if z != s], GRID)['thr']:.1f}" for s in CAL))
say()
say("Archive agreement at the fitted cuts, M against P99, pooled: season-block bootstrap of the HSS difference (5 seasons, 2000 draws)")
rng = np.random.default_rng(20261011)
def hss_at(T, seasons, thr): return calib.score(T, M, seasons, thr)["hss"]
obs = hss_at(tabs["P99"], CAL, cuts["P99"]) - hss_at(tabs["M"], CAL, cuts["M"])
bs = []
for _ in range(2000):
    pick = list(rng.choice(CAL, len(CAL)))
    # a resampled season list with repeats: score each season separately and pool the 2x2 counts
    cs = {}
    for nm in ("M", "P99"):
        a = b = c = dd = 0
        for s in pick:
            r = calib.score(tabs[nm], M, [s], cuts[nm]); a += r["a"]; b += r["b"]; c += r["c"]; dd += r["d"]
        n = a + b + c + dd; e_ = ((a + b) * (a + c) + (c + dd) * (b + dd)) / n
        cs[nm] = (a + dd - e_) / (n - e_)
    bs.append(cs["P99"] - cs["M"])
lo, hi = np.percentile(bs, [2.5, 97.5])
say(f"HSS(P99) - HSS(M) = {obs:+.3f}, 95% interval [{lo:+.3f}, {hi:+.3f}] (descriptive, not a decision rule)")
# tracks that flip in the calibration seasons
Ts = {n: tabs[n][tabs[n].season.isin(CAL)] for n in tabs}
fl = int(((Ts["M"]["index"] >= cuts["M"]) != (Ts["P99"]["index"] >= cuts["P99"])).sum())
say(f"calibration-season tracks whose label differs between M and P99: {fl} of {int((Ts['M']['index'] >= cuts['M']).sum())} M-positive")
json.dump({"floor_kt": FLOOR, "grid": "55-85 step 0.1", "cuts_kt": cuts}, open(os.path.join(HERE, "results/cuts.json"), "w"), indent=1)
open(os.path.join(HERE, "results/calibration.txt"), "w").write("\n".join(out) + "\n")
