"""Combine stage 1, label cases and classes, draw the control sample, write the power note.

Rules are in PREREGISTRATION.md. Reads $ERA5_WORK/greenland_jets/stage1/*.csv and the
intensity fix table; writes results/stage1_times.csv (all candidate times) and
results/sample_times.csv (cases + seeded control sample), and results/power.txt.
Only counts and the power calculation are produced here, before any ingredient is
related to the outcome.
"""
import os, glob
import numpy as np, pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(os.environ.get("ERA5_WORK", os.path.join(HERE, "work")), "greenland_jets")
RES = os.path.join(HERE, "results")
FIX = os.path.join(HERE, "..", "intensity", "results", "fixes_2004.csv.gz")
THR = 71.7
SEED = 20261008
out = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)


def main():
    os.makedirs(RES, exist_ok=True)
    s1 = pd.concat([pd.read_csv(p, dtype={"time": str}) for p in sorted(glob.glob(f"{WORK}/stage1/*.csv"))])
    s1["time"] = s1.time.astype(np.int64)
    f = pd.read_csv(FIX)
    a = f[f.basin == "atl"].copy()
    a["lon"] = np.where(a.lon > 180, a.lon - 360, a.lon)
    a = a[a.lat.between(55, 67) & a.lon.between(-50, -15)]
    ref = a.sort_values("msl").groupby("time").head(1)          # deepest low in the box at each time
    ref = ref[["time", "track", "lat", "lon", "msl", "speed", "heading", "dp12", "season"]].rename(
        columns={"lat": "ref_lat", "lon": "ref_lon", "msl": "ref_msl", "speed": "ref_speed",
                 "heading": "ref_heading", "dp12": "ref_dp12"})
    d = s1.merge(ref, on="time", how="left")
    d["case"] = d.G_T >= THR
    d["tip"] = d.case & d.max_lat.between(59, 61.5) & d.max_lon.between(310, 325)
    d["barrier"] = d.case & ~d.tip
    d["half"] = np.where(d.season % 2 == 0, "discovery", "heldout")
    d.to_csv(f"{RES}/stage1_times.csv", index=False)

    say("Greenland jets: stage 1 counts (ERA5 proxy; outcome G_T, not pipeline A's index)")
    say(f"candidate times {len(d)}, seasons {d.season.nunique()} ({d.season.min()}-{d.season.max()}), "
        f"bytes streamed {d.bytes.sum() / 1e9:.1f} GB")
    say(f"missing reference low after merge: {int(d.ref_msl.isna().sum())}")
    say(f"G_T >= {THR} kt (cases): {int(d.case.sum())} times; tip {int(d.tip.sum())}, barrier {int(d.barrier.sum())}")
    for h, g in d.groupby("half"):
        say(f"  {h}: times {len(g)}, cases {int(g.case.sum())} (tip {int(g.tip.sum())}, barrier {int(g.barrier.sum())}), "
            f"seasons {g.season.nunique()}")
    # episodes: cases joined when 24 h or less apart
    tt = pd.to_datetime(d.time.astype(str), format="%Y%m%d%H")
    d["_t"] = tt
    ep = 0
    for h, g in d[d.case].sort_values("_t").groupby("half"):
        gap = g._t.diff() > pd.Timedelta(hours=24)
        n = int(gap.sum()) + 1
        say(f"  {h}: case episodes (<=24 h apart joined) {n}")
    say("G_T quantiles over all candidate times (kt): " + ", ".join(
        f"{q}: {d.G_T.quantile(q):.1f}" for q in (0.5, 0.75, 0.9, 0.95, 0.99)))
    say("cases by season: " + " ".join(f"{int(k)}:{int(v)}" for k, v in d[d.case].groupby("season").size().items()))
    say("share of candidate times with an eligible point at all (n_ok>0): "
        f"{(d.n_ok > 0).mean():.3f}")

    # control sample: 3 per case, seeded, within each half, among Y_T = 0
    rng = np.random.default_rng(SEED)
    keep = []
    for h, g in d.groupby("half"):
        cases = g[g.case]
        ctrl = g[~g.case]
        m = min(len(ctrl), 3 * len(cases))
        pick = ctrl.sample(n=m, random_state=int(rng.integers(1 << 31)))
        say(f"  {h}: controls drawn {m} of {len(ctrl)} (sampling fraction {m / len(ctrl):.3f})")
        keep += [cases.assign(w=1.0), pick.assign(w=len(ctrl) / m)]
    smp = pd.concat(keep).drop(columns="_t")
    smp.to_csv(f"{RES}/sample_times.csv", index=False)
    say(f"sample times {len(smp)}; stage 2 pull about {len(smp) * 4.4 / 1000:.1f} GB")

    # power: odds ratio per SD detectable at 80% power in the discovery half, 3 controls per
    # case, normal approximation for a standardised covariate: var(beta) ~ 1 / (n_case*m/(m+1) * var(x)),
    # inflated by the design effect between case times and case episodes.
    g = d[(d.half == "discovery")]
    nc = int(g.case.sum())
    tt2 = g[g.case].sort_values("_t")._t
    ne = int((tt2.diff() > pd.Timedelta(hours=24)).sum()) + 1
    deff = nc / max(ne, 1)
    for label, n in (("case times", nc), ("case episodes", ne)):
        se = 1 / np.sqrt(n * 3 / 4)
        mde = (norm.ppf(0.975) + norm.ppf(0.8)) * se
        say(f"discovery half, {label} {n}: smallest log OR per SD at 80% power {mde:.3f} (OR {np.exp(mde):.2f})")
    say(f"design effect (case times per episode) {deff:.2f}")
    open(f"{RES}/power.txt", "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
