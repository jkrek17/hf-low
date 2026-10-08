"""Life cycle of pipeline A's HF-equivalent events: how long from genesis to
hurricane-force-equivalent winds, how long they last, and where in the storm's
life they occur. ERA5 proxy, not the archive.

usage: lifecycle.py            (reads results/, writes results/lifecycle.txt
                                and results/lifecycle_events.csv)

Reads only committed files: results/era5_hf_catalog.csv,
results/era5_hf_catalog_tracks.csv, ../tc_candidates.csv and the archive CSVs
under data/hf_lows/. No ERA5 access.

Definitions (each is a method choice; see the header of lifecycle.txt):
- HF fix: an in-domain 6-hourly fix with g800 >= 71.7 kt, the same rule that
  makes a track an event. Out-of-domain fixes do not count, as in apply.py.
- Genesis: the track's first fix (first 6-hourly MSLP minimum below 1010 hPa
  anywhere in 20-75N). It is "observed" only when that fix is >= 1000 hPa and
  north of 21N; otherwise the low was already mature or came in from the
  tropics when first detected, and lead times are left-censored.
- Onset: the first HF fix. Censored at onset when the first in-domain fix is
  already HF, or when an out-of-domain fix >= 71.7 kt precedes it.
- Deepening: 24 h MSLP fall in Bergerons, (p(t-24h) - p(t)) / 24 h scaled by
  sin 60 / sin(lat at mid-window); its time is the middle of the window.
- Episodes: runs of consecutive HF fixes; one non-HF fix ends a run.
- Tropical-cyclone link: any fix within 400 km of an IBTrACS point in
  ../tc_candidates.csv at the same time and basin. Those events are reported
  apart and left out of the main numbers, since their genesis is tropical.
- Main sample: seasons 2004-05 to 2025-26 (RECORD_START). Gust values drift
  upward before 2001, which moves threshold crossings earlier and lengthens
  durations, so 1979-80 to 2000-01 is shown separately and labelled.
- Uncertainty: 95% intervals from a season-block bootstrap (seasons resampled
  with replacement within basin, 2,000 draws, seed 1).
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
THR = 71.7
RE = 6371.0
RECORD_START = 2004
NBOOT = 2000


def gc(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * RE * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def hrs(a, b):
    return (a - b) / np.timedelta64(1, "h")


def load():
    C = pd.read_csv(os.path.join(RES, "era5_hf_catalog.csv"))
    E = C[C.role == "event"].set_index("track")
    T = pd.read_csv(os.path.join(RES, "era5_hf_catalog_tracks.csv"))
    T = T[T.track.isin(E.index)].copy()
    T["t"] = pd.to_datetime(T.time.astype(str), format="%Y%m%d%H")
    T["indom"] = T.basin.fillna("") != ""
    T = T.sort_values(["track", "t"]).reset_index(drop=True)
    return E, T


def tc_link(E, T):
    tc = pd.read_csv(os.path.join(HERE, "..", "tc_candidates.csv"))
    tc["t"] = pd.to_datetime(tc.time_utc).dt.tz_localize(None)
    tc["lon"] = tc.lon % 360
    m = T[["track", "t", "lat", "lon"]].merge(tc[["t", "basin", "lat", "lon"]], on="t", suffixes=("", "_tc"))
    m = m[m.basin == E.basin.reindex(m.track).values]
    near = gc(m.lat, m.lon, m.lat_tc, m.lon_tc) < 400
    return set(m.track[near.values])


def per_track(g):
    g = g.reset_index(drop=True)
    t, p, la, lo = g.t.values, g.msl.values, g.lat.values, g.lon.values
    hf = (g.indom & (g.g800 >= THR)).values
    i_on = int(np.argmax(hf))
    first_dom = int(np.argmax(g.indom.values))
    pre_out = bool(((~g.indom) & (g.g800 >= THR)).values[:i_on].any())
    # 24 h deepening, Bergerons, needs fixes 4 steps apart and contiguous in time
    berg = np.full(len(g), np.nan)
    for k in range(4, len(g)):
        if hrs(t[k], t[k - 4]) == 24:
            mid = la[k - 2]
            berg[k - 2] = (p[k - 4] - p[k]) / 24.0 * np.sin(np.radians(60)) / np.sin(np.radians(mid))
    i_md = int(np.nanargmax(berg)) if np.isfinite(berg).any() else -1
    i_mp = int(np.argmin(p))
    hf_dom = np.where(g.indom.values, g.g800.values, -1)
    i_pg = int(np.argmax(hf_dom))
    # episodes
    runs, k = [], 0
    while k < len(hf):
        if hf[k]:
            j = k
            while j + 1 < len(hf) and hf[j + 1] and hrs(t[j + 1], t[j]) == 6:
                j += 1
            runs.append((k, j))
            k = j + 1
        else:
            k += 1
    i_off = runs[-1][1]
    dist = sum(gc(la[a:b], lo[a:b], la[a + 1:b + 1], lo[a + 1:b + 1]).sum() for a, b in runs if b > a)
    nstep = sum(b - a for a, b in runs)
    return dict(
        gen_obs=bool(p[0] >= 1000 and la[0] > 21),
        p_gen=p[0],
        cens_dom=bool(i_on == first_dom),
        cens_out=pre_out,
        cens_gen=bool(i_on == 0),
        h_gen_on=hrs(t[i_on], t[0]),
        n_hf=int(hf.sum()),
        h_hf=6.0 * int(hf.sum()),
        h_span=hrs(t[i_off], t[i_on]),
        n_ep=len(runs),
        h_on_minp=hrs(t[i_on], t[i_mp]),
        h_pg_minp=hrs(t[i_pg], t[i_mp]),
        h_on_maxdeep=hrs(t[i_on], t[i_md]) if i_md >= 0 else np.nan,
        maxdeep=np.nanmax(berg) if i_md >= 0 else np.nan,
        frac_before_minp=float((hf & (np.arange(len(g)) < i_mp)).sum() / hf.sum()),
        frac_after_minp=float((hf & (np.arange(len(g)) > i_mp)).sum() / hf.sum()),
        lat_on=la[i_on], lat_minp=la[i_mp], lat_pg=la[i_pg],
        minp=p[i_mp], p_on=p[i_on],
        spd_hf=(dist / (6.0 * nstep) / 1.852) if nstep else np.nan,   # kt
        dist_hf=dist,
        t_on=t[i_on], t_off=t[i_off],
        mon_on=pd.Timestamp(t[i_on]).month,
    )


def archive():
    rows = []
    for b, f in (("atl", "HF_Data_-_Atl.csv"), ("pac", "HF_Data_-_Pac.csv")):
        x = pd.read_csv(os.path.join(ROOT, "data", "hf_lows", f), dtype={"ID": str, "date": str})
        x["t"] = pd.to_datetime(x.date, format="%Y%m%d%H", errors="coerce")
        x = x.dropna(subset=["t"])
        x["cat"] = x.Category.astype(str).str.upper().str.strip()
        x["key"] = b + ":" + x.ID.str.strip()
        rows.append(x[["key", "t", "cat"]])
    A = pd.concat(rows)
    h = A[A.cat == "HF"].groupby("key").t
    first = A.groupby("key").t.min()
    return pd.DataFrame({"a_on": h.min(), "a_off": h.max(), "a_nhf": h.size(), "a_first": first.reindex(h.min().index)})


def boot(df, stat, seed=1):
    """Season-block bootstrap within basin: returns (estimate, lo, hi)."""
    rng = np.random.default_rng(seed)
    est = stat(df)
    groups = {k: g for k, g in df.groupby(["basin", "season"])}
    by_b = {}
    for (b, s) in groups:
        by_b.setdefault(b, []).append((b, s))
    vals = []
    for _ in range(NBOOT):
        parts = []
        for b, keys in by_b.items():
            pick = rng.integers(0, len(keys), len(keys))
            parts += [groups[keys[i]] for i in pick]
        vals.append(stat(pd.concat(parts)))
    vals = np.array(vals, dtype=float)
    return est, np.nanpercentile(vals, 2.5), np.nanpercentile(vals, 97.5)


def fmt(x, nd=0):
    e, lo, hi = x
    return f"{e:.{nd}f} [{lo:.{nd}f}, {hi:.{nd}f}]"


def main():
    E, T = load()
    tcs = tc_link(E, T)
    L = pd.DataFrame({k: per_track(g) for k, g in T.groupby("track")}).T
    L.index.name = "track"
    L = L.join(E[["basin", "season", "gust800_kt", "archive_events"]])
    L["tc"] = L.index.isin(tcs)
    for c in L.columns:
        if c not in ("basin", "archive_events", "t_on", "t_off"):
            L[c] = pd.to_numeric(L[c])
    L["era"] = np.where(L.season >= RECORD_START, "2004+", np.where(L.season <= 2000, "1979-2000", "2001-03"))
    L["clean_lead"] = L.gen_obs & ~L.cens_dom & ~L.cens_out & ~L.tc

    out = []
    w = out.append
    w("Life cycle of pipeline A HF-equivalent events (ERA5 PROXY, not the archive)")
    w("lifecycle.py; inputs: results/era5_hf_catalog*.csv, ../tc_candidates.csv, data/hf_lows/*.csv")
    w(f"HF fix = in-domain 6-hourly fix with g800 >= {THR} kt. Durations are counts of 6-hourly fixes x 6 h.")
    w("Brackets: 95% season-block bootstrap interval (seasons resampled within basin, 2,000 draws).")
    w("Main sample: seasons 2004-05..2025-26, tropical-cyclone-linked events excluded.")
    w("")
    w("== Sample ==")
    for era in ("2004+", "2001-03", "1979-2000"):
        s = L[L.era == era]
        w(f"{era:10s} events {len(s):5d}  tc-linked {int(s.tc.sum()):4d}  "
          + "  ".join(f"{b} {int((s.basin == b).sum())}" for b in ("atl", "pac")))
    M = L[(L.era == "2004+") & ~L.tc]
    w(f"main sample (2004+, no TC link): {len(M)} events, {M.season.nunique()} seasons")
    w("")
    w("== Censoring (main sample) ==")
    for b in ("atl", "pac", "all"):
        s = M if b == "all" else M[M.basin == b]
        w(f"{b}: genesis observed (first fix >=1000 hPa, >21N) {s.gen_obs.mean():.3f}; "
          f"HF at first in-domain fix {s.cens_dom.mean():.3f}; out-of-domain HF before onset {s.cens_out.mean():.3f}; "
          f"usable for lead time {s.clean_lead.mean():.3f} (n={int(s.clean_lead.sum())})")
    w("")

    def q(col, p):
        return lambda d: np.nanpercentile(d[col].astype(float), p)

    def share(expr):
        return lambda d: float(np.mean(expr(d)))

    w("== Main numbers, 2004-05..2025-26, no TC link ==")
    for b in ("atl", "pac", "all"):
        s = M if b == "all" else M[M.basin == b]
        c = s[s.clean_lead]
        w(f"--- {b} (n={len(s)}; lead-time n={len(c)})")
        w(f"genesis -> first HF fix, h: median {fmt(boot(c, q('h_gen_on', 50)))}; "
          f"quartiles {np.percentile(c.h_gen_on, 25):.0f}/{np.percentile(c.h_gen_on, 75):.0f}; "
          f"share <= 24 h {fmt(boot(c, share(lambda d: d.h_gen_on <= 24)), 2)}")
        w(f"hours at HF: median {fmt(boot(s, q('h_hf', 50)))}; quartiles "
          f"{np.percentile(s.h_hf, 25):.0f}/{np.percentile(s.h_hf, 75):.0f}; 90th pct {np.percentile(s.h_hf, 90):.0f}; "
          f"share 6 h only {fmt(boot(s, share(lambda d: d.h_hf == 6)), 2)}; share >= 24 h {fmt(boot(s, share(lambda d: d.h_hf >= 24)), 2)}")
        w(f"first-to-last HF span, h: median {np.median(s.h_span):.0f}; "
          f"share with >1 HF episode {fmt(boot(s, share(lambda d: d.n_ep > 1)), 2)}")
        w(f"onset minus min-pressure time, h: median {fmt(boot(s, q('h_on_minp', 50)))}; "
          f"onset before min pressure {fmt(boot(s, share(lambda d: d.h_on_minp < 0)), 2)}; at {np.mean(s.h_on_minp == 0):.2f}")
        w(f"peak gust minus min-pressure time, h: median {fmt(boot(s, q('h_pg_minp', 50)))}; "
          f"peak gust before min pressure {fmt(boot(s, share(lambda d: d.h_pg_minp < 0)), 2)}; "
          f"at {np.mean(s.h_pg_minp == 0):.2f}; after {np.mean(s.h_pg_minp > 0):.2f}")
        w(f"onset minus max-deepening mid-time, h: median {fmt(boot(s, q('h_on_maxdeep', 50)))}")
        w(f"share of HF fixes before / after min-pressure fix: {s.frac_before_minp.mean():.2f} / {s.frac_after_minp.mean():.2f}")
        w(f"max 24 h deepening, Bergerons: median {np.nanmedian(s.maxdeep):.2f}; share >= 1 (bomb) {np.mean(s.maxdeep >= 1):.2f}")
        w(f"latitude at onset / min pressure / peak gust: median {np.median(s.lat_on):.1f} / "
          f"{np.median(s.lat_minp):.1f} / {np.median(s.lat_pg):.1f}")
        w(f"translation speed while HF, kt: median {np.nanmedian(s.spd_hf):.0f}; "
          f"distance covered while HF (events with >=2 consecutive HF fixes), km: median "
          f"{np.median(s.dist_hf[s.dist_hf > 0]):.0f}")
    w("")
    w("== Basin difference (Atlantic minus Pacific), main sample ==")
    for col, lab, sub in (("h_gen_on", "genesis->onset median h", "clean"), ("h_hf", "hours at HF median", "all"),
                          ("h_pg_minp", "peak gust - min pressure median h", "all")):
        def st(d, col=col, sub=sub):
            d = d[d.clean_lead] if sub == "clean" else d
            return np.nanmedian(d[d.basin == "atl"][col]) - np.nanmedian(d[d.basin == "pac"][col])
        w(f"{lab}: {fmt(boot(M, st))}")
    w("")
    w("== By month of onset, main sample, both basins ==")
    M2 = M.assign(part=M.mon_on.map(lambda m: "Sep-Nov" if m in (9, 10, 11) else "Dec-Feb" if m in (12, 1, 2)
                                    else "Mar-May" if m in (3, 4, 5) else "Jun-Aug"))
    for part in ("Sep-Nov", "Dec-Feb", "Mar-May", "Jun-Aug"):
        s = M2[M2.part == part]
        c = s[s.clean_lead]
        if len(s) < 10:
            w(f"{part}: n={len(s)} (too few)")
            continue
        w(f"{part}: n={len(s)}; genesis->onset median {fmt(boot(c, q('h_gen_on', 50)))} h (n={len(c)}); "
          f"hours at HF median {np.median(s.h_hf):.0f}; share >=24 h {np.mean(s.h_hf >= 24):.2f}; "
          f"peak gust - min p median {np.median(s.h_pg_minp):.0f} h")
    w("")
    w("== Tropical-cyclone-linked events, 2004+ (genesis is tropical; lead time not meaningful) ==")
    s = L[(L.era == "2004+") & L.tc]
    w(f"n={len(s)}; hours at HF median {np.median(s.h_hf):.0f}; onset latitude median {np.median(s.lat_on):.1f}; "
      f"onset months Aug-Nov share {np.mean(s.mon_on.isin([8, 9, 10, 11])):.2f}")
    w("")
    w("== Earlier seasons, 1979-80..2000-01, no TC link (PROXY; gust index drifts upward before 2001,")
    w("   which moves crossings earlier and lengthens durations; not comparable with 2004+) ==")
    s = L[(L.era == "1979-2000") & ~L.tc]
    c = s[s.clean_lead]
    w(f"n={len(s)}; genesis->onset median {np.median(c.h_gen_on):.0f} h (n={len(c)}); hours at HF median "
      f"{np.median(s.h_hf):.0f}; share >=24 h {np.mean(s.h_hf >= 24):.2f}; peak gust - min p median "
      f"{np.median(s.h_pg_minp):.0f} h")
    w("")
    w("== Check against the archive (events matched to one archive event, 2004+, no TC link) ==")
    A = archive()
    one = M[M.archive_events.fillna("").str.count(";").eq(0) & M.archive_events.notna()].copy()
    one = one.join(A, on="archive_events", how="inner")
    d_on = hrs(one.t_on.values.astype("datetime64[ns]"), one.a_on.values.astype("datetime64[ns]"))
    a_hours = 6.0 * one.a_nhf
    w(f"matched events n={len(one)}")
    w(f"ERA5 onset minus archive first HF fix, h: median {np.median(d_on):.0f}; quartiles "
      f"{np.percentile(d_on, 25):.0f}/{np.percentile(d_on, 75):.0f}; within +/-12 h {np.mean(np.abs(d_on) <= 12):.2f}")
    w(f"hours at HF, ERA5 vs archive (HF fixes x 6 h): medians {np.median(one.h_hf):.0f} vs {np.median(a_hours):.0f}; "
      f"Spearman rho {pd.Series(one.h_hf.values).rank().corr(pd.Series(a_hours.values).rank()):.2f}")
    a_lead = hrs(one.a_on.values.astype("datetime64[ns]"), one.a_first.values.astype("datetime64[ns]"))
    w(f"archive: first fix already HF {np.mean(a_lead == 0):.2f}; archive first fix minus ERA5 genesis, h: median "
      f"{np.median(hrs(one.a_first.values.astype('datetime64[ns]'), one.t_on.values.astype('datetime64[ns]') - one.h_gen_on.values.astype(float) * np.timedelta64(1, 'h'))):.0f}")

    txt = "\n".join(out) + "\n"
    open(os.path.join(RES, "lifecycle.txt"), "w").write(txt)
    cols = ["basin", "season", "era", "tc", "gen_obs", "cens_dom", "cens_out", "clean_lead", "p_gen", "h_gen_on",
            "n_hf", "h_hf", "h_span", "n_ep", "h_on_minp", "h_pg_minp", "h_on_maxdeep", "maxdeep", "lat_on",
            "lat_minp", "lat_pg", "minp", "spd_hf", "mon_on"]
    L[cols].round(2).to_csv(os.path.join(RES, "lifecycle_events.csv"))
    print(txt)


if __name__ == "__main__":
    main()
