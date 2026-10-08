"""Tests T1-T3 of PREREGISTRATION.md, the pre-2004 candidate list and the strongest-storm tables (ERA5 PROXY, pipeline A).

usage: python3 -I analyse.py      (run score.py first; reads committed files only; writes results/)
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
HFH = os.path.join(HERE, "..", "hf_history", "results")
GD = os.path.join(HERE, "..", "gust_depth", "results")
TI = os.path.join(HERE, "..", "tele_intensity", "results")
HF = 71.7
CAND = 0.5
rng = np.random.default_rng(20261008)
NB = 1000
out = []
w = out.append


def track_p(fx, col):
    """Largest P over a track's fixes that are not already HF (hf_now false); NaN if none."""
    g = fx[~fx.hf_now.astype(bool)].groupby("track")[col].max()
    return g


def auc(score, ob):
    ob = ob.astype(bool)
    if ob.sum() == 0 or (~ob).sum() == 0:
        return np.nan
    return mannwhitneyu(score[ob], score[~ob]).statistic / (ob.sum() * (~ob).sum())


def season_boot(df, stat, nb=NB):
    """90% interval of stat(df) resampling whole seasons."""
    seasons = df.season.unique()
    groups = {s: df[df.season == s] for s in seasons}
    vals = []
    for _ in range(nb):
        pick = rng.choice(seasons, len(seasons))
        vals.append(stat(pd.concat([groups[s] for s in pick])))
    return np.nanpercentile(vals, [5, 95])


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    r = p[o] * len(p) / (np.arange(len(p)) + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    return q


def main():
    A = pd.read_csv(os.path.join(HFH, "all_tracks.csv.gz"))
    fx = pd.read_csv(os.path.join(RES, "fix_probs_2004.csv.gz"), dtype={"time": str})
    pre = pd.read_csv(os.path.join(RES, "fix_probs_pre2004.csv.gz"), dtype={"time": str})
    mt = pd.read_csv(os.path.join(GD, "matches.csv"))
    life = pd.read_csv(os.path.join(HFH, "lifecycle_events.csv"))[["track", "tc"]]
    listed = set(mt.track.dropna().astype(int))

    # ------------------------------------------------------------------ 2004+ track table
    T = A[A.season >= 2004].copy()
    for m in ("F", "N", "S"):
        T["Ppre_" + m] = T.track.map(track_p(fx, "P_" + m))
        T["Pall_" + m] = T.track.map(fx.groupby("track")["P_" + m].max())
    T["listed"] = T.track.isin(listed)
    T["event"] = T.gust800_kt >= HF
    T = T.merge(life, on="track", how="left")
    T["tc"] = T.tc.fillna(False).astype(bool)
    T["margin"] = T.gust800_kt - HF
    w(f"2004-05..2025-26 pipeline A tracks: {len(T)}; proxy events {int(T.event.sum())}; archive-listed tracks {int(T.listed.sum())}; "
      f"with a P_pre: {int(T.Ppre_N.notna().sum())}")
    w(f"  listed among proxy events {T[T.event].listed.mean():.3f}; listed among sub-threshold tracks {T[~T.event].listed.mean():.4f} "
      f"({int(T[~T.event].listed.sum())} tracks)")

    # ------------------------------------------------------------------ T1
    w("\n== T1  Does P see archive-listed storms that the gust threshold misses?  (sub-threshold tracks, gust < 71.7 kt, with a P_pre)")
    sub = T[(~T.event) & T.Ppre_N.notna()].copy()
    w(f"  n tracks {len(sub)}, archive-listed {int(sub.listed.sum())} in {sub[sub.listed].season.nunique()} seasons")
    t1 = {}
    for m in ("N", "S", "F"):
        col = "Ppre_" + m
        a = auc(sub[col].values, sub.listed.values)
        lo, hi = season_boot(sub, lambda d: auc(d[col].values, d.listed.values))
        hi_r = sub[sub[col] >= CAND].listed.mean()
        lo_r = sub[sub[col] < CAND].listed.mean()
        n_hi = int((sub[col] >= CAND).sum())
        w(f"  {m}: AUC {a:.3f} [{lo:.3f}, {hi:.3f}]; listed rate P_pre >= {CAND}: {hi_r:.3f} (n {n_hi}, "
          f"{int(sub[sub[col] >= CAND].listed.sum())} listed) vs below: {lo_r:.4f}; ratio {hi_r / lo_r if lo_r > 0 else np.nan:.1f}")
        t1[m] = (a, lo, hi, hi_r, lo_r)
    for thr in (0.3, 0.7):
        s = sub[sub.Ppre_N >= thr]
        w(f"  sensitivity N P_pre >= {thr}: n {len(s)}, listed {int(s.listed.sum())} ({s.listed.mean():.3f})")
    a, lo, hi, hr, lr = t1["N"]
    if lo > 0.5 and lr > 0 and hr >= 2 * lr:
        v1 = "yes"
    elif hi <= 0.55:
        v1 = "no"
    else:
        v1 = "can't tell"
    w(f"  T1 verdict (N, pre-registered rule): {v1}")

    # ------------------------------------------------------------------ T2
    w("\n== T2  Unlisted high-P proxy events (N P_pre >= 0.5, gust >= 71.7 kt)")
    H = T[T.event & (T.Ppre_N >= CAND)].copy()
    H["north60"] = (H.basin == "atl") & (H.peak_lat >= 60)
    H["short"] = H.n_dom < 3
    w(f"  high-P proxy events {len(H)} (atl {int((H.basin=='atl').sum())}, pac {int((H.basin=='pac').sum())}); archive-listed {int(H.listed.sum())} "
      f"({H.listed.mean():.3f}); unlisted {int((~H.listed).sum())}")
    U = H[~H.listed].copy()
    w("  unlisted, removed in order:")
    rem = U.copy()
    for lab, col in (("(a) tropical-cyclone-linked", "tc"), ("(b) Atlantic gust peak at or north of 60N", "north60"),
                     ("(c) fewer than 3 in-domain fixes", "short")):
        mask = rem[col].astype(bool)
        w(f"    {lab}: {int(mask.sum())}")
        rem = rem[~mask]
    R = rem
    w(f"  RESIDUAL unlisted: {len(R)} tracks = {len(R) / len(H):.3f} of the high-P proxy events "
      f"(atl {int((R.basin=='atl').sum())}, pac {int((R.basin=='pac').sum())}); seasons with at least one {R.season.nunique()} of 22")
    per = R.groupby("season").size().reindex(range(2004, 2026), fill_value=0)
    w("  residual per season 2004..2025: " + " ".join(str(int(x)) for x in per.values))
    # same exclusion for the listed group, so the comparison is like for like
    Lg = H[H.listed & ~H.tc & ~H.north60 & ~H.short]
    w(f"  listed events after the same exclusions: {len(Lg)}")
    for thr in (0.3, 0.7):
        Hs = T[T.event & (T.Ppre_N >= thr)]
        Us = Hs[(~Hs.listed) & ~Hs.tc & ~((Hs.basin == 'atl') & (Hs.peak_lat >= 60)) & (Hs.n_dom >= 3)]
        w(f"  sensitivity P_pre >= {thr}: high-P events {len(Hs)}, residual unlisted {len(Us)}")
    Hs = T[T.event & (T.Ppre_S >= CAND)]
    Us = Hs[(~Hs.listed) & ~Hs.tc & ~((Hs.basin == 'atl') & (Hs.peak_lat >= 60)) & (Hs.n_dom >= 3)]
    w(f"  sensitivity model S (P_pre >= 0.5): high-P events {len(Hs)}, residual unlisted {len(Us)}")
    # POST HOC (not in the plan): a looser notion of "listed", any archive HF fix within 400 km of any 00/12 fix at the same time.
    import json
    d = json.load(open(os.path.join(HERE, "..", "..", "..", "docs", "data", "hf-lows.json")))
    arch = {}
    for r in d["lows"]:
        e = dict(zip(d["lowFields"], r))
        for f in e["fixes"]:
            f = dict(zip(d["fixFields"], f))
            if f["cat"] == "HF" and f["lat"] is not None and f["lon"] is not None:
                arch.setdefault(str(f["date"]), []).append((e["basin"], f["lat"], f["lon"] % 360))
    def hav(la1, lo1, la2, lo2):
        la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
        a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
        return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
    near = set()
    for r in fx.itertuples():
        for (b, la, lo) in arch.get(r.time, ()):
            if b == r.basin and hav(la, lo, r.lat, r.lon % 360) <= 400:
                near.add(r.track)
                break
    T["listed_any"] = T.track.isin(near)
    H["listed_any"] = H.track.map(T.set_index("track").listed_any)
    w("\n  POST HOC sensitivity: listed if ANY archive HF fix lies within 400 km of any 00/12 fix of the track at the same time")
    w(f"    proxy events listed: strict (peak-position match, PR 52) {T[T.event].listed.mean():.3f}, loose {T[T.event].listed_any.mean():.3f}; "
      f"high-P proxy events listed: strict {H.listed.mean():.3f}, loose {H.listed_any.mean():.3f}")
    U2 = H[~H.listed_any]
    for lab, col in (("tc", "tc"), ("north60", "north60"), ("short", "short")):
        U2 = U2[~U2[col].astype(bool)]
    w(f"    residual unlisted under the loose rule: {len(U2)} ({len(U2) / len(H):.3f} of high-P events); "
      f"gust margin median {U2.margin.median():.2f} kt vs {Lg.margin.median():.2f} kt for listed")
    both = R[R.track.isin(U2.track)]
    w(f"    unlisted under BOTH rules: {len(both)} (atl {int((both.basin=='atl').sum())}, pac {int((both.basin=='pac').sum())}); "
      f"with margin >= 10 kt: {int((both.margin >= 10).sum())}; margin >= 10 kt and minp <= 950: {int(((both.margin >= 10) & (both.minp <= 950)).sum())}")
    both.sort_values("gust800_kt", ascending=False).to_csv(os.path.join(RES, "residual_unlisted_both_rules_2004.csv"), index=False)
    R.sort_values("gust800_kt", ascending=False).to_csv(os.path.join(RES, "residual_unlisted_2004.csv"), index=False)

    # ------------------------------------------------------------------ T3
    w("\n== T3  Residual vs listed events (Mann-Whitney AUC, stratified-by-season permutation, 2,000 draws)")
    G = pd.concat([R.assign(grp=0), Lg.assign(grp=1)])
    G["nfix"] = G.n_dom
    rows = []
    for meas, lab in (("margin", "gust margin over 71.7 kt"), ("Ppre_N", "P_pre (N)"), ("nfix", "in-domain fixes"),
                      ("minp", "minimum pressure")):
        x = G[meas].values.astype(float)
        grp = G.grp.values

        obs = mannwhitneyu(x[grp == 0], x[grp == 1]).statistic / ((grp == 0).sum() * (grp == 1).sum())
        cnt = 0
        sidx = {s: np.where(G.season.values == s)[0] for s in G.season.unique()}
        for _ in range(2000):
            g2 = grp.copy()
            for ix in sidx.values():
                g2[ix] = rng.permutation(g2[ix])
            if (g2 == 0).sum() == 0:
                continue
            u = mannwhitneyu(x[g2 == 0], x[g2 == 1]).statistic / ((g2 == 0).sum() * (g2 == 1).sum())
            cnt += abs(u - 0.5) >= abs(obs - 0.5)
        p = (cnt + 1) / 2001
        rows.append((lab, np.median(x[grp == 0]), np.median(x[grp == 1]), obs, p))
    q = bh([r[4] for r in rows])
    for (lab, m0, m1, o, p), qq in zip(rows, q):
        w(f"  {lab}: residual median {m0:.2f}, listed median {m1:.2f}, P(residual > listed) {o:.3f}, p {p:.4f}, q {qq:.4f}")
    nsig = int((q < 0.05).sum())
    weaker = sum(1 for (lab, m0, m1, o, p), qq in zip(rows, q) if qq < 0.05 and (
        (lab in ("gust margin over 71.7 kt", "P_pre (N)", "in-domain fixes", "minimum pressure") and
         ((m0 < m1) if lab != "minimum pressure" else (m0 > m1)))))
    if len(R) < 30:
        v3 = "can't tell (fewer than 30 residual tracks)"
    elif nsig == 0:
        v3 = "consistent with missed observations (indistinguishable from listed events)"
    elif weaker >= 2:
        v3 = "looks like a proxy artefact (weaker or shorter)"
    else:
        v3 = "can't tell"
    w(f"  T3 verdict (pre-registered rule): {v3}; {nsig} of 4 measures at q < 0.05, {weaker} in the weaker direction")

    # ------------------------------------------------------------------ headline
    if v1 == "yes" or (len(R) >= 30 and v3.startswith("consistent")):
        head = "yes"
    elif v1 == "no" and (len(R) < 10 or v3.startswith("looks like")):
        head = "no"
    else:
        head = "can't tell"
    w(f"\nHEADLINE (pre-registered combination of T1-T3): {head}")

    # ------------------------------------------------------------------ calibration for pre-2004: listed rate by S P_pre, 2004+ all tracks
    w("\n== Archive-listed rate by S P_pre, all 2004+ tracks with a P_pre (calibration for reading pre-2004 candidates)")
    Z = T[T.Ppre_S.notna()].copy()
    for lo_, hi_ in ((0, .1), (.1, .3), (.3, .5), (.5, .7), (.7, .9), (.9, 1.01)):
        s = Z[(Z.Ppre_S >= lo_) & (Z.Ppre_S < hi_)]
        w(f"  S P_pre {lo_:.1f}-{min(hi_, 1):.1f}: tracks {len(s):6d}, listed {int(s.listed.sum()):5d} ({s.listed.mean():.3f}); proxy events {int(s.event.sum())}")
    Zc = Z[Z.event]
    w(f"  among proxy events only: S P_pre >= 0.5 listed {Zc[Zc.Ppre_S >= 0.5].listed.mean():.3f}, below {Zc[Zc.Ppre_S < 0.5].listed.mean():.3f}")
    # like-for-like calibration: the pre-2004 sample is events plus matched null cases, so use the same kind of tracks after 2004
    cat_all = pd.read_csv(os.path.join(HFH, "era5_hf_catalog.csv"))
    Zn = Z[Z.track.isin(cat_all[cat_all.role == "null_case"].track)]
    Ze = Z[Z.track.isin(cat_all[cat_all.role == "event"].track)]
    w(f"  catalog-type tracks after 2004 (same selection as the pre-2004 sample): null cases with S P_pre >= 0.5: {int((Zn.Ppre_S >= 0.5).sum())}, "
      f"archive-listed {Zn[Zn.Ppre_S >= 0.5].listed.mean():.3f}; events with S P_pre >= 0.5: {int((Ze.Ppre_S >= 0.5).sum())}, listed {Ze[Ze.Ppre_S >= 0.5].listed.mean():.3f}")
    T.to_csv(os.path.join(RES, "tracks_2004_summary.csv.gz"), index=False)

    # ------------------------------------------------------------------ pre-2004 candidates
    w("\n== Pre-2004 catalog tracks (seasons 1979-80..2003-04): S P_pre >= 0.5")
    cat = pd.read_csv(os.path.join(HFH, "era5_hf_catalog.csv"))
    pp = pre[~pre.hf_now.astype(bool)].groupby("track").P_S.max().rename("Ppre_S")
    pa = pre.groupby("track").P_S.max().rename("Pall_S")
    Cc = A[A.season < 2004].merge(cat[["track", "role"]], on="track").merge(life, on="track", how="left")
    Cc = Cc.merge(pp, on="track", how="left").merge(pa, on="track", how="left")
    Cc["tc"] = Cc.tc.fillna(False).astype(bool)
    c1 = Cc[Cc.Ppre_S >= CAND]
    w(f"  catalog tracks {len(Cc)}; with a P_pre {int(Cc.Ppre_S.notna().sum())}; P_pre >= 0.5: {len(c1)} "
      f"(proxy events {int((c1.role=='event').sum())}, null cases {int((c1.role=='null_case').sum())}); tc-linked {int(c1.tc.sum())}")
    w(f"  of the {int((Cc.role=='event').sum())} pre-2004 proxy events, {int(((Cc.role=='event') & (Cc.Ppre_S >= CAND)).sum())} reach P_pre >= 0.5 "
      f"and {int(((Cc.role=='event') & Cc.Ppre_S.isna()).sum())} have no P_pre (first in-domain fix already HF)")
    w(f"  null cases with P_pre >= 0.5 (model false alarms by construction): {int(((c1.role=='null_case')).sum())}")
    # after the 111 GB environment pull (approved by Jason): model N on the same pre-2004 fixes
    ppn = pre[~pre.hf_now.astype(bool)].groupby("track").P_N.max().rename("Ppre_N")
    Cc = Cc.merge(ppn, on="track", how="left")
    c2 = Cc[Cc.Ppre_N >= CAND]
    w(f"  model N (environment, after the pull): P_pre >= 0.5 on {len(c2)} tracks (proxy events {int((c2.role=='event').sum())}, "
      f"null cases {int((c2.role=='null_case').sum())}); tc-linked {int(c2.tc.sum())}; events with no P_pre {int(((Cc.role=='event') & Cc.Ppre_N.isna()).sum())}")
    Zn2 = Z[Z.track.isin(cat_all[cat_all.role == "null_case"].track)]
    Ze2 = Z[Z.track.isin(cat_all[cat_all.role == "event"].track)]
    w(f"  2004+ catalog-type tracks, N P_pre >= 0.5: null cases {int((Zn2.Ppre_N >= 0.5).sum())}, archive-listed {Zn2[Zn2.Ppre_N >= 0.5].listed.mean():.3f}; "
      f"events {int((Ze2.Ppre_N >= 0.5).sum())}, listed {Ze2[Ze2.Ppre_N >= 0.5].listed.mean():.3f}")
    w(f"  expected archive-listed among the {int((c2.role=='null_case').sum())} pre-2004 null cases at N P_pre >= 0.5: "
      f"{(c2.role=='null_case').sum() * Zn2[Zn2.Ppre_N >= 0.5].listed.mean():.0f}; among events {(c2.role=='event').sum() * Ze2[Ze2.Ppre_N >= 0.5].listed.mean():.0f}")
    w(f"  correlation of track-level P_pre, N vs S, pre-2004: {Cc[['Ppre_N','Ppre_S']].corr().iloc[0,1]:.3f}")
    Cc.sort_values("Ppre_S", ascending=False).to_csv(os.path.join(RES, "tracks_pre2004_summary.csv"), index=False)

    # ------------------------------------------------------------------ strongest storms
    w("\n== Strongest storms, ERA5 PROXY pipeline A, seasons 1979-80..2025-26 (75,087 tracks)")
    ev = pd.read_csv(os.path.join(TI, "event_table.csv"))[["track", "clim", "anom_clim", "depth_ring", "tmin", "lat_min", "lon_min"]]
    S = A.merge(ev, on="track", how="left").merge(life, on="track", how="left")
    S["tc"] = S.tc.where(S.tc.notna(), np.nan)
    S["gust_flag"] = np.where(S.season <= 2000, "pre-2001 gust drifts upward", "")
    # Named only where date and position match a well-known storm; inferred by the analyst, not from any data file.
    NAMED = {128894: "Braer storm, 10-11 Jan 1993", 388835: "Storm Dennis, 15 Feb 2020", 415097: "Hurricane Fiona (post-tropical), 24 Sep 2022",
             339400: "Bering Sea bomb (ex-Typhoon Nuri), 8 Nov 2014"}
    rows = []
    for b in ("atl", "pac"):
        sb = S[S.basin == b]
        t1_ = sb.nsmallest(25, "minp").assign(rank_by="minp", rank=range(1, 26))
        t2_ = sb[sb.anom_clim.notna()].nsmallest(25, "anom_clim").assign(rank_by="depth_vs_monthly_clim", rank=range(1, 26))
        n_noev = int(t1_.anom_clim.isna().sum())
        w(f"  {b}: top 25 by minp include {n_noev} track(s) that are not catalog events (no climatology-adjusted depth)")
        w(f"  {b}: top 25 by minp: {int((t1_.season <= 2000).sum())} before 2001-02, {int((t1_.season >= 2004).sum())} from 2004-05; "
          f"tropical-cyclone-linked {int((t1_.tc == True).sum())}; top 25 by depth vs climatology: tc-linked {int((t2_.tc == True).sum())}")
        rows += [t1_, t2_]
    top = pd.concat(rows)
    top["date_minp"] = top.tmin.astype("Int64").astype(str).str[:8]
    top["lon_minp"] = np.where(top.basin == "atl", ((top.lon_min + 180) % 360) - 180, top.lon_min)
    top["name_inferred"] = top.track.map(NAMED)
    top = top[["basin", "rank_by", "rank", "track", "season", "date_minp", "minp", "anom_clim", "depth_ring", "gust800_kt", "gust_flag", "tc",
               "name_inferred", "lat_min", "lon_minp", "peak_time", "peak_lat", "peak_lon", "n_dom"]]
    top.to_csv(os.path.join(RES, "strongest_storms.csv"), index=False)
    open(os.path.join(RES, "analysis.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
