"""Candidate historic storms, 1979-80 to 2003-04: the strongest ERA5 pipeline A events before the archive is complete.

usage: python3 -I research/era5/climatology_atlas/historic.py      (writes results/historic_storms.csv, historic.txt)

PROXY, not the archive, and nothing here is confirmed. Reads committed files only (hf_history catalog).
Ranking is WITHIN the pre-2004 sample, by two separate measures, because the gust index drifts upward before 2001 at
fixed depth (decision 1: pre-2001 gust-based levels are not comparable with later seasons):
  - minimum ERA5 MSLP along the track (depth; no gust-style drift found, but ERA5 pressure over the open ocean in the
    1980s was constrained by fewer observations, so depths may be too shallow rather than too deep);
  - the pipeline A 800 km gust index (within-era ranking only; compare storms with each other, not with 2004+ values).
Seasons 2001-02 to 2003-04 are also listed, but the archive is incomplete there, so a storm missing from the archive
is not evidence against it.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HFH = os.path.join(HERE, "..", "hf_history", "results")
RES = os.path.join(HERE, "results")


def main():
    C = pd.read_csv(os.path.join(HFH, "era5_hf_catalog.csv"))
    E = C[(C.role == "event")].copy()
    LC = pd.read_csv(os.path.join(HFH, "lifecycle_events.csv"))[["track", "tc"]]
    E = E.merge(LC, on="track", how="left")
    E["lon"] = np.where(E.basin == "atl", ((E.peak_lon + 180) % 360) - 180, E.peak_lon % 360)
    E["date"] = E.peak_time.astype(str).str[:8]
    hist = E[E.season <= 2003].copy()
    late = E[E.season >= 2004]
    out = []
    w = out.append
    w("Candidate historic storms from the ERA5 PROXY (pipeline A), seasons 1979-80 to 2003-04. Nothing here is confirmed.")
    w("Rank within the pre-2004 sample only. Gust ranks are within-era (pre-2001 gust drifts upward at fixed depth).")
    w(f"Events: {len(hist)} (atl {int((hist.basin=='atl').sum())}, pac {int((hist.basin=='pac').sum())}); 1979-2000: {int((hist.season<=2000).sum())}.")
    w("\nMedian ERA5 minimum pressure (hPa) by era, as a depth sanity check (descriptive):")
    for b in ("atl", "pac"):
        w(f"  {b}: 1979-2000 {E[(E.basin==b)&(E.season<=2000)].minp.median():.1f} (n={int(((E.basin==b)&(E.season<=2000)).sum())});"
          f" 2001-2003 {E[(E.basin==b)&(E.season.between(2001,2003))].minp.median():.1f};"
          f" 2004+ {late[late.basin==b].minp.median():.1f}")
    rows = []
    for b in ("atl", "pac"):
        h = hist[hist.basin == b]
        for how, col, asc in (("depth", "minp", True), ("gust", "gust800_kt", False)):
            top = h.sort_values(col, ascending=asc).head(15)
            w(f"\n-- {b} top 15 by {how} (ERA5 min MSLP hPa / gust index kt / date of peak gust / position of peak gust)")
            for _, r in top.iterrows():
                w(f"   {r.date}  minp {r.minp:.0f}  gust {r.gust800_kt:.0f}  ({r.peak_lat:.1f}N, {r.lon:.1f}E)  tc-linked={bool(r.tc)}  track {int(r.track)}")
                rows.append(dict(basin=b, ranked_by=how, track=int(r.track), season=int(r.season), date=r.date, minp=r.minp,
                                 gust800_kt=r.gust800_kt, lat=r.peak_lat, lon=r.lon, tc_linked=bool(r.tc)))
    pd.DataFrame(rows).assign(caveat="ERA5 proxy, not observation; candidate only; gust ranks within-era").to_csv(os.path.join(RES, "historic_storms.csv"), index=False, float_format="%.1f")
    open(os.path.join(RES, "historic.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
