"""Duration labels for RA-22 (ERA5 proxy, pipeline A, and the archive). Plan: PREREGISTRATION.md.

Writes results/track_nhf.csv (one row per pipeline A HF-equivalent track: n_hf = in-domain 6-hourly fixes with
g800 >= 71.7 kt, run = longest run of consecutive 6-hourly HF fixes) and, in work/ (ignored), the derived track tables
that the unchanged PR 14 / PR 63 / PR 64 scripts read:

  all_tracks_k{k}.csv.gz        HF_k: a track keeps its gust800_kt only if n_hf >= k (otherwise 0.0): k = 1, 2, 3
  all_tracks_run{k}.csv.gz      HF_k by longest consecutive run >= k (secondary S1)
  all_tracks_gcut{k}.csv.gz     count-matched peak-gust cut (secondary S3): the N_k highest-gust tracks per basin over
                                2004-05..2025-26 keep a gust of 99.0 kt, the rest 0.0, where N_k = tracks with HF_k in
                                that basin and period (a count, not an outcome)

Demoted tracks stay in the file, so they remain cyclones in every denominator; only the HF flag changes.
usage: labels.py WORKDIR
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ERA = os.path.abspath(os.path.join(HERE, ".."))
THR = 71.7
FIRST, LAST = 2004, 2025


def build():
    cat = pd.read_csv(os.path.join(ERA, "hf_history/results/era5_hf_catalog_tracks.csv"), dtype={"time": str})
    ev = pd.read_csv(os.path.join(ERA, "hf_history/results/lifecycle_events.csv"))
    cat["t"] = pd.to_datetime(cat.time, format="%Y%m%d%H")
    cat = cat[cat.track.isin(set(ev.track))].sort_values(["track", "t"])
    cat["indom"] = cat.basin.fillna("") != ""
    cat["hf"] = cat.indom & (cat.g800 >= THR)
    rows = []
    for tr, g in cat.groupby("track", sort=False):
        hf = g.hf.values
        t = g.t.values
        best = cur = 0
        for j in range(len(hf)):
            if hf[j]:
                cur = cur + 1 if (j > 0 and hf[j - 1] and (t[j] - t[j - 1]) == np.timedelta64(6, "h")) else 1
                best = max(best, cur)
            else:
                cur = 0
        rows.append(dict(track=tr, n_hf=int(hf.sum()), run=best))
    out = pd.DataFrame(rows)
    chk = out.merge(ev[["track", "n_hf"]], on="track", suffixes=("", "_lc"))
    assert len(out) == len(ev) == 4157, (len(out), len(ev))
    assert (chk.n_hf == chk.n_hf_lc).all(), "n_hf disagrees with lifecycle_events.csv"
    return out


def derive(work, nh):
    T = pd.read_csv(os.path.join(ERA, "hf_history/results/all_tracks.csv.gz"), dtype={"start": str, "peak_time": str})
    T = T.merge(nh, on="track", how="left")
    hf1 = T.gust800_kt >= THR
    assert hf1.sum() == len(nh) and T.loc[hf1, "n_hf"].notna().all()
    base = T.drop(columns=["n_hf", "run"])
    for k in (1, 2, 3):
        for name, col in (("k", "n_hf"), ("run", "run")):
            if name == "run" and k == 1:
                continue
            x = base.copy()
            x.loc[hf1 & ~(T[col] >= k), "gust800_kt"] = 0.0
            x.to_csv(os.path.join(work, f"all_tracks_{name}{k}.csv.gz"), index=False)
        # count-matched peak-gust cut
        x = base.copy()
        x["gust800_kt"] = 0.0
        for b in ("atl", "pac"):
            sel = (T.basin == b) & T.season.between(FIRST, LAST)
            n_k = int((hf1 & sel & (T.n_hf >= k)).sum())
            top = T[sel].sort_values("gust800_kt", ascending=False).index[:n_k]
            x.loc[top, "gust800_kt"] = 99.0
            # tracks outside 2004-2025 keep their original flag at k = 1 (they are never analysed)
        x.to_csv(os.path.join(work, f"all_tracks_gcut{k}.csv.gz"), index=False)


if __name__ == "__main__":
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    nh = build()
    nh.to_csv(os.path.join(HERE, "results/track_nhf.csv"), index=False)
    derive(work, nh)
    print(nh.n_hf.value_counts().sort_index().head(6).to_dict(), "longest-run mean", nh.run.mean())
