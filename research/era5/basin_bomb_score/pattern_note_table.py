"""Table behind the weekly 'pattern loading' note on Outlook live (pipeline A, ERA5 proxy; GEFS-fed z+u index, seasons 2004-05 to 2018-19, Oct-Apr weeks).

Low / typical / high = bottom fifth / middle three fifths / top fifth of the 450 pooled weeks by the GEFS-fed index (index_zu.csv, `<basin>_gefs`).
The two cut points are the 20th and 80th percentile of those 450 values (the live note applies them to new index values).
Per class: number of weeks, HF lows per week (peak gust >= 71.7 kt, week of the peak), bombs per week (24 hPa/24 h scaled by sin(lat)/sin(60), week of first bomb fix),
share of weeks with at least one HF low / bomb. Descriptive: the weights behind the index were fitted on seasons 2004-05 to 2014-15 (11 of the 15), so these are
mostly in-sample; PR 277 gives the clean 2015-16 to 2018-19 check. Writes results/pattern_note_table.csv and results/pattern_note_cuts.json. Committed files only."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bomblib as B
import datetime as dt

def hf_weekly(basin, df):
    tr = pd.read_csv(os.path.join(B.ERA, "hf_history/results/all_tracks.csv.gz")); tr = tr[(tr.season >= 2004) & (tr.basin == basin) & (tr.gust800_kt >= 71.7)]
    cnt = {}
    for d in pd.to_datetime(tr.peak_time.astype(str), format="%Y%m%d%H").dt.date:
        s = d.year if d.month >= 10 else d.year - 1; k = (d - dt.date(s, 10, 1)).days // 7; cnt[(s, k)] = cnt.get((s, k), 0) + 1
    return np.array([cnt.get((s, k), 0) for s, k in zip(df.season, df.week)], float)

rows, cuts = [], {}
for b in ("atl", "pac"):
    base = B.weekly(b); gi = pd.read_csv(os.path.join(B.ERA, "gfs_test/results/index_zu.csv"))
    df = base.drop(columns=["idx"]).merge(gi[["season", "week", f"{b}_gefs"]], on=["season", "week"], how="inner").reset_index(drop=True)
    assert len(df) == 450
    x = df[f"{b}_gefs"].values; lo, hi = np.percentile(x, [20, 80]); cuts[b] = dict(low_below=float(lo), high_from=float(hi))
    cls = np.where(x < lo, "low", np.where(x >= hi, "high", "typical"))
    hf = hf_weekly(b, df); bomb = df.bomb.values.astype(float)
    for c in ("low", "typical", "high"):
        m = cls == c
        rows.append(dict(basin=b, cls=c, weeks=int(m.sum()), hf_per_week=hf[m].mean(), bomb_per_week=bomb[m].mean(),
                         share_weeks_hf=(hf[m] > 0).mean(), share_weeks_bomb=(bomb[m] > 0).mean()))
T = pd.DataFrame(rows)
T.round(3).to_csv(os.path.join(B.HERE, "results/pattern_note_table.csv"), index=False)
json.dump(cuts, open(os.path.join(B.HERE, "results/pattern_note_cuts.json"), "w"), indent=1)
print(T.round(3).to_string()); print(cuts)
