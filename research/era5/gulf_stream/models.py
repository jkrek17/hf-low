"""Wind and wave models of the pre-registration, as functions of a station-time table."""
import numpy as np, pandas as pd
import samples as S
from stats_lib import Boot
BINS = [(-9, 0), (3, 6), (6, 9), (9, 99)]      # reference bin 0-3
def wind_xy(d, stations):
    d = d[d.station.isin(stations)]
    X = [(d.station == s).astype(float).values for s in stations] + [np.log(d.vg.values)]
    for lo, hi in BINS: X.append(((d.dT.values >= lo) & (d.dT.values < hi)).astype(float))
    return d, np.column_stack(X), np.log(d.v10.values), len(stations)
def wind_contrast(k, kind):
    def f(b):
        i = k + 1   # first bin index
        return {"high": 0.5 * (b[i + 2] + b[i + 3]) , "top_minus_mid": b[i + 3] - b[i + 2]}[kind]
    return f
def wave_xy(d, stations, instream_interaction=True):
    d = d[d.station.isin(stations) & d.sector.isin(["opposed", "following"])]
    opp = (d.sector == "opposed").astype(float).values; ins = d.instream.astype(float).values
    X = [(d.station == s).astype(float).values for s in stations]
    X += [(d.vbin == v).astype(float).values for v in ["30-40", "40+"]] + [d.dpd.values]
    if instream_interaction: X += [ins, opp, opp * ins]
    else: X += [opp]
    return d, np.column_stack(X), d.wvht_ft.values, len(stations)
