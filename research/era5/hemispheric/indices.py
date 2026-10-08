"""Named-index predictors for a window, same conventions as the additive teleconnection test.

NAO, PNA, AO: mean of the CPC daily index over the window. ONI: the monthly value of the window's
midpoint month. MJO: the two-coordinate longitude projection (MJO1 = amp cos(phase), MJO2 = amp sin(phase))
at the last pentad whose centre is <= midpoint + 2 days (no further lag; the window supplies it).
Source: docs/data/teleconnections.json (CPC files retrieved 2026-10-07).
"""
import json, datetime as dt
import numpy as np


class Indices:
    def __init__(self, path):
        t = json.load(open(path))
        self.d0 = dt.date.fromisoformat(t["start"])
        self.daily = {k: np.array([np.nan if v is None else v for v in t[k]["values"]], float) for k in ("nao", "pna", "ao")}
        oy, om = map(int, t["oni"]["start"].split("-"))
        self.oni0 = (oy, om)
        self.oni = np.array([np.nan if v is None else v for v in t["oni"]["values"]], float)
        m = t["mjo"]
        self.mdates = [dt.datetime.strptime(str(x), "%Y%m%d").date() for x in m["dates"]]
        self.phase = np.array([np.nan if v is None else v for v in m["eofPhase"]], float)
        self.amp = np.array([np.nan if v is None else v for v in m["eofAmplitude"]], float)

    def mjo_row(self, day):
        lo, hi = 0, len(self.mdates) - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if self.mdates[mid] <= day + dt.timedelta(days=2):
                lo = mid
            else:
                hi = mid - 1
        return lo

    def window(self, start, end):
        """dict of indices for days start..end (dates, inclusive)."""
        days = [start + dt.timedelta(days=i) for i in range((end - start).days + 1)]
        ix = [(d - self.d0).days for d in days]
        out = {k: float(np.nanmean(self.daily[k][ix])) for k in ("nao", "pna", "ao")}
        mid = days[len(days) // 2]
        out["oni"] = float(self.oni[(mid.year - self.oni0[0]) * 12 + mid.month - self.oni0[1]])
        i = self.mjo_row(mid)
        out["mjo1"] = float(self.amp[i] * np.cos(self.phase[i]))
        out["mjo2"] = float(self.amp[i] * np.sin(self.phase[i]))
        return out
